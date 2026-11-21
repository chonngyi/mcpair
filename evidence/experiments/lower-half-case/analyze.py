"""T068: read frozen statistical units; no transport, network, or subprocesses."""
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'experiments/cluster-accounting'))
from common import EDGES, SPECIES

REGIONS = ['4pi', '0-45', '45-90', '90-135', '135-180']
BANDS = {'middle': list(range(1, 8)), **{f'bin{i}': [i] for i in range(9)}}
METRICS = ['n', 'B', 'N_total']
KEYS = ['g4-default', 'g4-incl42', 'fluka-default', 'fluka-coalesce', 'g4-cluster-off']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, allow_nan=False,
                               separators=(',', ':')) + '\n')


def table(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)


def moments(values):
    x = np.asarray(values, float)
    if x.ndim != 2 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError('Invalid statistical units')
    mean = x.mean(axis=0)
    centered = x - mean
    cov = np.einsum('bi,bj->ij', centered, centered) / (len(x) * (len(x) - 1))
    return mean, cov


def ratio(a, b, va, vb, cab=0.0):
    if b == 0 or b - 2 * math.sqrt(vb) <= 0 <= b + 2 * math.sqrt(vb):
        return dict(value=None, mc_se=None, low2se=None, high2se=None,
                    status='zero_denominator' if b == 0 else 'denominator_unresolved')
    q = a / b
    if a == 0 and va == 0:
        return dict(value=q, mc_se=None, low2se=None, high2se=None, status='zero_observed')
    var = va / b**2 + a*a*vb / b**4 - 2*a*cab / b**3
    if var < -1e-14:
        raise ValueError('Negative propagated variance')
    se = math.sqrt(max(0.0, var))
    return dict(value=q, mc_se=se, low2se=q-2*se, high2se=q+2*se, status='computed')


def change(a, b, va, vb, denominator=None):
    """Independent a-b; denominator=(mean,var) is independent of both when given."""
    delta, variance = a-b, va+vb
    se = math.sqrt(variance)
    if denominator is None:
        rel = ratio(a, b, va, vb)
        rel = {k: (v-1 if k in ('value', 'low2se', 'high2se') and v is not None else v)
               for k, v in rel.items()}
    else:
        d, vd = denominator
        rel = ratio(delta, d, variance, vd)
    both_zero = a == 0 and b == 0 and variance == 0
    return dict(delta=delta, delta_mc_se=se, delta_low2se=None if both_zero else delta-2*se,
                delta_high2se=None if both_zero else delta+2*se,
                zero_observed_component=(a == 0 and va == 0) or (b == 0 and vb == 0),
                **{'relative_'+k: v for k, v in rel.items()})


def load_units(root, key, expected_manifest):
    if sha(root/'MANIFEST.json') != expected_manifest:
        raise ValueError('Manifest identity mismatch')
    manifest = json.loads((root/'MANIFEST.json').read_text())
    name = 'task/reduction/result.json' if key == 'g4-cluster-off' else f'task/results/{key}/result.json'
    p = root/name
    if sha(p) != manifest['files'][name]['sha256']:
        raise ValueError('Count file identity mismatch')
    rec = json.loads(p.read_text())
    n, histories = (800, 10000) if key.startswith('g4') else (40, 200000)
    rows = rec['units']
    if (rec['key'] != key or rec['status'] != 'passed' or rec['histories_used'] != 8000000
            or len(rows) != n or [r['index'] for r in rows] != list(range(n))
            or any(r['histories'] != histories for r in rows)):
        raise ValueError('Statistical units mismatch')
    integers = key.startswith('g4')
    if rec['representation'] != ('integer_counts' if integers else 'integral_yield_per_primary'):
        raise ValueError('Unexpected representation: '+rec['representation'])
    x = np.asarray([r['counts' if integers else 'values'] for r in rows], float)
    if x.shape != (n, 37, 4, 9) or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('Invalid count tensor')
    if integers:
        if not np.equal(x, np.floor(x)).all():
            raise ValueError('Noninteger G4 counts')
        x /= histories
    return x, dict(key=key, manifest_sha256=expected_manifest, path=name,
                   sha256=sha(p), units=n, histories=histories*n, representation=rec['representation'])


def analyze(old, new, out):
    preservation = json.loads((HERE/'preservation.json').read_text())
    if sha(HERE/'METHODS.md') != preservation['method_sha256_before_analysis']:
        raise ValueError('Analysis method changed after initial pin')
    for p, h in preservation['files'].items():
        if sha(ROOT/p) != h:
            raise ValueError('Protected file changed: '+p)
    original = json.loads((ROOT/'experiments/cluster-kernel/results/ratios.json').read_text())
    expected = list(original['source_manifests'].values())
    # Identify the frozen manifests by their recorded source roles, not input directory names.
    old_hash = next(v for k, v in original['source_manifests'].items() if 't056-' in k)
    new_hash = next(v for k, v in original['source_manifests'].items() if 't060-' in k)
    assert len(expected) == 2
    out.mkdir(exist_ok=False)
    weights = np.zeros((3, len(SPECIES)))
    weights[0, 0] = 1
    for i, s in enumerate(SPECIES[2:], 2):
        weights[1, i] = s['A']-s['Z']
    weights[2] = weights[0]+weights[1]
    stats, yields, sparse, sources = {}, [], [], []
    labels = [f'{region}|bin{b}|{metric}' for region in REGIONS[1:]
              for b in range(9) for metric in METRICS]
    for key in KEYS:
        x, source = load_units(new if key == 'g4-cluster-off' else old, key,
                              new_hash if key == 'g4-cluster-off' else old_hash)
        sources.append(source)
        joint = np.einsum('usab,ms->uabm', x, weights).reshape(len(x), -1)
        jm, jc = moments(joint)
        dump(out/(key+'-covariance.json'), dict(labels=labels, mean=jm.tolist(), covariance=jc.tolist()))
        for ai, region in enumerate(REGIONS):
            regional = x.sum(axis=2) if ai == 0 else x[:, :, ai-1, :]
            for band, indices in BANDS.items():
                species = regional[:, :, indices].sum(axis=2)
                values = np.einsum('us,ms->um', species, weights)
                mu, cov = moments(values)
                stats[key, region, band] = (mu, cov)
                # Independent projection of full joint covariance must reproduce each summary.
                for mi, metric in enumerate(METRICS):
                    ix = [labels.index(f'{r}|bin{b}|{metric}')
                          for r in (REGIONS[1:] if ai == 0 else [region]) for b in indices]
                    assert math.isclose(mu[mi], float(jm[ix].sum()), rel_tol=1e-12, abs_tol=1e-16)
                    assert math.isclose(cov[mi, mi], float(jc[np.ix_(ix, ix)].sum()), rel_tol=1e-11, abs_tol=1e-22)
                    se = math.sqrt(max(0, cov[mi, mi]))
                    zero = mu[mi] == 0
                    yields.append(dict(setting=key, region=region, band=band, metric=metric,
                                       yield_per_primary=float(mu[mi]), empirical_mc_se=se,
                                       low2se=None if zero else float(mu[mi]-2*se),
                                       high2se=None if zero else float(mu[mi]+2*se),
                                       nonzero_units=int(np.count_nonzero(values[:, mi])),
                                       status='zero_observed' if zero else 'computed'))
                if key == 'g4-cluster-off':
                    contributing = species[:, weights[1] > 0].sum(axis=1)
                    count = int(round(float(contributing.sum()) * 10000))
                    bcount = int(round(float(values[:, 1].sum()) * 10000))
                    sparse.append(dict(region=region, band=band, B=float(mu[1]),
                                       empirical_mc_se=float(math.sqrt(cov[1, 1])),
                                       nonzero_batches=int(np.count_nonzero(values[:, 1])),
                                       contributing_ion_crossings=count, bound_neutron_weighted_count=bcount,
                                       poisson_approx_ion_upper95=2.995732273553991/8000000 if count == 0 else None,
                                       poisson_approx_B_upper95=7*2.995732273553991/8000000 if count == 0 else None,
                                       **{'B_over_n_'+k: v for k, v in ratio(mu[1], mu[0], cov[1, 1], cov[0, 0], cov[1, 0]).items()}))
    crosses, responses = [], []
    for region in REGIONS:
        for band in BANDS:
            for g in ['g4-default', 'g4-cluster-off']:
                for f in ['fluka-default', 'fluka-coalesce']:
                    a, ca = stats[g, region, band]
                    b, cb = stats[f, region, band]
                    for mi in [0, 2]:
                        crosses.append(dict(G4=g, FLUKA=f, region=region, band=band, metric=METRICS[mi],
                                            **ratio(a[mi], b[mi], ca[mi, mi], cb[mi, mi])))
            d, cd = stats['g4-default', region, band]
            o, co = stats['g4-cluster-off', region, band]
            i, ci = stats['g4-incl42', region, band]
            for mi, metric in enumerate(METRICS):
                for name, a, b, va, vb, den in [
                    ('cluster-off-minus-default', o[mi], d[mi], co[mi, mi], cd[mi, mi], None),
                    ('incl42-minus-default', i[mi], d[mi], ci[mi, mi], cd[mi, mi], None),
                    ('incl42-minus-cluster-off', i[mi], o[mi], ci[mi, mi], co[mi, mi], (d[mi], cd[mi, mi]))]:
                    responses.append(dict(contrast=name, region=region, band=band, metric=metric,
                                          common_default_yield=float(d[mi]),
                                          shared_absolute_response_covariance=float(cd[mi, mi]),
                                          **change(a, b, va, vb, den)))
    for name, rows in [('yields', yields), ('four-cells', crosses), ('g4-responses', responses), ('cluster-off-B', sparse)]:
        table(out/(name+'.csv'), rows)
    # Reproduce preexisting means/SE and the original reviewed four cells.
    old_yields = list(csv.DictReader((ROOT/'experiments/cluster-accounting/results/yields.csv').open()))
    compared = 0
    for row in old_yields:
        if row['metric'] not in ('n', 'bound_n', 'N_total'):
            continue
        band = {'below': 'bin0', 'above': 'bin8'}.get(row['band'], row['band'])
        mi = {'n': 0, 'bound_n': 1, 'N_total': 2}[row['metric']]
        mu, cov = stats[row['setting'], row['region'], band]
        assert math.isclose(float(row['yield_per_primary']), mu[mi], rel_tol=1e-11, abs_tol=1e-16)
        assert math.isclose(float(row['mc_se']), math.sqrt(cov[mi, mi]), rel_tol=1e-11, abs_tol=1e-16)
        compared += 1
    matched = 0
    for row in original['ratios']:
        if row['numerator'] == 'g4-incl42':
            continue
        v = next(r for r in crosses if r['G4'] == row['numerator'] and r['FLUKA'] == row['denominator']
                 and r['region'] == '4pi' and r['band'] == 'middle' and r['metric'] == 'n')
        assert math.isclose(v['value'], row['ratio'], rel_tol=1e-12)
        assert math.isclose(v['mc_se'], row['mc_se'], rel_tol=1e-12)
        matched += 1
    assert matched == 4
    dump(out/'verification.json', dict(status='passed', method_sha256=sha(HERE/'METHODS.md'),
         script_sha256=sha(__file__), common_sha256=sha(ROOT/'experiments/cluster-accounting/common.py'),
         input_sources=sources, protected_files_unchanged=len(preservation['files']),
         reproduced_t056_yield_and_se_rows=compared, reproduced_t060_four_cells=matched,
         covariance_projection_checks=750, new_histories=0, k8_jobs=0, model_calls=0,
         review_status='pending', output_rows=dict(yields=len(yields),four_cells=len(crosses),
         g4_responses=len(responses),cluster_off_B=len(sparse))))
    print(json.dumps(dict(rows=[len(yields),len(crosses),len(responses),len(sparse)],
                         old_yield_rows=compared, four_cells=matched)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--t056', type=Path, required=True)
    p.add_argument('--t060', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    analyze(args.t056, args.t060, args.out)
