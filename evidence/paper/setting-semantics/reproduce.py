"""Offline reproduction of a T070 snapshot. NumPy only; no native runs or network.

Checks the package inventory before loading data, recomputes all response tables
and joint covariances from statistical units, and scans the retained pair records.
This is author verification, not an independent implementation or kernel replay.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def check_inventory(root):
    manifest = read(root / 'SNAPSHOT.json')
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    require(actual == set(manifest['files']) | {'SNAPSHOT.json'}, 'Package inventory mismatch')
    for name, item in manifest['files'].items():
        path = root / name
        require(not path.is_symlink() and path.resolve().is_relative_to(root.resolve()), 'Unsafe path: ' + name)
        require(path.stat().st_size == item['bytes'] and digest(path) == item['sha256'], 'Hash mismatch: ' + name)
    for role, subset in manifest.get('original_archive_subsets', {}).items():
        prefix = 'data/' + role + '/'
        original = read(root / prefix / 'MANIFEST.json')
        require(manifest['files'][prefix+'MANIFEST.json']['sha256'] == subset['manifest_sha256'],
                'Archive manifest identity: ' + role)
        require(len(original['files']) == subset['original_manifest_files'] and
                len(original['files'])-len(subset['selected_files']) == subset['omitted_files'],
                'Archive subset count: ' + role)
        for name in subset['selected_files']:
            require(manifest['files'][prefix+name]['sha256'] == original['files'][name]['sha256'],
                    'Original archive binding: ' + role + '/' + name)
    return manifest


def compare_csv(generated, expected):
    with generated.open(newline='') as f:
        left = list(csv.DictReader(f))
    with expected.open(newline='') as f:
        right = list(csv.DictReader(f))
    require(len(left) == len(right), 'Row count: ' + str(expected))
    for index, (a, b) in enumerate(zip(left, right), 2):
        require(a.keys() == b.keys(), 'Columns: ' + str(expected))
        for key in a:
            if a[key] == b[key]:
                continue
            try:
                x, y = float(a[key]), float(b[key])
            except ValueError:
                raise ValueError(f'{expected.name}:{index}:{key}: categorical mismatch') from None
            require(math.isfinite(x) and math.isfinite(y) and
                    math.isclose(x, y, rel_tol=1e-11, abs_tol=1e-20),
                    f'{expected.name}:{index}:{key}: {x} != {y}')
    return len(left)


def response_tables(out, *, load_units=None):
    a = module('frozen_analysis', ROOT / 'experiments/lower-half-case/analyze.py')
    if load_units is None:
        original = read(ROOT / 'experiments/cluster-kernel/results/ratios.json')
        hashes = {role: next(v for k, v in original['source_manifests'].items() if role + '-' in k)
                  for role in ('t056', 't060')}

        def load_units(key):
            role = 't060' if key == 'g4-cluster-off' else 't056'
            return a.load_units(ROOT / 'data' / role, key, hashes[role])
    weights = np.zeros((3, len(a.SPECIES)))
    weights[0, 0] = 1
    weights[1] = [0, 0] + [s['A'] - s['Z'] for s in a.SPECIES[2:]]
    weights[2] = weights[0] + weights[1]
    stats, yields, sparse, sources = {}, [], [], []
    labels = [f'{region}|bin{b}|{metric}' for region in a.REGIONS[1:]
              for b in range(9) for metric in a.METRICS]
    covariance_cells = 0
    for key in a.KEYS:
        x, source = load_units(key)
        sources.append(source)
        joint = np.einsum('usab,ms->uabm', x, weights).reshape(len(x), -1)
        jm, jc = a.moments(joint)
        expected = read(ROOT / f'experiments/lower-half-case/results/{key}-covariance.json')
        require(expected['labels'] == labels, 'Covariance labels')
        require(np.allclose(jm, expected['mean'], rtol=1e-11, atol=1e-20), 'Joint means')
        require(np.allclose(jc, expected['covariance'], rtol=1e-10, atol=1e-24), 'Joint covariance')
        covariance_cells += jc.size
        a.dump(out / (key + '-covariance.json'), dict(labels=labels, mean=jm.tolist(), covariance=jc.tolist()))
        for ai, region in enumerate(a.REGIONS):
            regional = x.sum(axis=2) if ai == 0 else x[:, :, ai-1, :]
            for band, indices in a.BANDS.items():
                species = regional[:, :, indices].sum(axis=2)
                values = np.einsum('us,ms->um', species, weights)
                mu, cov = a.moments(values)
                stats[key, region, band] = mu, cov
                for mi, metric in enumerate(a.METRICS):
                    ix = [labels.index(f'{r}|bin{b}|{metric}')
                          for r in (a.REGIONS[1:] if ai == 0 else [region]) for b in indices]
                    require(math.isclose(mu[mi], jm[ix].sum(), rel_tol=1e-12, abs_tol=1e-16), 'Mean projection')
                    require(math.isclose(cov[mi, mi], jc[np.ix_(ix, ix)].sum(), rel_tol=1e-11, abs_tol=1e-22), 'Covariance projection')
                    se = math.sqrt(max(0, cov[mi, mi]))
                    yields.append(dict(setting=key, region=region, band=band, metric=metric,
                        yield_per_primary=float(mu[mi]), empirical_mc_se=se,
                        low2se=None if mu[mi] == 0 else float(mu[mi]-2*se),
                        high2se=None if mu[mi] == 0 else float(mu[mi]+2*se),
                        nonzero_units=int(np.count_nonzero(values[:, mi])),
                        status='zero_observed' if mu[mi] == 0 else 'computed'))
                if key == 'g4-cluster-off':
                    count = int(round(species[:, weights[1] > 0].sum() * 10000))
                    bcount = int(round(values[:, 1].sum() * 10000))
                    sparse.append(dict(region=region, band=band, B=float(mu[1]),
                        empirical_mc_se=float(math.sqrt(cov[1, 1])),
                        nonzero_batches=int(np.count_nonzero(values[:, 1])),
                        contributing_ion_crossings=count, bound_neutron_weighted_count=bcount,
                        poisson_approx_ion_upper95=2.995732273553991/8000000 if count == 0 else None,
                        poisson_approx_B_upper95=7*2.995732273553991/8000000 if count == 0 else None,
                        **{'B_over_n_'+k: v for k, v in a.ratio(mu[1], mu[0], cov[1, 1], cov[0, 0], cov[1, 0]).items()}))
    crosses, responses = [], []
    for region in a.REGIONS:
        for band in a.BANDS:
            for g in ('g4-default', 'g4-cluster-off'):
                for f in ('fluka-default', 'fluka-coalesce'):
                    u, cu = stats[g, region, band]
                    v, cv = stats[f, region, band]
                    for mi in (0, 2):
                        crosses.append(dict(G4=g, FLUKA=f, region=region, band=band, metric=a.METRICS[mi],
                            **a.ratio(u[mi], v[mi], cu[mi, mi], cv[mi, mi])))
            d, cd = stats['g4-default', region, band]
            o, co = stats['g4-cluster-off', region, band]
            i, ci = stats['g4-incl42', region, band]
            for mi, metric in enumerate(a.METRICS):
                for name, u, v, vu, vv, den in (
                    ('cluster-off-minus-default', o[mi], d[mi], co[mi, mi], cd[mi, mi], None),
                    ('incl42-minus-default', i[mi], d[mi], ci[mi, mi], cd[mi, mi], None),
                    ('incl42-minus-cluster-off', i[mi], o[mi], ci[mi, mi], co[mi, mi], (d[mi], cd[mi, mi]))):
                    responses.append(dict(contrast=name, region=region, band=band, metric=metric,
                        common_default_yield=float(d[mi]), shared_absolute_response_covariance=float(cd[mi, mi]),
                        **a.change(u, v, vu, vv, den)))
    counts = {}
    for name, rows in [('yields', yields), ('four-cells', crosses), ('g4-responses', responses), ('cluster-off-B', sparse)]:
        generated = out / (name + '.csv')
        a.table(generated, rows)
        counts[name] = compare_csv(generated, ROOT / f'experiments/lower-half-case/results/{name}.csv')
        if name == 'g4-responses':
            for line, row in enumerate(rows, 2):
                zero = row['zero_observed_component']
                if zero:
                    for field in ('delta_low2se', 'delta_high2se', 'relative_low2se', 'relative_high2se'):
                        row[field] = None
                row['interval_display_status'] = 'withheld_zero_observed_component' if zero else 'original'
                row['original_csv_row'] = line
        shown = out / (name + '-display.csv')
        a.table(shown, rows)
        compare_csv(shown, ROOT / f'experiments/lower-half-case/presentation-20261114/{name}.csv')
    return dict(rows=counts, covariance_cells=covariance_cells, covariance_projection_checks=750,
                sources=sources, reuse='Frozen moments/ratio/change/unit loader; new portable orchestration')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(__debug__, 'Python -O is unsupported')
    require(not args.out.resolve().is_relative_to(ROOT), 'Write outputs outside the immutable snapshot')
    manifest = check_inventory(ROOT)
    calls = read(ROOT / 'generation-metadata.json')['samples']
    expected_slots = {f'{s}-r{r}' for s in ('Fe-M01', 'Fe-C01', 'Pb-M01', 'Pb-C01') for r in (1, 2)}
    require(len(calls) == 8 and {c['selected_fields']['slot'] for c in calls} == expected_slots,
            'Generation metadata slots')
    for call in calls:
        rec = call['selected_fields']
        require(rec['model'] == 'gpt-6.1-sol' and rec['effort'] == 'high' and rec['passed'],
                'Generation model/status')
        prompt = ROOT / 'experiments/four-condition-inputs/consistency/results/formal-v2' / rec['slot'] / 'prompt.txt'
        require(digest(prompt) == rec['prompt_sha256'], 'Generation prompt identity')
    args.out.mkdir(parents=True, exist_ok=False)
    response = response_tables(args.out)
    verifier = module('writing_evidence', HERE / 'verify_evidence.py')
    evidence = verifier.verify(ROOT / 'data/transport', source_commit=manifest['scientific_baseline'],
                               include_review_history=False)
    (args.out / 'evidence.json').write_text(json.dumps(evidence, ensure_ascii=False, indent=2)+'\n')
    report = dict(status='passed', snapshot_manifest_sha256=digest(ROOT / 'SNAPSHOT.json'),
        files_checked=len(manifest['files']), python=platform.python_version(), numpy=np.__version__,
        response=response, initialization_samples=len(evidence['initialization_samples']),
        recorded_generation_calls=len(calls),
        transport_records=sum(g['records'] for g in evidence['transport'].values()),
        withheld_response_intervals=evidence['withheld_rows'], new_histories=0, model_calls=0,
        kernel_replayed=False, all_trajectory_reduction=False, independent_review=False)
    (args.out / 'verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','files_checked','transport_records','initialization_samples')}))


if __name__ == '__main__':
    main()
