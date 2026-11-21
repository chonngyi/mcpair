"""Recompute the preprint tables, formal native example Fe-M01/g000-f00, and excerpts using only package files."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import struct
import sys
import types

import numpy as np

ROOT = Path(__file__).resolve().parent


def load_json(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_units(key):
    filename = f'units/{key}.npz'
    with np.load(ROOT / filename, allow_pickle=False) as archive:
        array = archive['per_primary']
    n = 800 if key.startswith('g4-') else 40
    require(array.shape == (n, 37, 4, 9) and np.isfinite(array).all() and (array >= 0).all(),
            'Invalid unit tensor: ' + key)
    return array, {'file': filename, 'sha256': sha(ROOT / filename), 'units': n}


def load_settings_sources():
    settings_dir = ROOT / 'settings-source'
    if str(settings_dir) not in sys.path:
        sys.path.insert(0, str(settings_dir))
    import usryield
    import score_reduction

    pkg = types.ModuleType('settings_source')
    pkg.__path__ = [str(settings_dir)]
    sys.modules['settings_source'] = pkg
    base_stub = types.ModuleType('settings_source.base_kernel')
    base_stub.sha = sha

    def _close(a, b, rtol=6e-4, atol=1e-10):
        if a is None or b is None:
            return False
        if type(a) is bool or type(b) is bool:
            return type(a) is type(b) and a == b
        if isinstance(a, dict) and isinstance(b, dict):
            return set(a) == set(b) and all(_close(a[k], b[k], rtol, atol) for k in a)
        if type(a) in (int, float) and type(b) in (int, float):
            return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=rtol, abs_tol=atol)
        return type(a) is type(b) and a == b

    base_stub.close = _close
    sys.modules['settings_source.base_kernel'] = base_stub
    spec = importlib.util.spec_from_file_location('settings_source.scoring', settings_dir / 'scoring.py')
    scoring = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scoring)
    return usryield, scoring, score_reduction


def parse_tsv(path):
    return dict(line.split('\t', 1) for line in path.read_text().splitlines() if '\t' in line)


def verify_formal_index_and_files(formal_dir, units_data, prov):
    idx_path = formal_dir / 'kernel-index.json'
    require(sha(idx_path) == units_data['kernel_index_sha256'],
            'kernel-index.json SHA256 mismatch with units/Fe-M01.json')
    index = load_json(idx_path)
    require(len(index) == 800, 'Fe-M01 kernel-index expected 800 pair entries')
    entry0 = index[0]
    require(entry0['g4_batch'] == 0 and entry0['fluka_stream'] == 0
            and entry0['passed'] is True and entry0['rows'] == 22
            and entry0['record'] == 'results/Fe-M01/pairs/p000/kernel.json'
            and entry0['sha256'] == prov['pair_kernel_excerpt']['source_sha256'],
            'kernel-index[0] -> p000/kernel.json binding mismatch')
    require(prov['unit_id'] == 'Fe-M01/g000-f00' and prov['array_index'] == 0
            and prov['pair_kernel_excerpt']['passed']['value'] is True,
            'p000 pair provenance status mismatch')

    rt_files = prov['pair_kernel_excerpt']['runtime_files']['value']
    for rel, rt_key in [
        ('example/formal-g000-f00/fluka/input.inp', 'fluka/input.inp'),
        ('example/formal-g000-f00/fluka/input001.out', 'fluka/input001.out'),
        ('example/formal-g000-f00/fluka/input001_fort.21', 'fluka/input001_fort.21'),
        ('example/formal-g000-f00/geant4/exits.csv', 'geant4/output/exits.csv'),
        ('example/formal-g000-f00/geant4/primaries.csv', 'geant4/output/primaries.csv'),
        ('example/formal-g000-f00/geant4/identity.tsv', 'geant4/output/identity.tsv'),
        ('example/formal-g000-f00/geant4/incl-physics.tsv', 'geant4/output/incl-physics.tsv'),
        ('example/formal-g000-f00/geant4/cuts-energy.tsv', 'geant4/output/cuts-energy.tsv'),
        ('example/formal-g000-f00/geant4/rng_seed.txt', 'geant4/output/rng_seed.txt'),
        ('example/geometry.gdml', 'geant4/geometry.gdml'),
        ('example/beam.mac', 'geant4/beam.mac'),
        ('example/physics.txt', 'geant4/physics.txt'),
        ('example/fluka.inp', 'fluka/original.inp'),
    ]:
        require(sha(ROOT / rel) == rt_files[rt_key],
                f'SHA256 mismatch vs p000 runtime.files[{rt_key}]: {rel}')
    return idx_path


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def verify_formal_fluka(formal_dir, problem, expected_yield, out, usryield, scoring):
    fluka_inp_text = (formal_dir / 'fluka/input.inp').read_text(encoding='ascii')
    require('PRECISIO' in fluka_inp_text and 'COALESCE' in fluka_inp_text
            and '72690000.' in fluka_inp_text and 'START        200000.' in fluka_inp_text,
            'Formal FLUKA u000 input.inp card verification failed')

    fluka_out_text = (formal_dir / 'fluka/input001.out').read_text(encoding='ascii')
    require('FLRM64 INITIALIZED           0           0       90000         726    72690000' in fluka_out_text
            and 'Total number of primaries run:         200000 for a weight of: 2.000000E+05' in fluka_out_text,
            'Formal FLUKA u000 input001.out seed or primary/weight check failed')
    ds = scoring.summary(fluka_out_text)
    echoes, controls = scoring.echoed_yields(fluka_out_text)
    require(len(ds) == 7 and not controls and set(echoes) == {d['name'] for d in ds}
            and all(d['particle'] == 8 and d['unit'] == 21 and (d['from_region'], d['to_region']) == (1, 2)
                    and d['normalization'] == 1.0 and d['cross_section_mb'] == 1.0 for d in ds)
            and all(int(echoes[d['name']]['second'][5]) % 100 == 3 for d in ds)
            and '124.0' in fluka_inp_text,
            'Formal FLUKA u000 initialized USRYIELD summary mismatch')

    parsed = usryield.read((formal_dir / 'fluka/input001_fort.21').read_text(encoding='ascii'))
    require(parsed['primaries'] == 200000 and parsed['weight'] == 200000.0 and len(parsed['detectors']) == 7,
            'Formal FLUKA fort.21 header or detector count mismatch')
    td_by_name = {td['name']: td for td in parsed['detectors']}
    for d in ds:
        td = td_by_name.get(d['name'])
        require(td is not None and scoring.equivalent(d, {k: v for k, v in td.items() if k != 'values'} | {'unit': 21}),
                'FLUKA input001.out summary vs fort.21 tally mismatch: ' + d['name'])

    edges = problem['scoring']['energy_edges_MeV']
    angles = problem['scoring']['angle_edges_deg']
    expected_cells = [[a, b, c, d] for a, b in zip(edges, edges[1:]) for c, d in zip(angles, angles[1:])]
    keys = {scoring.native_cell_key(c): c for c in expected_cells}
    fluka_factor = scoring.expected(problem, 'fluka')['S07']['factor']
    require(fluka_factor == 0.001, 'Unexpected FLUKA GeV->MeV factor')

    all_grid_cells, cell_rows, contribs = [], [], []
    for detector in parsed['detectors']:
        grid = scoring.variable_cells(detector)
        require(len(grid) == len(detector['values']) == 4, 'Unexpected detector grid length: ' + detector['name'])
        for angle_idx, (raw_cell, native_val) in enumerate(zip(grid, detector['values']), 1):
            all_grid_cells.append(raw_cell)
            lo, hi, a, b = keys[scoring.native_cell_key(raw_cell)]
            dE = hi - lo
            solid = 2 * math.pi * (math.cos(math.radians(a)) - math.cos(math.radians(b)))
            converted = native_val * fluka_factor
            contrib = converted * dE * solid
            contribs.append(contrib)
            cell_rows.append({
                'cell_index': len(cell_rows) + 1,
                'detector': detector['name'],
                'angle_bin': angle_idx,
                'energy_low_MeV': f'{lo:.5f}',
                'energy_high_MeV': f'{hi:.5f}',
                'energy_width_MeV': repr(dE),
                'angle_low_deg': str(a),
                'angle_high_deg': str(b),
                'solid_angle_sr': repr(solid),
                'native_value_per_primary_GeV_sr': repr(native_val),
                'unit_conversion_GeV_to_MeV': repr(fluka_factor),
                'converted_value_per_primary_MeV_sr': repr(converted),
                'integral_contribution_per_primary': repr(contrib),
            })

    cov_status, cov_detail = scoring.coverage(all_grid_cells, problem)
    require(cov_status['coverage'] == 'exact_once' and cov_detail['cell_count'] == 28
            and cov_detail['required_cells'] == 28
            and not cov_detail['missing'] and not cov_detail['overlap'] and not cov_detail['unmatched'],
            'FLUKA 28-cell exact_once coverage check failed')

    fluka_cells_csv = out / 'fe-m01-g000-f00-fluka-28-cells.csv'
    write_csv(fluka_cells_csv, cell_rows)
    fluka_recomputed = math.fsum(contribs)
    fluka_sequential_sum = sum(contribs)
    require(math.isclose(fluka_recomputed, expected_yield, rel_tol=1e-12, abs_tol=0.0),
            'FLUKA stream 0 recomputed sum mismatch with units/Fe-M01.json fluka_unit_yields[0]')
    return {
        'stream_id': 0,
        'seed': 72690000,
        'primaries': parsed['primaries'],
        'total_weight': parsed['weight'],
        'scoring_card': 'USRYIELD code=124.0 ixa=3 unit=+21 (per-primary double-differential yield)',
        'native_normalization_field': {'user_normalization': 1.0, 'adopted_cross_section_mb': 1.0},
        'unit_conversion_factor_GeV_to_MeV': fluka_factor,
        'per_cell_note': 'Native fort.21 values are already per primary per GeV per sr; 0.001 converts energy density from /GeV to /MeV only.',
        'cell_1_fe01_angle1': {
            'energy_interval_MeV': [35.18064, 41.25585],
            'angle_interval_deg': [0, 45],
            'native_value': parsed['detectors'][0]['values'][0],
            'integral_contribution': contribs[0],
        },
        'recomputed_sum_fsum': fluka_recomputed,
        'recomputed_sum_sequential': fluka_sequential_sum,
        'units_fe_m01_array_value': expected_yield,
        'abs_diff_vs_units_array': abs(fluka_recomputed - expected_yield),
        'tolerances': {
            'float_summation_rel_tol': 1e-12,
            'sequential_vs_fsum_abs_diff': abs(fluka_sequential_sum - fluka_recomputed),
            'native_fort21_print_precision': '5 significant digits (.4E formatting, <=5e-5 rounding scale per printed cell); distinct from float summation rounding.',
        },
        'breakdown_28_cells_csv': fluka_cells_csv.name,
    }


def verify_formal_geant4(formal_dir, problem, expected_yield, out, score_reduction):
    id_tsv = parse_tsv(formal_dir / 'geant4/identity.tsv')
    incl_tsv = parse_tsv(formal_dir / 'geant4/incl-physics.tsv')
    m1_excerpt = load_json(formal_dir / 'geant4/m1-readback-excerpt.json')
    m1_map = dict(entry['text'].split('\t', 1) for entry in m1_excerpt['lines'] if '\t' in entry['text'])
    require(id_tsv['seed'] == '70690000' and id_tsv['requested_events'] == '10000'
            and id_tsv['physics'] == 'QGSP_INCLXX_HP' and incl_tsv['after_cluster'] == '0'
            and m1_map['seed'] == '70690000' and m1_map['requested_events'] == '10000'
            and m1_map['target_proton_cut_mm_bits'] == '0000000000000000',
            'Formal Geant4 u000 host metadata check failed')

    with (formal_dir / 'geant4/primaries.csv').open(newline='') as f:
        prim_rows = list(csv.DictReader(f))
    n_prim = len(prim_rows)
    require(n_prim == 10000 and [int(r['event']) for r in prim_rows] == list(range(10000))
            and all(int(r['pdg']) == 2212 for r in prim_rows),
            'Formal Geant4 u000 primaries.csv completeness check failed')

    def unbits(row, name):
        return struct.unpack('>d', bytes.fromhex(row[name + '_bits']))[0]

    edges = problem['scoring']['energy_edges_MeV']
    angles = problem['scoring']['angle_edges_deg']
    counts = [[0] * (len(angles) - 1) for _ in range(len(edges) - 1)]
    direct_crossings = 0
    total_exits_rows = 0
    with (formal_dir / 'geant4/exits.csv').open(newline='') as f:
        for row in csv.DictReader(f):
            total_exits_rows += 1
            require(0 <= int(row['event']) < n_prim, 'Out-of-range event ID in exits.csv')
            cell = score_reduction.cell(
                int(row['pdg']),
                unbits(row, 'energy_MeV'),
                [unbits(row, k) for k in ('dx', 'dy', 'dz')],
                problem['scoring'],
            )
            if cell is not None:
                direct_crossings += 1
                counts[cell[0]][cell[1]] += 1

    g4_direct = direct_crossings / n_prim
    g4_cell_rows, g4_contribs = [], []
    for ei, (lo, hi) in enumerate(zip(edges, edges[1:])):
        dE = hi - lo
        for ai, (a, b) in enumerate(zip(angles, angles[1:])):
            solid = 2 * math.pi * (math.cos(math.radians(a)) - math.cos(math.radians(b)))
            cnt = counts[ei][ai]
            diff_val = cnt * score_reduction.yield_factor(n_prim, (lo, hi), (a, b))
            contrib = diff_val * dE * solid
            g4_contribs.append(contrib)
            g4_cell_rows.append({
                'cell_index': len(g4_cell_rows) + 1,
                'energy_low_MeV': f'{lo:.5f}',
                'energy_high_MeV': f'{hi:.5f}',
                'energy_width_MeV': repr(dE),
                'angle_low_deg': str(a),
                'angle_high_deg': str(b),
                'solid_angle_sr': repr(solid),
                'window_exit_count': cnt,
                'diff_yield_per_primary_MeV_sr': repr(diff_val),
                'integral_contribution_per_primary': repr(contrib),
            })

    g4_cells_csv = out / 'fe-m01-g000-f00-geant4-28-cells.csv'
    write_csv(g4_cells_csv, g4_cell_rows)
    g4_binned_integrated = math.fsum(g4_contribs)
    require(math.isclose(g4_direct, expected_yield, rel_tol=1e-12, abs_tol=0.0)
            and math.isclose(g4_binned_integrated, expected_yield, rel_tol=1e-12, abs_tol=0.0),
            'Geant4 batch 0 two-path recomputation mismatch with g4_unit_yields[0]')
    return {
        'batch_id': 0,
        'seed': 70690000,
        'primaries': n_prim,
        'exits_csv_rows': total_exits_rows,
        'window_neutron_exits': direct_crossings,
        'direct_crossing_over_primaries': g4_direct,
        'binned_28_cells_integrated': g4_binned_integrated,
        'units_fe_m01_array_value': expected_yield,
        'abs_diff_direct_vs_array': abs(g4_direct - expected_yield),
        'abs_diff_binned_vs_array': abs(g4_binned_integrated - expected_yield),
        'breakdown_28_cells_csv': g4_cells_csv.name,
    }


def verify_formal_example(out):
    usryield, scoring, score_reduction = load_settings_sources()
    problem = load_json(ROOT / 'contracts/problem-fe.json')
    units_data = load_json(ROOT / 'units/Fe-M01.json')
    formal_dir = ROOT / 'example/formal-g000-f00'
    prov = load_json(formal_dir / 'pair-p000-provenance.json')
    idx_path = verify_formal_index_and_files(formal_dir, units_data, prov)
    return {
        'unit_id': 'Fe-M01/g000-f00',
        'array_index': 0,
        'statistical_independence': prov['statistical_independence'],
        'kernel_index_sha256': sha(idx_path),
        'p000_kernel_source_sha256': prov['pair_kernel_excerpt']['source_sha256'],
        'fluka': verify_formal_fluka(formal_dir, problem, units_data['fluka_unit_yields'][0], out, usryield, scoring),
        'geant4': verify_formal_geant4(formal_dir, problem, units_data['g4_unit_yields'][0], out, score_reduction),
        'scope_note': 'Native cell/CSV recomputation covers only this single Fe-M01/g000-f00 formal pair (Geant4 batch 0 and FLUKA stream 0); remaining 799 Geant4 batches and 39 FLUKA streams remain unit-level yield array checks.',
    }


def formal_ratios():
    plan = load_json(ROOT / 'records/formal-transport-plan.json')
    methods_path = ROOT / plan['methods_document']
    require(sha(methods_path) == plan['methods_sha256'], 'T069 METHODS.md hash mismatch')
    methods_text = methods_path.read_text()
    require('1.3147499288226958' in methods_text and '0.0143609212300579' in methods_text
            and '1.315±0.029' in methods_text and '1.10' in methods_text,
            'T069 METHODS.md reference or criteria text missing')
    freeze1 = load_json(ROOT / 'reviews/T069/offline-transport-20261114/freeze.json')
    freeze2 = load_json(ROOT / 'reviews/T069/repair-20261114/freeze.json')
    require(plan['source_commit'] == 'da22a9bbc74fb4e5c3a98df9be5d29f77ac5a185'
            and plan['kernel_version'] == freeze1['kernel_version'] == 't069/1-transport-unreviewed'
            and plan['previous_checker_commit'] == freeze1['previous_checker'] == '66d3acb04ae834e3dd4579d2239ede9a8bac17ba'
            and plan['recheck_kernel_version'] == freeze2['kernel_version'] == 't069/2-exact-seed-unreviewed'
            and freeze2['previous_checker'] == plan['source_commit']
            and plan['recheck_source_commit'] == 'e2fb2d95b7cfbde9584b866f0c20b8516b24e66c',
            'T069 freeze commit or kernel version mismatch')
    require(freeze1['files']['experiments/recorded-input-transport/METHODS.md'] == plan['methods_sha256']
            and freeze2['files']['experiments/recorded-input-transport/METHODS.md'] == plan['methods_sha256'],
            'T069 METHODS.md freeze binding')
    sel_ref = plan['selected_reference']
    both_on = [r for r in plan['reference_table']['ratios']
               if r['numerator'] == 'g4-default' and r['denominator'] == 'fluka-coalesce']
    require(len(both_on) == 1
            and math.isclose(sel_ref['ratio'], both_on[0]['ratio'], rel_tol=1e-15)
            and math.isclose(sel_ref['mc_se'], both_on[0]['mc_se'], rel_tol=1e-15)
            and sel_ref['source_sha256'] == '8612e1108a5ddf435aa1b63b40e0f32433ed142ab138a2f38eb104e0c6ace210',
            'Pre-selected both-ON reference mismatch')

    results = {}
    for name in ('Fe-M01', 'Fe-C01'):
        data = load_json(ROOT / 'units' / (name + '.json'))
        g, f = data['g4_unit_yields'], data['fluka_unit_yields']
        require((len(g), len(f)) == (800, 40), 'Statistical units: ' + name)
        require(data['source_commit'] == plan['source_commit']
                and data['archived_paired_setting_records'] == 800, 'Transport metadata: ' + name)
        gm, fm = statistics.mean(g), statistics.mean(f)
        ratio = gm / fm
        se = math.sqrt(statistics.variance(g) / len(g) / fm**2 +
                       gm**2 * statistics.variance(f) / len(f) / fm**4)
        for field, actual in [('ratio', ratio), ('mc_se', se),
                              ('low2se', ratio - 2*se), ('high2se', ratio + 2*se)]:
            require(math.isclose(actual, data[field], rel_tol=1e-12), name + ': ' + field)
        ref = data['reference']
        require(ref == {k: sel_ref[k] for k in ('ratio', 'mc_se', 'source_sha256')},
                'Reference binding: ' + name)
        compatible = abs(ratio - ref['ratio']) <= 2 * math.hypot(se, ref['mc_se'])
        require(compatible and ratio - 2*se > 1.10, 'Reference criteria: ' + name)
        results[name] = {
            'g4_mean': gm,
            'fluka_mean': fm,
            'ratio': ratio,
            'two_mcse': 2*se,
            'low2se': ratio - 2*se,
            'reference_compatible': compatible,
            'above_1_10': ratio - 2*se > 1.10,
        }
    return results


def check_records():
    failed = {}
    for sample in ('Fe-M01', 'Fe-C01', 'Pb-M01', 'Pb-C01'):
        first, last = [load_json(ROOT / 'records' / f'{sample}-r{r}.json') for r in (1, 2)]
        require(len(first['table']) == len(last['table']) == 22 and not first['passed'] and last['passed'],
                'Round status: ' + sample)
        require(all(row['passed'] for row in last['table']), 'Final row: ' + sample)
        require(len(last['intrinsic']) == 19 and all(not row['confirmed'] for row in last['intrinsic']),
                'Intrinsic register: ' + sample)
        failed[sample] = [row['id'] for row in first['table'] if not row['passed']]
        require('P26' in failed[sample], 'Recoil range cut P26 in round 1: ' + sample)
        changed = sorted(k for k, v in first['answer_field_sha256'].items()
                         if v != last['answer_field_sha256'][k])
        require(changed == ['answers_json', 'differences_json'], 'Changed answer fields: ' + sample)
    final_example = load_json(ROOT / 'records/Fe-M01-r2.json')['answer_field_sha256']
    for key, filename in [('geometry_gdml', 'geometry.gdml'), ('beam_mac', 'beam.mac'),
                          ('physics_txt', 'physics.txt'), ('host_options_json', 'host-options.json'),
                          ('fluka_inp', 'fluka.inp')]:
        require(sha(ROOT / 'example' / filename) == final_example[key], 'Example file: ' + filename)
    fe_m01_inp_lines = (ROOT / 'example/fluka.inp').read_text().splitlines()
    fe_m01_uy = [(i, l) for i, l in enumerate(fe_m01_inp_lines, 1) if l.startswith('USRYIELD')]
    require(len(fe_m01_uy) == 14 and [i for i, _ in fe_m01_uy] == list(range(35, 49))
            and all('124.0' in l and '21.0' in l for _, l in fe_m01_uy[0::2])
            and all('4.0' in l and '3.0&' in l for _, l in fe_m01_uy[1::2]),
            'Fe-M01 example/fluka.inp lines 35-48 (7 detectors of code 124.0, 4 angle bins each)')

    readback = load_json(ROOT / 'example/readback.json')
    incl_lines = [entry['text'] for entry in readback['incl-physics.tsv']]
    m1_lines = [entry['text'] for entry in readback['m1-readback.tsv']]
    fluka_lines = [entry['text'] for entry in readback['input001.out']]
    require('mode\tdefault' in incl_lines and 'configuration_verified\t1' in incl_lines
            and 'before_cluster\t0' in incl_lines and 'after_cluster\t0' in incl_lines,
            'Fe-M01 incl-physics.tsv cluster readback excerpt')
    require('target_gamma_cut_mm_bits\t3fb999999999999a' in m1_lines
            and 'target_electron_cut_mm_bits\t3fb999999999999a' in m1_lines
            and 'target_positron_cut_mm_bits\t3fb999999999999a' in m1_lines
            and 'target_proton_cut_mm_bits\t0000000000000000' in m1_lines
            and 'material_temperature_K_bits\t4072800000000000' in m1_lines
            and len([s for s in m1_lines if s.startswith('neutron_process_')]) == 6,
            'Fe-M01 m1-readback.tsv EM cut, proton cut and neutron process excerpt')
    require(any('COALESCE' in s for s in fluka_lines)
            and any('Ae =' in s and '0.611000' in s for s in fluka_lines)
            and any('Ap =' in s and '0.100000' in s for s in fluka_lines)
            and any('1.000E-14 GeV' in s for s in fluka_lines)
            and any('adopted cross section (if any):  1.0000E+00 mb' in s for s in fluka_lines),
            'Fe-M01 FLUKA input001.out readback excerpt')

    groups = load_json(ROOT / 'records/early-explicit-request-summary.json')['groups']
    totals = [g for g in groups if g['scope'] == 'all_complete_texts']
    require({g['group'] for g in totals} == {'M', 'C'} and len(totals) == 2
            and all(g['n'] == g['declared_same'] == g['P03_archived_readback_confirmed'] == 4
                    for g in totals), 'Earlier eight explicit requests')

    early = load_json(ROOT / 'records/early-Fe-scoring.json')
    old_phys = early['old_physics_record']
    old_table_path = ROOT / old_phys['bundled_table']
    old_table_text = old_table_path.read_text()
    require(sha(old_table_path) == old_phys['bundled_table_sha256']
            and old_phys['passed'] is True
            and old_phys['physics_row_count'] == old_phys['physics_passed_count'] == 13
            and old_phys['intrinsic_count'] == 19
            and old_phys['readback_errors'] == []
            and len(old_phys['physics_rows']) == 13
            and all(r['passed'] for r in old_phys['physics_rows']),
            'Early input M02-r5 old 13/13 physics kernel record')
    expected_p_ids = ('P01', 'P02', 'P03', 'P22', 'P23', 'P24', 'P25', 'P26', 'P28', 'P29', 'P30', 'P31', 'P32')
    expected_intr_ids = tuple(f'P{i:02d}' for i in range(4, 22)) + ('P27',)
    require('通过：True；满足行数：13/13。' in old_table_text
            and '读回失败：`[]`' in old_table_text
            and all(f'| {pid} ' in old_table_text for pid in expected_p_ids)
            and all(f'| {iid} |' in old_table_text for iid in expected_intr_ids),
            'Early input M02-r5 bound table.md 13 physics + 19 intrinsic rows')
    rows = {row['id']: row['sides']['fluka'] for row in early['replay_rows']}
    require(all(rows[k]['conforms_to_question'] for k in ('S02', 'S03', 'S04', 'S09')),
            'Early input particle, geometry or bins')
    require(rows['S05']['effective'] == 'd2sigma_dE_dOmega' and not rows['S05']['conforms_to_question']
            and rows['S06']['effective']['cross_section_factor'] == 20.534, 'Early input observable')
    old_cards = early['first_detector_input_cards']
    require(len(old_cards) == 2 and [c['line'] for c in old_cards] == [35, 36]
            and '1401.0' in old_cards[0]['text'] and '-21.0' in old_cards[0]['text']
            and old_cards[1]['text'].endswith('1.0&'),
            'Early input first detector cards (1401.0, unit -21, ixa=1)')
    native = [line['text'] for line in early['native_output_excerpt']]
    require(any('2.0534E+01 mb' in line for line in native) and '每初级' in early['model_statement'],
            'Early input evidence excerpt')

    r1 = early['continuation_round_1']
    r1_cards = r1['first_detector_input_cards']
    require(not r1['passed'] and r1['fluka_exit'] == 2
            and r1['usryield_line_count'] == 56
            and len(r1_cards) == 2 and [c['line'] for c in r1_cards] == [35, 36]
            and '2401.0' in r1_cards[0]['text'] and '21.0' in r1_cards[0]['text']
            and r1_cards[1]['text'].endswith('3.0&')
            and all(v == 'no_readback' for v in r1['scoring_rows_status'].values())
            and any('Bad integer for item 2 in list input' in line['text']
                    for line in r1['native_error_excerpt']),
            'Fe-M02 continuation round 1 already changed cards (2401.0, unit +21, ixa=3) and stopped at geometry header')

    r2 = early['continuation_round_2']
    r2_rows = {row['id']: row for row in r2['scoring_rows']}
    require(r2['passed'] and r2['initialized_detector_count'] == 28
            and r2['all_detectors_unit_21_and_factor_1'] is True
            and r2['usryield_line_count'] == 56
            and r2['usryield_lines_sha256'] == r1['usryield_lines_sha256']
            and r2['first_detector_input_cards'] == r1_cards
            and r2['usryield_cards_identical_to_round_1'] is True
            and r2['only_fluka_inp_line_changed_from_round_1']['line'] == 7
            and set(r2_rows) == {f'S0{i}' for i in range(1, 10)}
            and all(row['passed'] and row['aligned'] for row in r2_rows.values()),
            'Fe-M02 continuation round 2 identical 56 USRYIELD lines and passed scoring rows')
    require(r2_rows['S05']['sides']['fluka']['effective'] == 'd2N_dE_dOmega_per_primary'
            and r2_rows['S06']['sides']['fluka']['effective']['cross_section_factor'] == 1.0
            and r2_rows['S07']['sides']['fluka']['effective']['factor'] == 0.001
            and r2_rows['S08']['sides']['fluka']['effective'] == 'formatted_USRYIELD',
            'Fe-M02 continuation round 2 effective FLUKA values')
    r2_native = [line['text'] for line in r2['native_output_excerpt']]
    require(any('USRYIELD   2401.' in line and '21.00' in line for line in r2_native)
            and any('3.000      &' in line for line in r2_native)
            and any('user normalization:  1.0000E+00, adopted cross section (if any):  1.0000E+00 mb' in line
                    for line in r2_native)
            and any('data will be printed on unit   21' in line for line in r2_native),
            'Fe-M02 continuation round 2 native output excerpt')
    return failed


def verify_figures(responses):
    fig_dir = ROOT / 'paper/setting-semantics/figures/responses'
    with (fig_dir / 'data.csv').open(newline='') as f:
        figure = list(csv.DictReader(f))
    with (fig_dir / 'appendix-b1-data.csv').open(newline='') as f:
        b1_fig = list(csv.DictReader(f))
    with (fig_dir / 'appendix-b2-data.csv').open(newline='') as f:
        b2_fig = list(csv.DictReader(f))
    require(len(figure) == 20 and len(b1_fig) == 10 and len(b2_fig) == 10,
            'Figure point counts (total 20, appendix B1 10, appendix B2 10)')
    require([r for r in figure if r['figure'] == 'appendix-B1'] == b1_fig
            and [r for r in figure if r['figure'] == 'appendix-B2'] == b2_fig,
            'Appendix B1/B2 figure split consistency')
    require(all(r['contrast'] == 'cluster-off-minus-default' for r in b1_fig)
            and all(r['contrast'] == 'incl42-minus-default' for r in b2_fig),
            'Appendix B1/B2 contrast assignment')
    for point in figure:
        row = responses[int(point['original_csv_row']) - 2]
        require(row['band'] == 'middle' and all(row[k] == point[k] for k in ('region', 'metric', 'contrast')),
                'Figure point identity')
        for k, field, factor in [('response_percent', 'relative_value', 100),
                                 ('two_mcse_percentage_points', 'relative_mc_se', 200)]:
            require(math.isclose(float(point[k]), factor*float(row[field]), rel_tol=1e-10), 'Figure: ' + k)
    return {'total': len(figure), 'appendix_B1': len(b1_fig), 'appendix_B2': len(b2_fig)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    out = parser.parse_args().out.resolve()
    require(not out.is_relative_to(ROOT) and not out.exists(), 'Choose a new directory outside the package')
    manifest = load_json(ROOT / 'MANIFEST.json')
    for filename, expected in manifest['files'].items():
        require(sha(ROOT / filename) == expected, 'File hash: ' + filename)
    ratios = formal_ratios()
    failed = check_records()
    path = ROOT / 'paper/setting-semantics/reproduce.py'
    spec = importlib.util.spec_from_file_location('tables', path)
    tables = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tables)
    out.mkdir(parents=True)
    formal_example = verify_formal_example(out)
    summary = tables.response_tables(out, load_units=read_units)
    with (out / 'g4-responses-display.csv').open(newline='') as f:
        responses = list(csv.DictReader(f))
    fig_counts = verify_figures(responses)
    report = {
        'status': 'passed',
        'files_checked': len(manifest['files']),
        'formal_example_fe_m01_g000_f00': formal_example,
        'formal_transport': ratios,
        'pre_transport_reference': {
            'cell': '两边都开：G4 默认 / FLUKA COALESCE',
            'ratio': 1.3147499288226958,
            'mc_se': 0.0143609212300579,
            'source_commit': 'da22a9bbc74fb4e5c3a98df9be5d29f77ac5a185',
            'recheck_source_commit': 'e2fb2d95b7cfbde9584b866f0c20b8516b24e66c',
            'archived_paired_records_total': 1600,
        },
        'response_rows': summary['rows'],
        'response_rows_total': sum(summary['rows'].values()),
        'covariance_cells': summary['covariance_cells'],
        'figure_points': fig_counts['total'],
        'figure_points_split': fig_counts,
        'first_round_failed_rows': failed,
        'fe_m01_readback_excerpts': 'verified (START 1.0 initialization: incl-physics.tsv after_cluster=0, m1-readback.tsv EM cut 0.1 mm & proton cut 0 mm, FLUKA COALESCE, Ae=0.611 MeV/Ap=0.1 MeV, and USRYIELD)',
        'early_scoring_case': 'verified (M02-r5 old 13/13 physics pass table, replay 20.534 mb cross-section, round-1 geometry-header stop, round-2 unit 21 and factor 1.0 mb yield confirmation)',
        'scope': 'Native per-cell/CSV recomputation for Fe-M01/g000-f00 only; unit-level reproduction and retained record check for remaining 799 G4 batches and 39 FLUKA streams; no new Monte Carlo or full kernel replay.',
    }
    (out / 'verification.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
