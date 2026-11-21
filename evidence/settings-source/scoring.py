"""Scoring semantics from native artifacts, with an explicit question contract.

No transport or spectra are evaluated here. USRYIELD summary/tally regexes are
reused from the installed-format reader. Input declarations never supply an
effective value. Print-resolution comparisons are kept separate from integers.
"""
import copy
import hashlib
import json
import math
from pathlib import Path
import re

from . import base_kernel as base
import usryield

TITLES = {
    'S01': '计分安装', 'S02': '能量×极角格子覆盖，无缺口无重叠',
    'S03': '计分粒子', 'S04': '穿越面、方向与极角参考轴',
    'S05': '结果量与产额类型', 'S06': '实际归一',
    'S07': '到题面结果量的换算', 'S08': '输出可读性',
    'S09': '两个变量的编码、含义与数值单位',
}
G4_BINARY = '0fcc6efaf81480236abe5a9b47c86bdb7cb01ff924d3fbc3c17afff5ff1e33d6'
EXITS = 'event,track,parent,pdg,energy_MeV_bits,x_mm_bits,y_mm_bits,z_mm_bits,dx_bits,dy_bits,dz_bits,time_ns_bits'


def expected(problem, side):
    s = problem['scoring']
    return {
        'S01': True,
        'S02': {'energy_edges_MeV': s['energy_edges_MeV'],
                'angle_edges_deg': s['angle_edges_deg'], 'coverage': 'exact_once'},
        'S03': s['particle'],
        'S04': s['surface'],
        'S05': s['quantity'],
        'S06': {'per_primary': True, 'per_energy_width': True, 'per_solid_angle': True,
                'user_normalization': 1.0, 'cross_section_factor': 1.0},
        'S07': {'native_energy_unit': 'MeV' if side == 'geant4' else 'GeV',
                'output_unit': s['unit'], 'factor': 1.0 if side == 'geant4' else 0.001},
        'S08': 'event_CSV_plus_fixed_reducer' if side == 'geant4' else 'formatted_USRYIELD',
        'S09': {'energy': 'kinetic', 'angle': 'laboratory_polar',
                'energy_readback_unit': 'MeV' if side == 'geant4' else 'GeV',
                'angle_readback_unit': 'direction_cosine' if side == 'geant4' else 'radian'},
    }


def equivalent(a, b, tol=6e-5):
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(equivalent(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(equivalent(a[k], b[k], tol) for k in a)
    return base.close(a, b, tol, 1e-10)


def summary(raw):
    """Actual initialized estimators, including negative-unit binary output."""
    heads = list(usryield.HEAD.finditer(raw)); result = []
    for i, h in enumerate(heads):
        block = raw[h.end():heads[i+1].start() if i+1 < len(heads) else len(raw)]
        f, s, n = usryield.FIRST.search(block), usryield.SECOND.search(block), usryield.NORM.search(block)
        unit = re.search(r'data will be printed on unit\s+(-?\d+)', block)
        if not all((f, s, n, unit)):
            raise ValueError('incomplete initialized estimator '+h[2].strip())
        kinds = dict(usryield.KIND.findall(block))
        result.append({'number': int(h[1]), 'name': h[2].strip(), 'particle': int(h[3]),
            'from_region': int(h[4]), 'to_region': int(h[5]), 'unit': int(unit[1]),
            'normalization': float(n[1]), 'cross_section_mb': float(n[2]),
            'first': {'scale': f[1], 'low': float(f[2]), 'high': float(f[3]),
                      'bins': int(f[4]), 'quantity': kinds.get('1st')},
            'second': {'low': float(s[1]), 'high': float(s[2]), 'quantity': kinds.get('2nd')}})
    if len({d['number'] for d in result}) != len(result) or len({d['name'] for d in result}) != len(result):
        raise ValueError('duplicate runtime detector identity')
    return result


def echoed_yields(raw):
    result = {}; pending = None; controls = []
    for line in raw.splitlines():
        if 'Next control card' not in line:
            continue
        words = line.split('*****')[-1].split()
        # Split at the known marker instead: the trailing spaces are not data.
        words = re.split(r'Next control card\s+\*+\s*', line, maxsplit=1)[-1].split()
        if not words:
            continue
        key = words[0]
        if key in ('AUXSCORE', 'DCYSCORE', 'USERWEIG', 'USRICALL'):
            controls.append(line.strip())
        if key != 'USRYIELD':
            if pending: raise ValueError('unpaired runtime USRYIELD echo')
            continue
        if len(words) != 8:
            raise ValueError('unrecognized runtime USRYIELD echo')
        values = [float(v.replace('D', 'E')) for v in words[1:7]]
        name = words[7]
        if name == 'BEAMDEF':
            controls.append(line.strip()); continue
        if name == '&':
            if pending is None: raise ValueError('orphan runtime continuation')
            first, name = pending
            if name in result: raise ValueError('redefined detector not supported')
            result[name] = {'first': first, 'second': values}; pending = None
        else:
            if pending: raise ValueError('missing runtime continuation')
            pending = (values, name)
    if pending: raise ValueError('missing runtime continuation')
    return result, controls


def variable_cells(d):
    f, s = d['first'], d['second']
    kinds = (f['quantity'], s['quantity'])
    angle, energy = 'Laboratory Angle (radians)', 'Laboratory Kinetic Energy'
    if set(kinds) != {angle, energy}:
        raise ValueError('unsupported variables: '+repr(kinds))
    n = f['bins']
    if not 1 <= n <= 10000 or not f['low'] < f['high'] or not s['low'] < s['high']:
        raise ValueError('invalid runtime bin limits/count')
    if f['scale'] == 'linear':
        edges = [f['low']+(f['high']-f['low'])*i/n for i in range(n+1)]
    else:
        raise ValueError('unsupported bin scale')
    cells = []
    for lo, hi in zip(edges, edges[1:]):
        e, a = ((lo, hi), (s['low'], s['high'])) if kinds[0] == energy else ((s['low'], s['high']), (lo, hi))
        cells.append([e[0]*1000, e[1]*1000, math.degrees(a[0]), math.degrees(a[1])])
    return cells


def native_cell_key(cell):
    """FLUKA initializes limits in GeV/radians, printing five significant digits.

    Equality is at that printed precision, not a relative tolerance band.
    Values reconstructed from grouped bins must meet the same rule.
    """
    return tuple(format(v, '.4E') for v in
                 (cell[0]/1000, cell[1]/1000, math.radians(cell[2]), math.radians(cell[3])))


def coverage(cells, problem):
    """Recognize exact cells in any order/grouping at native print resolution.

    A coarser cell is not split by assumption; extra/finer cells are not silently
    discarded. This preserves the question's actual bin boundaries.
    """
    s = problem['scoring']; expected_cells = [[a,b,c,d]
        for a,b in zip(s['energy_edges_MeV'], s['energy_edges_MeV'][1:])
        for c,d in zip(s['angle_edges_deg'], s['angle_edges_deg'][1:])]
    counts = [0]*len(expected_cells); unmatched = []
    for cell in cells:
        matches = [i for i, want in enumerate(expected_cells) if native_cell_key(cell) == native_cell_key(want)]
        if len(matches) == 1: counts[matches[0]] += 1
        else: unmatched.append(cell)
    detail = {'cell_count': len(cells), 'required_cells': len(expected_cells),
              'missing': [c for c,n in zip(expected_cells,counts) if n == 0],
              'overlap': [c for c,n in zip(expected_cells,counts) if n > 1], 'unmatched': unmatched}
    ok = not any(detail[k] for k in ('missing','overlap','unmatched'))
    return (expected(problem,'fluka')['S02'] if ok else {'coverage':'mismatch','cell_count':len(cells)}), detail


def fluka(directory, problem, physics_evidence, execution=None):
    values = {k: None for k in TITLES}; detail = {'errors': []}
    rawpath = directory/'input001.out'; raw = rawpath.read_text(errors='replace')
    if not physics_evidence or base.sha(rawpath) != physics_evidence['raw_sha256']:
        raise ValueError('no complete, bound FLUKA run')
    ds = summary(raw); echoes, controls = echoed_yields(raw)
    detail.update(summary=ds, echoed_yields=echoes, unsupported_controls=controls, raw_sha256=base.sha(rawpath))
    values['S01'] = bool(ds)
    if not ds:
        if 'No user yield' not in raw: detail['errors'].append('no initialized USRYIELD table')
        return values, detail
    if set(echoes) != {d['name'] for d in ds}:
        raise ValueError('runtime echo/initialized-detector inventory differs')
    particles = sorted({d['particle'] for d in ds})
    values['S03'] = 'neutron' if particles == [8] and not controls else {'particle_ids':particles,'unread_controls':bool(controls)}
    regions = physics_evidence['em_tables']['regions']; surface = problem['scoring']['surface']
    source = regions[surface['from']]['region_id']; dest = regions[surface['to']]['region_id']
    boundary = all(d['from_region']==source and d['to_region']==dest for d in ds)
    if boundary and not controls:
        # Base question checks establish TARGET geometry, surrounding vacuum and beam axis.
        values['S04'] = copy.deepcopy(surface)
    else:
        values['S04'] = {'region_pairs': sorted([list(p) for p in {(d['from_region'],d['to_region']) for d in ds}]),
                         'unread_controls':bool(controls)}
    cells = []; variables_ok = not controls; kinds = []; echo_issues = []
    for d in ds:
        e = echoes[d['name']]; f, c = e['first'], e['second']
        if not all(math.isfinite(x) for x in f+c) or not all(x.is_integer() for x in [f[0],f[1],f[2],f[3],f[4],c[2],c[5]]):
            raise ValueError('nonfinite/nonintegral runtime scoring control')
        # Decode supported positive linear or simple negative-first-variable form.
        code = int(f[0]); ia = round(code/100); ie = code-ia*100
        ixa = int(c[5]) % 100; ixm = int(c[5]) // 100
        kinds.append(ixa)
        if {abs(ie),abs(ia)} not in ({1,14},{1,24}) or ixm != 0:
            variables_ok = False; echo_issues.append(d['name']+': unsupported variable/material code')
        order = (abs(ie),abs(ia))
        want_quantities = tuple('Laboratory Kinetic Energy' if v==1 else 'Laboratory Angle (radians)' for v in order)
        if want_quantities != (d['first']['quantity'],d['second']['quantity']): variables_ok = False
        for observed, original in zip([d['particle'],d['unit'],d['from_region'],d['to_region']], f[1:5]):
            if observed != original: variables_ok = False
        # Echo units are chosen by codes; initialized limits are radians/GeV.
        converted = [c[1],c[0],c[4],c[3]]
        if order[0]==24: converted[:2] = [math.radians(x) for x in converted[:2]]
        if order[1]==24: converted[2:] = [math.radians(x) for x in converted[2:]]
        limits = [d['first']['low'],d['first']['high'],d['second']['low'],d['second']['high']]
        if not equivalent(limits, converted,6e-4) or d['first']['bins']!=int(c[2]): variables_ok=False
        if d['first']['scale'] != ('linear' if ie>0 else 'logar.'): variables_ok=False
        try: cells.extend(variable_cells(d))
        except ValueError as err: variables_ok=False; echo_issues.append(str(err))
    values['S02'], detail['coverage'] = coverage(cells,problem)
    detail['variable_errors'] = echo_issues
    if variables_ok: values['S09'] = expected(problem,'fluka')['S09']
    else: values['S09'] = {'recognized':False}
    distinct_kinds = sorted(set(kinds)); detail['ixa'] = distinct_kinds
    values['S05'] = ('d2N_dE_dOmega_per_primary' if distinct_kinds==[3] else
                     'd2sigma_dE_dOmega' if distinct_kinds==[1] else {'ixa':distinct_kinds})
    norms = {(d['normalization'],d['cross_section_mb']) for d in ds}
    if len(norms)==1 and variables_ok and not controls:
        norm, cross = next(iter(norms))
        values['S06'] = {'per_primary':True,'per_energy_width':True,'per_solid_angle':True,
                         'user_normalization':norm,'cross_section_factor':cross}
    if distinct_kinds==[3] and variables_ok and not controls:
        values['S07'] = expected(problem,'fluka')['S07']
    units = sorted({d['unit'] for d in ds}); detail['output_units'] = units
    if any(u<0 for u in units):
        values['S08'] = 'unformatted_unsupported'
    else:
        try:
            tallies = []
            for u in units:
                assert 21<=u<=99, 'unsupported output unit'
                p = directory/f'input001_fort.{u}'
                parsed = usryield.read(p.read_text(encoding='ascii'))
                count = 1 if execution is None else execution['fluka']['histories']
                assert parsed['primaries']==count and equivalent(parsed['weight'],count), 'primary normalization'
                for d in parsed['detectors']:
                    d.pop('values'); d['unit']=u; tallies.append(d)
            assert len(tallies)==len(ds), 'tally detector count'
            for d in ds:
                found = [t for t in tallies if t['name']==d['name'] and t['unit']==d['unit']]
                assert len(found)==1 and equivalent(d,found[0]), 'summary/tally mismatch: '+d['name']
            values['S08'] = 'formatted_USRYIELD'
            detail['tally_files'] = {f'input001_fort.{u}':base.sha(directory/f'input001_fort.{u}') for u in units}
        except (OSError,UnicodeError,ValueError,AssertionError) as err:
            values['S08'] = 'unreadable'; detail['errors'].append(repr(err))
    return values, detail


def geant4(directory, problem, physics_evidence, runtime, host, execution=None):
    """Zero-event schema + fixed binary/source/reducer identity, not sampled counts."""
    if not physics_evidence: raise ValueError('no complete Geant4 question readback')
    plan = json.loads((directory.parent/'plan.json').read_text())
    if plan['g4_sha256'] != G4_BINARY or host['binary_sha256'] != G4_BINARY:
        raise ValueError('unregistered Geant4 scorer binary')
    count = 0 if execution is None else execution['geant4']['events']
    if runtime['argv'][0] != plan['g4_binary'] or runtime['argv'][4]!='Target' or runtime['argv'][5]!=str(count):
        raise ValueError('host scorer target/event binding')
    reducer = Path(__file__).with_name('score_reduction.py')
    if base.sha(reducer) != host['reducer_sha256']: raise ValueError('host reducer changed')
    if base.sha(base.ROOT/host['source_relative_path']) != host['source_sha256']:
        raise ValueError('registered host source changed')
    output = directory/'output'
    if execution is not None:
        from . import transport_scope
        complete = transport_scope.g4_csv(output, count)
        return expected(problem,'geant4'), {'basis':'fixed binary/source and host reducer; complete event CSVs',
            'host':host,'binning':expected(problem,'geant4')['S02'], 'transport_completeness':complete,
            'scope':'C++ writes unbinned crossings; Python applies unchanged contract bins and normalization. Spectrum estimation is separate.'}
    for name, header in [('exits.csv',EXITS),('events.csv','event,steps,exits,neutron_exits,edep_MeV_bits'),
                         ('primaries.csv','event,pdg,x_mm_bits,y_mm_bits,z_mm_bits,px_MeV_bits,py_MeV_bits,pz_MeV_bits')]:
        if (output/name).read_text().splitlines()!=[header]: raise ValueError('zero-event scorer output schema')
    return expected(problem,'geant4'), {'basis':'fixed binary/source and host reducer; zero-event CSV schemas',
        'host':host,'binning':expected(problem,'geant4')['S02'],
        'scope':'C++ writes unbinned crossings; Python applies contract bins and normalization. No event response measured.'}


def append(record, answer, root, problem, host, execution=None):
    declared = json.loads(answer.get('answers_json','{}')); actual={}; evidence={}
    for side, reader in [('fluka',fluka),('geant4',geant4)]:
        try:
            args=[Path(root)/side,problem,record['evidence'].get(side)]
            if side=='geant4': args += [record['runtime'][side],host]
            actual[side],evidence[side]=reader(*args,execution=execution)
        except Exception as err:
            actual[side]={};evidence[side]={'error':repr(err)}
    record['non_scoring_passed']=record['passed']
    for pid,title in TITLES.items():
        row={'id':pid,'title':title,'classification':'question prerequisite','align':True,'sides':{}}
        for side in ('geant4','fluka'):
            effective=actual[side].get(pid); claimed=declared.get(pid,{}).get(side) if isinstance(declared.get(pid),dict) else None
            tolerance = 1e-12 if pid == 'S02' else 6e-5
            confirmed=equivalent(effective,claimed,tolerance); conforms=equivalent(effective,expected(problem,side)[pid],tolerance)
            row['sides'][side]={'effective':effective,'answer':claimed,'confirmed':confirmed,
                'conforms_to_question':conforms,'basis':'fixed host' if side=='geant4' else 'native output',
                'status':'no_readback' if effective is None else 'question_mismatch' if not conforms else 'answer_missing_or_mismatch' if not confirmed else 'confirmed'}
        row['passed']=all(s['confirmed'] and s['conforms_to_question'] for s in row['sides'].values())
        row['aligned']=all(s['conforms_to_question'] for s in row['sides'].values())
        row['registered_difference']=None
        record['table'].append(row)
    record['evidence']['scoring']=evidence
    record['passed']=record['non_scoring_passed'] and all(r['passed'] for r in record['table'])
    record['kernel_version']='t067/2-offline-repair-unreviewed'
    if execution is not None:
        from . import transport_scope
        record['kernel_version'] = transport_scope.VERSION
        record['execution_contract'] = execution
    return record
