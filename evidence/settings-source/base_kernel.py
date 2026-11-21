"""T064 parameterized reader derived from the unchanged T063 candidate.
Question constants are explicit arguments; no input repair or automatic answers.
"""
if not __debug__:
    raise RuntimeError('artifact kernel requires assertion checks; Python -O is unsupported')

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import struct
import sys
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'src'), str(ROOT/'adapters/geant4'), str(ROOT/'adapters/fluka')]
from pilot import precheck
import runtime_evidence as g4read
from .hp_reader import hp_coverage
from .material_readback import fluka_target
from .input_scope import (gdml_topology_scope, gps_direction_scope, fluka_field_coverage,
                          macro_lines, fudgem_echo_prerequisite, fluka_text_scope,
                          blank_input_prerequisites)
from . import mass_readback
from electromagnetic import electromagnetic_cuts

CHECKLIST = json.loads((HERE/'checklist.json').read_text())
FIELDS = {'geometry_gdml': 'geant4/geometry.gdml', 'beam_mac': 'geant4/beam.mac',
          'physics_txt': 'geant4/physics.txt', 'fluka_inp': 'fluka/original.inp'}
COVERED = {'TITLE','DEFAULTS','BEAM','BEAMPOS','GEOBEGIN','GEOEND','MATERIAL','COMPOUND',
           'ASSIGNMA','LOW-PWXS','LOW-MAT','PHOTONUC','PHYSICS','PART-THR','EMFCUT',
           'RANDOMIZ','START','STOP','USRYIELD','USRBDX','AUXSCORE'}
NUM = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cards(text):
    result = []
    inside = title = False
    for line in text.splitlines():
        if title:
            title = False
            continue
        if not line.strip() or line.startswith('*'):
            continue
        key = line[:10].strip().upper()[:8]
        if key == 'TITLE':
            title = True
        if key == 'GEOBEGIN':
            inside = True
        elif key == 'GEOEND':
            inside = False
        elif inside:
            continue
        result.append((key, [line[i:i+10].strip() for i in range(10,70,10)], line[70:78].strip(), line))
    return result


def validate(answer):
    errors = []
    try:
        assert set(answer) == set(FIELDS) | {'host_options_json','answers_json','differences_json'}
        assert all(isinstance(v,str) and len(v.encode()) <= 256*1024 for v in answer.values())
        opts = json.loads(answer['host_options_json'])
        assert isinstance(opts,dict) and set(opts) <= {'range_cut_mm','gamma_nuclear','incl_physics'}
        if 'range_cut_mm' in opts:
            x = opts['range_cut_mm']
            assert type(x) in (int,float) and math.isfinite(x) and 0 < x <= 1000
        assert opts.get('gamma_nuclear','default') in ('default','off')
        assert opts.get('incl_physics','default') in ('default','incl42','cluster-off')
        for name in ('answers_json','differences_json'):
            assert isinstance(json.loads(answer[name]),dict)
        macro_lines(answer['beam_mac'])  # guard the legacy precheck's broader splitlines alphabet
        fluka_text_scope(answer['fluka_inp'], cards(answer['fluka_inp']))
        gf = {FIELDS[k]: answer[k] for k in ('geometry_gdml','beam_mac','physics_txt')}
        for code, files in [('geant4',gf), ('fluka',{'fluka/input.inp':answer['fluka_inp']})]:
            verdict = precheck.check(files, code=code)
            if not verdict['allowed']:
                errors.append(code+': '+str(verdict['reason']))
        assert answer['physics_txt'].strip() == 'QGSP_INCLXX_HP', '題面物理表'
        gdml_topology_scope(answer['geometry_gdml'])
        gps_direction_scope(answer['beam_mac'])
        lines = answer['fluka_inp'].splitlines()
        nonempty = [l for l in lines if l.strip() and not l.startswith('*')]
        keys = [k for k,_,_,_ in cards(answer['fluka_inp'])]
        for k in ('START','STOP','GEOBEGIN','GEOEND','DEFAULTS'):
            assert keys.count(k) == 1, 'unique '+k
        assert keys[-2:] == ['START','STOP'], 'START/STOP must be final cards'
        assert nonempty[-1].strip().upper() == 'STOP', 'final STOP record must contain only the STOP keyword and optional trailing spaces; no WHAT fields or SDUM'
        # Close the historical precheck's geometry-block bypass; no arbitrary
        # control card can hide between GEOBEGIN and GEOEND.
        b = next((i for i,l in enumerate(nonempty) if l[:10].strip()=='GEOBEGIN'), None)
        e = next((i for i,l in enumerate(nonempty) if l[:10].strip()=='GEOEND'), None)
        assert b is not None and e is not None, 'exact GEOBEGIN/GEOEND keyword required'
        geo = nonempty[b+1:e]
        ends = [i for i,l in enumerate(geo) if l.strip().upper()=='END']
        assert len(ends)==2 and ends[-1]==len(geo)-1 and ends[0]>1, 'geometry structure'
        assert re.fullmatch(r'\s*\d+\s+\d+.*',geo[0]), 'geometry header'
        for line in geo[1:ends[0]]:
            words = line.split()
            assert words and words[0] in ('SPH','RCC','RPP'), 'inline primitive only'
            assert re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,7}',words[1])
            assert len(words)=={'SPH':6,'RCC':9,'RPP':8}[words[0]], 'body field count'
            assert all(math.isfinite(float(v.replace('D','E'))) for v in words[2:])
        for line in geo[ends[0]+1:ends[1]]:
            assert re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,7}\s+\d+\s+[ A-Za-z0-9_+\-()|]+',line.strip()), 'region expression'
        for key, whats, sdum, line in cards(answer['fluka_inp']):
            assert len(line)<=80, 'card longer than 80 columns'
            if key=='DEFAULTS': assert sdum=='PRECISIO', '題面 DEFAULTS'
            if key=='RANDOMIZ':
                assert float(whats[0] or 1)==1 and 0<=float(whats[1] or 0)<9e8, 'RANDOMIZ requires WHAT(1) blank or 1 and WHAT(2) blank or 0 <= seed < 9e8 (FLUKA 4-5.2 manual 7.22.65)'
            if key in ('USRYIELD','USRBDX') and sdum!='&':
                assert 21 <= abs(float(whats[2])) <= 99, 'scoring output unit'
            if key=='GEOBEGIN': assert sdum=='COMBNAME', 'COMBNAME required'
        coverage = [key for key,_,_,_ in cards(answer['fluka_inp']) if key not in COVERED]
        coverage += ['PHYSICS '+s for k,w,s,l in cards(answer['fluka_inp']) if k=='PHYSICS' and s!='COALESCE']
        coverage += fluka_field_coverage(cards(answer['fluka_inp']))
        if errors:
            return {'safe':False,'errors':errors,'coverage_issues':coverage}
        return {'safe':True,'errors':[],'coverage_issues':coverage,'options':opts}
    except (AssertionError,KeyError,ValueError,TypeError,IndexError) as error:
        return {'safe':False,'errors':errors+[str(error) or 'invalid output/host options'],'coverage_issues':[]}


def tsv(path):
    pairs = [line.split('\t',1) for line in path.read_text().splitlines()]
    assert all(len(p)==2 for p in pairs) and len(dict(pairs))==len(pairs), 'invalid TSV'
    return dict(pairs)


def read_g4(directory, answer, problem, execution=None):
    mat,geo,beam=problem["material"],problem["geometry"],problem["beam"]
    result, evidence = {}, {}
    run = directory/'output'
    d = tsv(run/'m1-readback.tsv')
    count = 0 if execution is None else execution['geant4']['events']
    assert 'geant4-11-04-patch-03' in d['geant4_version'] and d['requested_events']==str(count)
    assert d['macro_sha256']==hashlib.sha256(answer['beam_mac'].encode()).hexdigest()
    unbits = lambda k: struct.unpack('>d',bytes.fromhex(d[k]))[0]
    if execution is None:
        assert 'TRANSPORT_COMPLETED 0' in (directory/'driver.log').read_text()
    else:
        from . import transport_scope
        assert re.findall(r'(?m)^TRANSPORT_COMPLETED (\d+)\s*$', (directory/'driver.log').read_text()) == [str(count)], 'Geant4 completion count'
        assert d['seed'] == str(execution['geant4']['seed']), 'Geant4 seed readback'
        evidence['execution'] = transport_scope.bind_runtime(directory, 'geant4', execution, answer)
    plan = json.loads((directory.parent/'plan.json').read_text())
    runtime = json.loads((directory.parent/'runtime.json').read_text())
    assert runtime['plan_sha256'] == sha(directory.parent/'plan.json'), 'native plan identity'
    for key in ('geometry_gdml','beam_mac','physics_txt'):
        name = FIELDS[key]; digest = hashlib.sha256(answer[key].encode()).hexdigest()
        assert plan['files'][name] == runtime['geant4']['input_hashes'][name] == digest, 'native input binding: '+name
        assert sha(directory.parent/name) == digest, 'retained input changed: '+name
    topology = gdml_topology_scope(answer['geometry_gdml'])
    angular_scope = gps_direction_scope(answer['beam_mac'])
    expected_geometry={'density_g_cm3_bits':mat['density_g_cm3'],'target_radius_cm_bits':geo['radius_cm'],
                       'volume_mm3_bits':math.pi*geo['radius_cm']**2*geo['length_cm']*1000,
                       'target_length_cm_bits':geo['length_cm'],'world_radius_cm_bits':geo['world_radius_cm'],
                       **{'target_'+a+'_cm_bits':v for a,v in zip('xyz',geo['centre_cm'])},
                       **{'target_axis_'+a+'_bits':v for a,v in zip('xyz',geo['axis'])},
                       **{'isotope_'+str(mat['Z'])+'_'+str(i['A'])+'_fraction_bits':i['atom_fraction'] for i in mat['isotopes']}}
    expected_isotopes={k for k in expected_geometry if k.startswith('isotope_')}
    assert {k for k in d if re.fullmatch(r'isotope_\d+_\d+_fraction_bits',k)}==expected_isotopes, 'question isotope inventory'
    masses=mass_readback.geant4(d,unbits,problem)
    assert masses['conforms'], 'question isotope molar masses / atom number density'
    assert d['material']==problem['interface']['g4_material'], 'target material binding'
    for key,want in expected_geometry.items():
        assert close(unbits(key),want,1e-9), 'question geometry mismatch: '+key
    assert d['world_material']==problem['interface']['g4_world_material'] and d['target_shape']=='cylinder'
    gps={}
    for line in macro_lines(answer['beam_mac']):
        words=line.split()
        if words:gps[words[0]]=words[1:]
    assert gps.get('/gps/particle')==[beam['particle']]
    assert gps.get('/gps/pos/type')==['Point'] and gps.get('/gps/ene/type')==['Mono']
    energy=gps['/gps/ene/mono'];factor={'MeV':1.,'GeV':1000.,'keV':.001}[energy[1]]
    assert close(float(energy[0])*factor,beam['energy_MeV'],1e-9), 'question beam energy'
    pos=gps['/gps/pos/centre'];factor={'cm':1.,'mm':.1,'m':100.}[pos[3]]
    assert all(close(float(v)*factor,w,1e-9) for v,w in zip(pos[:3],beam['position_cm']))
    direction=[float(v) for v in gps['/gps/direction']]
    assert len(direction)==3 and sum(v*v for v in direction)>0
    norm=math.sqrt(sum(v*v for v in direction))
    assert all(close(v/norm,w,1e-9) for v,w in zip(direction,beam['direction']))
    incl = tsv(run/'incl-physics.tsv')
    assert incl['configuration_verified']=='1' and incl['phase']=='idle_before_first_BeamOn'
    assert incl['schema_version']=='t060-incl/1' and incl['command_status']=='0'
    result['P01'] = d['physics'] if any('HadronPhysicsINCLXX' in v for k,v in d.items() if k.startswith('physics_constructor_type_')) else None
    result['P02'] = incl['mode']
    result['P03'] = {'0':True,'1':False}.get(incl.get('after_cluster'))
    processes = [v for k,v in d.items() if re.fullmatch(r'neutron_process_\d+',k)]
    result['P23'] = any('killer' in p.lower() or 'cut' in p.lower() for p in processes) if processes else None
    result['P24'] = {p+'_mm':unbits('target_'+p+'_cut_mm_bits') for p in ('gamma','electron','positron')}
    result['P26'] = unbits('target_proton_cut_mm_bits')
    dataset=Path(d.get('neutron_data_path','').replace(chr(92),'/')).name
    data={'dataset':dataset} if dataset.startswith('G4NDL') else None
    hp=hp_coverage(d,unbits,mat)
    result['P28'] = data['dataset'] if data else None
    result['P29'] = unbits('material_temperature_K_bits')
    result['P30'] = 'HP_'+mat['symbol']+'_0_20MeV_dispatch' if hp['complete'] else None
    for key, process in [('P31','photonuclear'),('P32','electronuclear')]:
        state, error = g4read.process_switch_detail(d,process)
        if state:
            result[key] = (state['particles']['gamma']['enabled'] if key=='P31' else
                           {p:v['enabled'] for p,v in state['particles'].items()})
        else: result[key] = None
        evidence[key] = {'detail':state,'unread':error}
    evidence.update(identity=d, incl=incl, hp=hp, question_mass=masses, fixed_host_input_scope={'topology':topology,'angular':angular_scope},
                    raw_sha256={p.name:sha(p) for p in run.iterdir() if p.is_file()})
    return result, evidence


def read_fluka(directory, answer, problem, execution=None):
    mat,beam=problem["material"],problem["beam"]
    plan=json.loads((directory.parent/'plan.json').read_text())
    runtime=json.loads((directory.parent/'runtime.json').read_text())
    assert runtime['plan_sha256']==sha(directory.parent/'plan.json'), 'native plan identity'
    original=answer['fluka_inp']
    fluka_text_scope(original, cards(original))
    if execution is None:
        executed='\n'.join(line[:10]+'1.0'.rjust(10)+line[20:] if line[:10].strip().upper()=='START' else line
                           for line in original.splitlines())+'\n'
        count = 1
    else:
        from . import transport_scope
        count = execution['fluka']['histories']
        executed = transport_scope.fluka_input(original, count, execution['fluka']['seed'])
        execution_evidence = transport_scope.bind_runtime(directory, 'fluka', execution, answer)
    for name,text in [('fluka/original.inp',original),('fluka/input.inp',executed)]:
        digest=hashlib.sha256(text.encode()).hexdigest()
        assert plan['files'][name]==runtime['fluka']['input_hashes'][name]==digest, 'native input binding: '+name
        assert sha(directory.parent/name)==digest, 'retained input changed: '+name
    raw = (directory/'input001.out').read_text(errors='replace')
    assert 'FLUKA Version 4-5.2' in raw or re.search(r'FLUKA Version\s+4-5\.2\s',raw)
    assert 'End of the output associated with the input' in raw and 'End of FLUKA v4 run' in raw
    counts = re.findall(r'Total number of primaries run:\s*(\d+)',raw)
    assert counts==[str(count)], 'FLUKA history count not exactly one' if execution is None else 'FLUKA history count differs from execution contract'
    if execution is not None:
        weights = re.findall(r'Total number of primaries run:\s*\d+\s+for a weight of:\s*('+NUM+')', raw)
        assert len(weights) == 1 and close(float(weights[0]), count, 6e-7), 'FLUKA completed primary weight'
        echoes = [line for line in raw.splitlines() if 'Next control card' in line and re.search(r'\bRANDOMIZ\b', line)]
        assert len(echoes) == 1, 'FLUKA RANDOMIZ echo inventory'
        values = echoes[0].split('RANDOMIZ', 1)[1].split()
        assert len(values) == 6 and float(values[0]) == 1 and all(float(v) == 0 for v in values[2:]), 'FLUKA RANDOMIZ echo controls'
        assert close(float(values[1]), execution['fluka']['seed'], 6e-5), 'FLUKA seed echo at printed resolution'
        execution_evidence['exact_seed'] = transport_scope.fluka_seed(raw, execution['fluka']['seed'])
        starts = [line for line in raw.splitlines() if 'Next control card' in line and re.search(r'\bSTART\b', line)]
        assert len(starts) == 1 and close(float(starts[0].split('START', 1)[1].split()[0]), count, 6e-5), 'FLUKA START echo count'
    beam_match=re.search(r'Average beam kinetic energy:\s+('+NUM+')',raw)
    assert beam_match and close(float(beam_match[1])*1000,beam['energy_MeV'],6e-7), 'question beam energy'
    assert re.search(r'Beam particle:\s+'+re.escape(beam['particle'].upper())+r'\b',raw), 'question beam particle'
    for label,wanted in [('Beam hit position',beam['position_cm']),('Beam direction cosines',beam['direction'])]:
        m=re.search(label+r'\s*:\s*('+NUM+r')\s+('+NUM+r')\s+('+NUM+r')',raw)
        assert m and all(close(float(v),w,6e-7) for v,w in zip(m.groups(),wanted)), 'question '+label
    for label in ('Momentum deviation at FWHM', 'Beam spot FWHM X-width',
                  'Beam spot FWHM Y-width', 'Beam FWHM angular divergence'):
        values=re.findall(re.escape(label)+r'[^:\n]*:\s*('+NUM+r')',raw)
        assert len(values)==1 and float(values[0].replace('D','E'))==0, 'question zero '+label
    inp = cards(answer['fluka_inp'])
    echoed = [l for l in raw.splitlines() if 'Next control card' in l]
    defaults = [l.split()[-1] for l in echoed if re.search(r'\bDEFAULTS\b',l)]
    result = {'P01':defaults[0] if defaults==['PRECISIO'] else None}
    physics = [(w,s) for k,w,s,l in inp if k=='PHYSICS' and s=='COALESCE']
    pe = [l.split('PHYSICS',1)[1].split() for l in echoed if re.search(r'\bPHYSICS\b',l) and l.split()[-1]=='COALESCE']
    cluster = False if not physics and not pe and result['P01']=='PRECISIO' else None
    if physics and len(physics)==len(pe) and result['P01']=='PRECISIO':
        cluster=False
        for (w,s), echo in zip(physics,pe):
            try:
                values=[float(x or 0) for x in w]
                assert s=='COALESCE' and echo[-1]=='COALESCE' and len(echo)==7
                assert all(close(a,float(b),6e-4) for a,b in zip(values,echo[:6])) and values[1:]==[0.]*5
                if values[0]: cluster=values[0]>0
            except (AssertionError,ValueError): cluster=None;break
    result['P03']=cluster
    cuts={p:float(v.replace('D','E'))*1000 for p,v in re.findall(r'Cut-off kinetic energy for (PROTON|NEUTRON)\s+transport:\s+('+NUM+') GeV',raw)}
    result['P22']=cuts.get('PROTON');result['P23']=cuts.get('NEUTRON')
    em=electromagnetic_cuts(raw)
    for region,setting in problem['interface']['fluka_exterior_materials'].items():
        outside=(em or {}).get('regions',{}).get(region,{})
        assert outside.get('fluka_material_id')==setting['id'] and outside.get('fluka_material')==setting['name'], 'question exterior material: '+region
    target=(em or {}).get('regions',{}).get(problem['interface']['target_region'],{})
    binding=fluka_target(raw,target,problem)
    masses=mass_readback.fluka(raw,binding,problem)
    assert masses['conforms'], 'question isotope molar masses / atom number density'
    material=target.get('emf_material')
    result['P24']=(em or {}).get('production_by_material',{}).get(material)
    result['P25']={k:target[k] for k in ('prompt','delayed')} if all(k in target for k in ('prompt','delayed')) else None
    rows=binding['pointwise_rows']
    data_ok=binding['pointwise_complete']
    result['P28']=rows[0]['dataset'] if data_ok else None
    result['P29']=rows[0]['temperature_K'] if data_ok else None
    result['P30']='pointwise' if data_ok else None
    mid=target.get('fluka_material_id')
    details={}
    for pid,label,sdum in [('P31','Gamma',''),('P32','Electron','ELECTNUC')]:
        activated=[(int(m),int(c)) for m,c in re.findall(label+r'\s+Photo-nuclear int\. activated for media #\s*(\d+)\s+with code\s+(\d+)',raw)]
        pcs=[(w,s) for k,w,s,l in inp if k=='PHOTONUC' and (s=='ELECTNUC')==(pid=='P32')]
        # An activation line is actual evidence; an absent line is usable only
        # with a complete supported deck, known defaults and no positive card.
        target_codes=[c for m,c in activated if m==mid]
        if target_codes: value=True if target_codes==[1111] else None
        elif not pcs and result['P01']=='PRECISIO': value=False
        elif pcs and all(float(w[0] or 0)<=0 for w,s in pcs) and result['P01']=='PRECISIO': value=False
        else: value=None
        result[pid]=value;details[pid]={'activations':activated,'target_fluka_material_id':mid,
            'unsupported_target_codes':target_codes if target_codes and target_codes!=[1111] else []}
    evidence = {'raw_sha256':sha(directory/'input001.out'),'em_tables':em,'neutron_rows':rows,'question_mass':masses,
                'material_binding':binding,'switch_detail':details,'cluster_support':'accepted PHYSICS echo + complete run + 4-5.2/PRECISIO rule'}
    if execution is not None:
        evidence['execution'] = execution_evidence
    return result, evidence


def close(a,b,rtol=6e-4,atol=1e-10):
    if a is None or b is None: return False
    if type(a) is bool or type(b) is bool: return type(a) is type(b) and a==b
    if isinstance(a,dict) and isinstance(b,dict):
        return set(a)==set(b) and all(close(a[k],b[k],rtol,atol) for k in a)
    if type(a) in (int,float) and type(b) in (int,float):
        return math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=rtol,abs_tol=atol)
    return type(a) is type(b) and a==b


def explicit(answer,pid,side):
    if side=='geant4':
        opt=json.loads(answer['host_options_json']);gd=ET.fromstring(answer['geometry_gdml'])
        if pid=='P01':return True
        if pid in ('P02','P03'):return 'incl_physics' in opt
        if pid=='P24':return 'range_cut_mm' in opt
        if pid=='P29':return any(n.tag.rsplit('}',1)[-1]=='T' for n in gd.iter())
        if pid=='P31':return 'gamma_nuclear' in opt
        return False
    cs=cards(answer['fluka_inp']);keys={k for k,w,s,l in cs}
    if pid=='P01':return True
    if pid=='P03':return any(k=='PHYSICS' and s=='COALESCE' and float((w[0] or '0').upper().replace('D','E'))!=0 for k,w,s,l in cs)
    if pid in ('P22','P23'):
        target={'P22':1,'P23':8}[pid]
        def particle(x):
            if x.upper() in ('PROTON','NEUTRON'):return {'PROTON':1,'NEUTRON':8}[x.upper()]
            return int(float(x))
        for k,w,s,l in cs:
            if k=='PART-THR' and w[0]:
                try:
                    lo=particle(w[1] or '1');hi=particle(w[2]) if w[2] else lo
                    if target in range(lo,hi+1,int(float(w[3] or 1))):return True
                except (ValueError,OverflowError):pass
        return False
    if pid=='P24':return any(k=='EMFCUT' and s=='PROD-CUT' for k,w,s,l in cs)
    if pid=='P25':return any(k=='EMFCUT' and s!='PROD-CUT' for k,w,s,l in cs)
    if pid=='P28':return any(k=='LOW-PWXS' and s for k,w,s,l in cs)
    if pid=='P29':
        try:return any(k=='LOW-PWXS' and w[2] and float(w[2].replace('D','E'))>0 for k,w,s,l in cs)
        except (ValueError,OverflowError):return None  # Named coverage error from validate().
    if pid=='P30':return 'LOW-PWXS' in keys
    if pid in ('P31','P32'):return any(k=='PHOTONUC' and (s=='ELECTNUC')==(pid=='P32') for k,w,s,l in cs)
    return False


def assess(answer, runtime, root, problem, generation_ok=True, execution=None):
    safety=validate(answer);actual={};evidence={};errors=[]
    for side,reader in [('geant4',read_g4),('fluka',read_fluka)]:
        try:
            assert runtime.get(side,{}).get('exit')==0, 'native run failed or not submitted'
            actual[side],evidence[side]=reader(Path(root)/side,answer,problem,execution)
        except Exception as error:
            actual[side]={};errors.append(side+': '+repr(error))
    try:answers=json.loads(answer['answers_json']);differences=json.loads(answer['differences_json'])
    except (ValueError,KeyError):answers={};differences={}
    table=[]
    for definition in CHECKLIST['per_input']:
        pid=definition['id'];row=dict(definition,sides={})
        for side in definition['sides']:
            value=actual.get(side,{}).get(pid);declared=answers.get(pid,{}).get(side) if isinstance(answers.get(pid),dict) else None
            confirmed=close(value,declared,0 if pid=='P29' else 1e-9 if side=='geant4' else 6e-4,0.051 if pid=='P29' else 0.0)
            row['sides'][side]={'effective':value,'answer':declared,'explicit':explicit(answer,pid,side) if safety['safe'] else None,
                                'confirmed':confirmed,'status':'confirmed' if confirmed else 'no_readback' if value is None else 'answer_missing_or_mismatch'}
        row['aligned']=None;row['registered_difference']=None
        if definition['align'] and set(row['sides'])=={'geant4','fluka'}:
            a=row['sides']['geant4']['effective'];b=row['sides']['fluka']['effective']
            if pid=='P32' and isinstance(a,dict): row['aligned']=all(close(v,b) for v in a.values())
            else:row['aligned']=close(a,b,0,0.051 if pid=='P29' else 0)
            reason=differences.get(pid)
            if row['aligned'] is False and isinstance(reason,str) and reason.strip():row['registered_difference']=reason
        row['passed']=all(v['confirmed'] for v in row['sides'].values()) and (not definition['align'] or row['aligned'] or bool(row['registered_difference']))
        table.append(row)
    intrinsic=[{'id':pid,'status':'registered','confirmed':False,'basis':'T058 program-specific boundary; no per-input confirmation claimed'} for pid in CHECKLIST['intrinsic']]
    passed=bool(generation_ok and safety['safe'] and not safety['coverage_issues'] and not errors and all(r['passed'] for r in table))
    inp=cards(answer.get('fluka_inp',''))
    fluka_evidence=evidence.get('fluka',{})
    target=fluka_evidence.get('em_tables',{}).get('regions',{}).get(problem['interface']['target_region'],{})
    # Only consume raw output whose hash/binding has passed the existing reader.
    raw=(Path(root)/'fluka/input001.out').read_text(errors='replace') if fluka_evidence else ''
    prerequisites={'fluka_fudgem':fudgem_echo_prerequisite(inp,raw,target),
                   'fluka_blank_fields':blank_input_prerequisites(inp)}
    passed=passed and prerequisites['fluka_fudgem']['input_scope_satisfied']
    try:
        prerequisites['fluka_text']=fluka_text_scope(answer.get('fluka_inp',''), cards(answer.get('fluka_inp','')))
    except (ValueError, TypeError) as error:
        prerequisites['fluka_text']={'basis':'input_text_only','native_readback':False,
            'input_scope_satisfied':False,'error':str(error),
            'independent_review':'T067 离线修复；作为正式依据前待合并审查'}
    return {'input_prerequisites':prerequisites,'passed':passed,'generation_ok':generation_ok,'safety':safety,'runtime':runtime,'readback_errors':errors,
            'table':table,'intrinsic':intrinsic,'evidence':evidence,'new_spectrum_estimate':False}


def feedback(record):
    short={k:v for k,v in record.items() if k not in ('evidence',)}
    return ('输入范围前提另列于 input_prerequisites；FUDGEM 核对输入文本加原生输入回显，无生效值读回；留空角落只查输入文本。\n'
            '以下是本轮 kernel 记录；effective 只来自原生读回，null 表示没有证据。请修正格式/运行失败。'
            '对逐份清单每项，逐侧决定你要取什么，并在 answers_json 中按 P编号 -> geant4/fluka 填数值、布尔或记录中的对象结构。'
            '可以选择保留当前读回值，但必须明确回答。可对齐项若选择保留不同，在 differences_json 的该 P编号写理由。'
            '程序固有登记项只告知，不要求回答。内部 ID/类名不填。请重新交出全部七字段与完整输入，不只交补丁。\n'
            +json.dumps(short,ensure_ascii=False,indent=2))
