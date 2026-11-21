"""Fixed T056 counting contract and standard-library helpers."""
import hashlib
import json
import math
from pathlib import Path

EDGES = [15., 35.18064, 41.25585, 48.38017, 56.73476, 66.53208,
         78.02125, 91.49445, 107.29429, 200.]
ANGLES = [0., 45., 90., 135., 180.]
SOLID = [2*math.pi*(math.cos(math.radians(a))-math.cos(math.radians(b)))
         for a, b in zip(ANGLES[:-1], ANGLES[1:])]
LIGHT = {(2, 1): ('DEUTERON', -3), (3, 1): ('TRITON', -4),
         (3, 2): ('3-HELIUM', -5), (4, 2): ('4-HELIUM', -6)}
SPECIES = [dict(name='n', A=1, Z=0, particle='NEUTRON', pid=8),
           dict(name='p', A=1, Z=1, particle='PROTON', pid=1)]
for A in range(2, 9):
    for Z in range(1, A+1):
        particle, pid = LIGHT.get((A, Z), ('ALL-PART', 201))
        SPECIES.append(dict(name=f'A{A}Z{Z}', A=A, Z=Z, particle=particle, pid=pid,
                            filter=None if (A, Z) in LIGHT else -(A*100000+Z*100)))
AUDIT = [dict(name=f'checkA{A}', A=A, Z=0, particle='ALL-PART', pid=201,
              filter=-A*100000) for A in range(2, 9)]
AUDIT += [dict(name=f'checkA{A}Z{Z}', A=A, Z=Z, particle='ALL-PART', pid=201,
               filter=-(A*100000+Z*100)) for A, Z in LIGHT]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, allow_nan=False, separators=(',', ':'))+'\n')


def pdg_species(pdg):
    if pdg == 2112:
        return (1, 0)
    if pdg == 2212:
        return (1, 1)
    if 1000000000 <= pdg < 1010000000:
        return ((pdg//10) % 1000, (pdg//10000) % 1000)
    return None


def energy_bin(energy, A):
    # Compare total energies, exactly as FLUKA cards do; no Jacobian or division rounding.
    return next((i for i, (lo, hi) in enumerate(zip(EDGES[:-1], EDGES[1:]))
                 if A*lo <= energy < A*hi), None)


def angle_bin(theta):
    return min(3, int(theta/45)) if 0 <= theta <= 180 else None
