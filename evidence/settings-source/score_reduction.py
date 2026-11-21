"""Fixed host mapping of unbinned Geant4 crossings to question cells.

This module is not executed on new events in T065. Tests use synthetic records.
Unit-weight analog histories only: the registered host stores no track weight.
"""
import math


def cell(pdg, energy_MeV, direction, scoring):
    if pdg != 2112: return None
    if not math.isfinite(energy_MeV) or not all(math.isfinite(x) for x in direction):
        raise ValueError('nonfinite scoring record')
    if len(direction)!=3 or not math.isclose(sum(x*x for x in direction),1,rel_tol=1e-9):
        raise ValueError('non-unit direction')
    axis=scoring['surface']['axis']
    theta=math.degrees(math.acos(max(-1,min(1,sum(x*y for x,y in zip(direction,axis))))))
    energies=scoring['energy_edges_MeV']; angles=scoring['angle_edges_deg']
    ei=next((i for i,(lo,hi) in enumerate(zip(energies,energies[1:])) if lo<=energy_MeV<hi),None)
    ai=next((i for i,(lo,hi) in enumerate(zip(angles,angles[1:])) if lo<=theta<hi or i==len(angles)-2 and theta==hi),None)
    return (ei,ai) if ei is not None and ai is not None else None


def yield_factor(primary_count, energy_pair_MeV, angle_pair_deg):
    if primary_count<=0: raise ValueError('no spectrum from zero primaries')
    lo,hi=energy_pair_MeV; a,b=angle_pair_deg
    width=hi-lo; solid=2*math.pi*(math.cos(math.radians(a))-math.cos(math.radians(b)))
    if width<=0 or solid<=0: raise ValueError('invalid scoring cell')
    return 1/(primary_count*width*solid)
