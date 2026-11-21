"""Read formatted USRYIELD output (FLUKA 4-5.2 manual, section 9.5.9).

One file holds the detectors of one cycle. Each detector block gives the binning of the first
quantity (many intervals), the single interval of the second quantity, and a vector of plain
double-differential yields per primary. Anything that does not match the documented layout raises.
"""
import re
import math

TOTAL = re.compile(r"Total number of particles followed\s+(\d+), for a total weight of\s+(\S+)")
HEAD = re.compile(r'Yield n\.\s+(\d+)\s+"([^"]*)"\s*, generali[sz]ed particle n\.\s+(-?\d+), '
                  r"from region n\.\s+(-?\d+) to region n\.\s+(-?\d+)")
FIRST = re.compile(r"(linear|logar\.) 1st variable binning from\s+(\S+) to\s+(\S+)\s+(\d+) bins")
SECOND = re.compile(r"2nd variable ranges from\s+(\S+) to\s+(\S+)")
KIND = re.compile(r"(1st|2nd) variable is: (.+?)\s*$", re.M)
DATA = "Data follow in a vector A(ie)"
NORM = re.compile(r"user normalization:\s+(\S+), adopted cross section \(if any\):\s+(\S+) mb")


def read(text):
    total = TOTAL.search(text)
    if not total:
        raise ValueError("No primary count in USRYIELD output")
    heads = list(HEAD.finditer(text))
    if not heads:
        raise ValueError("No USRYIELD detector in output")
    detectors = []
    for i, head in enumerate(heads):
        block = text[head.end():heads[i + 1].start() if i + 1 < len(heads) else len(text)]
        first, second = FIRST.search(block), SECOND.search(block)
        if not first or not second or DATA not in block:
            raise ValueError("Detector %s does not follow the documented layout" % head[2].strip())
        kinds = dict(KIND.findall(block))
        normalization = NORM.search(block)
        lines = block.split(DATA, 1)[1].splitlines()[1:]
        values = [float(token) for line in lines for token in line.split()
                  if line.strip() and line.strip() != "1"]
        bins = int(first[4])
        if len(values) != bins:
            raise ValueError("Detector %s has %d values for %d bins" % (head[2].strip(), len(values), bins))
        if not all(math.isfinite(v) and v >= 0 for v in values):
            raise ValueError("Nonfinite or negative plain yield")
        detectors.append({"number": int(head[1]), "name": head[2].strip(), "particle": int(head[3]),
                          "normalization": float(normalization[1]) if normalization else None,
                          "cross_section_mb": float(normalization[2]) if normalization else None,
                          "from_region": int(head[4]), "to_region": int(head[5]),
                          "first": {"scale": first[1], "low": float(first[2]), "high": float(first[3]),
                                    "bins": bins, "quantity": kinds.get("1st")},
                          "second": {"low": float(second[1]), "high": float(second[2]), "quantity": kinds.get("2nd")},
                          "values": values})
    return {"primaries": int(total[1]), "weight": float(total[2]), "detectors": detectors}
