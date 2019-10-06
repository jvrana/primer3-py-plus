import re
from typing import Tuple

import primer3

from primer3plus.utils import reverse_complement as rc


class PrimerResult:
    def __init__(
        self, name, start, anneal, overhang, strand, meta, thermo_settings=None
    ):
        self.name = name
        self.start = start
        self.anneal = anneal
        self.overhang = overhang
        self.strand = strand
        self._tm = None
        self._thermo = None
        self._thermo_tm = None
        if thermo_settings is None:
            thermo_settings = {}
        self._thermo_settings = thermo_settings
        self.meta = meta

    @property
    def tm(self):
        if self._tm is None:
            self._tm = primer3.calcTm(self.anneal)
        return self._tm

    def thermo(self):
        settings = self._thermo_settings
        if self._thermo is None:
            if len(self.sequence) > 60:
                warning = "sequence length greater than 60. Thermo results are limited to 60bp."
            else:
                warning = ""
            self._thermo = {
                "hairpin": primer3.calcHairpin(self._safe_sequence, **settings),
                "homodimer": primer3.calcHomodimer(self._safe_sequence, **settings),
                "annealing": primer3.calcHeterodimer(
                    self.anneal, rc(self.anneal), **settings
                ),
                "sequence": primer3.calcHeterodimer(
                    self._safe_sequence, rc(self._safe_sequence), **settings
                ),
                "warning": warning,
            }
        return self._thermo

    def thermo_tm(self):
        x = {k: v.tm for k, v in self.thermo().items() if not isinstance(v, str)}
        x["warning"] = self.thermo()["warning"]
        return x

    @property
    def anneal_length(self):
        return len(self.sequence)

    @property
    def sequence(self):
        return self.anneal + self.overhang

    @property
    def _safe_sequence(self):
        return self.sequence[-60:]

    def _thermowarning(self):
        if len(self.sequence) > 60:
            return (
                "sequence length greater than 60. Thermo results are limited to 60bp."
            )
        return ""

    @property
    def top_strand_slice(self):
        if self.strand == 1:
            return self.start, self.start + self.anneal_length
        elif self.strand == -1:
            return self.start - self.anneal_length, self.start


class PairResult:
    def __init__(self, p1, p2, size, meta):
        self.p1 = p1
        self.p2 = p2
        self.size = size
        self.meta = meta

    def thermo(self):
        return {
            "heterodimer": primer3.calcHeterodimer(
                self.p1._safe_sequence, self.p2._safe_sequence
            ),
            "left": self.p1.thermo(),
            "right": self.p2.thermo(),
        }

    def thermo_tm(self):
        return {
            "heterodimer": primer3.calcHeterodimerTm(
                self.p1._safe_sequence, self.p2._safe_sequence
            ),
            "left": self.p1.thermo_tm(),
            "right": self.p2.thermo_tm(),
        }


def parse_primer3_results(results_dict: dict) -> Tuple[dict, dict]:
    """Parse the primer3 results. To JSON format.

    :param results_dict: the raw results dictionary
    :return: parse dictionary
    """
    num_pairs = results_dict["PRIMER_PAIR_NUM_RETURNED"]
    num_left = results_dict["PRIMER_LEFT_NUM_RETURNED"]
    num_right = results_dict["PRIMER_RIGHT_NUM_RETURNED"]

    pairs = {}
    other = {}
    for i in range(max([num_pairs, num_left, num_right])):
        pairs.setdefault(i, {})
    key_pattern = r"PRIMER_(?P<label>[a-zA-Z]+)_(?P<pair_id>\d+)_(?P<key>.+)"
    location_pattern = r"PRIMER_(?P<label>[a-zA-Z]+)_(?P<pair_id>\d+)\s*$"
    for k in results_dict:
        m = re.match(key_pattern, k)
        loc_m = re.match(location_pattern, k)
        if m:
            groupdict = m.groupdict()
            pair_id = int(groupdict["pair_id"])
            label = groupdict["label"]
            key = groupdict["key"]
            pairdict = pairs[pair_id]
            pairdict.setdefault(label, {})
            pairdict[label][key] = results_dict[k]
        elif loc_m:
            groupdict = loc_m.groupdict()
            pair_id = int(groupdict["pair_id"])
            label = groupdict["label"]
            pairdict = pairs[pair_id]
            pairdict.setdefault(label, {})
            pairdict[label]["location"] = results_dict[k]
        else:
            other[k] = results_dict[k]
    return pairs, other


def to_pair_result(data):
    left_data = data["LEFT"]
    right_data = data["RIGHT"]
    left = PrimerResult(
        name="",
        start=left_data["location"][0],
        anneal=left_data["SEQUENCE"],
        overhang=left_data["OVERHANG"],
        strand=1,
        meta=left_data,
    )
    right = PrimerResult(
        name="",
        start=right_data["location"][0],
        anneal=right_data["SEQUENCE"],
        overhang=left_data["OVERHANG"],
        strand=-1,
        meta=right_data,
    )
    pair = PairResult(left, right, data["PAIR"]["PRODUCT_SIZE"], meta=data["PAIR"])
    return pair
