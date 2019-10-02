"""To run designs, you copy the underlying parameters and change certain
values.

There should be a number of presets available to run certain design
tasks  There should also be a relaxation parameter (a subclass of
Design?)  Async option  Design with overhangs
"""
from __future__ import annotations

import itertools
import re
import webbrowser
from collections import Counter
from functools import wraps
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union

import primer3

from .interfaces import AllParameters
from .interfaces import ParameterAccessor
from .results import parse_primer3_results
from primer3plus.constants import DOCURL
from primer3plus.log import logger
from primer3plus.params import default_boulderio


def _summarize_reasons(reasons):
    reason_dict = {}
    for reason in reasons:
        for k, v in reason.items():
            if "EXPLAIN" in k:
                for m in re.finditer(r"\s*([\w\s\-]+)\s+(\d+)", v):
                    reason_token = m.group(1)
                    num = int(m.group(2))
                    reason_dict.setdefault(k, Counter())[reason_token] += num
    return {k: dict(v) for k, v in reason_dict.items()}


def combine_results(res):
    """Combine and sort results. Combine all explainations.

    :param results: :type results: :return: :rtype:
    """
    all_pairs = []
    all_reasons = []
    for r in res:
        all_pairs += list(r[0].values())
        all_reasons.append(r[1])
    sorted_pairs = sorted(all_pairs, key=lambda x: x["PAIR"]["PENALTY"])
    explain = _summarize_reasons(all_reasons)
    return sorted_pairs, explain


class DesignPresets:
    """
    Interface for setting design parameters. This is typically accessed from
    a :class:`Design <primer3plus.Design>` instance's
    :meth:`Design <primer3plus.Design.set>` method. As in:

    .. code-block::

        design = Design()
        design.presets.left_sequence("AGGGAGATAGATA")
        design.run()
    """

    def __init__(self, design):
        """
        Initializes a new interface from a :class:`~primer3plus.design.Design`.

        :param design: The design
        """
        self._design = design

    def _interval_from_sequences(
        self, template: str, target: str
    ) -> Union[None, Tuple[int, int]]:
        if isinstance(target, str):
            matches = self._get_index_of_match(template, target)
            if not matches:
                print("Target not in template")
                return None
            if len(matches) > 1:
                print("More than one target found")
                return None
            return matches[0]

    @staticmethod
    def _get_index_of_match(template: str, sequence: str) -> List[Tuple[int, int]]:
        matches = []
        for m in re.finditer(sequence, template, re.IGNORECASE):
            matches.append((m.start(0), m.end(0)))
        return matches

    def update(self, update: Dict[str, Any]):
        """Update an arbitrary parameter"""
        self._design.params.update(update)
        return self

    def task(self, task: str) -> DesignPresets:
        """This tag tells primer3 what task to perform.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_TASK

        :param task: the task name
        :return self
        """
        self.update({"PRIMER_TASK": task})
        return self

    def as_cloning_task(self) -> DesignPresets:
        """
        Set the design as a cloning task.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_TASK

        :return: self
        """
        return self.task("pick_cloning_primers")

    def as_generic_task(self) -> DesignPresets:
        """
        Set the design as a generic task.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_TASK

        :return: self
        """
        return self.task("generic")

    def template(self, template: str) -> DesignPresets:
        """
        Set the template sequence for the design. This sets the 'SEQUENCE_TEMPLATE'
        parameter.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_TEMPLATE

        :param template: the template sequence
        :return: self
        """
        self.update({"SEQUENCE_TEMPLATE": template})
        return self

    # TODO: set_iterations, set_num_return, set_force_return, set_gradient
    def primer_num_return(self, n: int) -> DesignPresets:
        """
        Set the number of primers to return for the design task.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_NUM_RETURN

        :param n: number of primers to return
        :return: self
        """
        return self.update({"PRIMER_NUM_RETURN": n})

    def product_size(
        self, interval: Union[Tuple[int, int], List[Tuple[int, int]]], opt=None
    ) -> DesignPresets:
        """
        Set the product size. Optionally include the optimal size.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PRODUCT_SIZE_RANGE
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PRODUCT_OPT_SIZE

        :param interval: a tuple of <min>,<max> or a list of such tuples
        :param opt: optional product size as an int.
        :return: self
        """
        if isinstance(interval, tuple):
            interval = [interval]
        if opt is not None:
            return self.update(
                {"PRIMER_PRODUCT_SIZE_RANGE": interval, "PRIMER_PRODUCT_OPT_SIZE": opt}
            )
        return self.update({"PRIMER_PRODUCT_SIZE_RANGE": interval})

    def pair_region_list(
        self, region_list: List[Tuple[int, int, int, int]]
    ) -> DesignPresets:
        """
        The list of regions from which to design primers.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER_PAIR_OK_REGION_LIST

        :param region_list: list of regions
        :return: self
        """
        return self.update({"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": region_list})

    def left_sequence(self, primer: str) -> DesignPresets:
        """The sequence of a left primer to check and around which to design
        right primers and optional internal oligos. Must be a substring of
        SEQUENCE_TEMPLATE.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER

        :param primer: :type primer: :return: :rtype:
        """
        return self.update({"SEQUENCE_PRIMER": primer, "PRIMER_PICK_RIGHT_PRIMER": 1})

    def right_sequence(self, primer: str) -> DesignPresets:
        """The sequence of a right primer to check and around which to design
        left primers and optional internal oligos. Must be a substring of the
        reverse strand of SEQUENCE_TEMPLATE.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER_REVCOMP
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER

        :param primer: primer sequence
        :return: self
        """
        return self.update(
            {"SEQUENCE_PRIMER_REVCOMP": primer, "PRIMER_PICK_LEFT_PRIMER": 1}
        )

    def pick_left_only(self) -> DesignPresets:
        """
        Design only the left primer.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER

        :return: self
        """
        return self.update(
            {"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0}
        )

    def pick_right_only(self) -> DesignPresets:
        """
        Design only the right primer.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER

        :return: self
        """
        return self.update(
            {"PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1}
        )

    def internal_sequence(self, primer: str) -> DesignPresets:
        """The sequence of an internal oligo to check and around which to
        design left and right primers. Must be a substring of
        SEQUENCE_TEMPLATE.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_INTERNAL_OLIGO
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_INTERNAL_OLIGO


        :param primer: :type primer: :return: :rtype:
        """
        return self.update(
            {"SEQUENCE_INTERNAL_OLIGO": primer, "PRIMER_PICK_INTERNAL_OLIGO": 1}
        )

    def primers(self, p1: str, p2: str) -> DesignPresets:
        """
        Set the left and right primer sequences.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER_REVCOMP
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER

        :param p1:
        :param p2:
        :return:
        """
        if p1:
            self.left_sequence(p1)
        if p2:
            self.right_sequence(p2)
        return self

    def _parse_interval(
        self, interval: Union[str, Tuple[int, int], List[Tuple[int, int]]]
    ) -> List[Tuple[int, int]]:
        if isinstance(interval, str):
            interval = self._interval_from_sequences(
                self._design.params["SEQUENCE_TEMPLATE"], interval
            )
        if isinstance(interval, tuple):
            interval = [interval]
        return interval

    def included(
        self, interval: Union[str, Tuple[int, int], List[Tuple[int, int]]]
    ) -> DesignPresets:
        """
        Specify interval from which primers must be selected.
        A sub-region of the given sequence in which to pick primers. For
        example, often the first dozen or so bases of a sequence are vector,
        and should be excluded from consideration. The value for this parameter
        has the form.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_INCLUDED_REGION

        <start>,<length> where <start> is the index of the first base to consider, and
        <length> is the number of subsequent bases in the primer-picking region.

        :param interval: One of the following: the sequence of the target region,
                         a tuple of the interval of <start>,<length>, or a list of
                         tuples of <start>,<length>
        :return: self
        """
        return self.update({"SEQUENCE_INCLUDED_REGION": self._parse_interval(interval)})

    def target(
        self, interval: Union[str, Tuple[int, int], List[Tuple[int, int]]]
    ) -> DesignPresets:
        """
        Specify the interval that designed primers must flank.
        If one or more targets is specified then a legal primer pair must
        flank at least one of them. A target might be a simple sequence repeat
        site (for example a CA repeat) or a single-base-pair polymorphism, or
        an exon for resequencing. The value should be a space-separated list
        of <start>,<length> pairs where <start> is the index of the first base of a
        target,and <length> is its length.  See also PRIMER_INSIDE_PENALTY,
        PRIMER_OUTSIDE_PENALTY.
        PRIMER_TASK=pick_sequencing_primers. See PRIMER_TASK for more information.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_TEMPLATE
        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_TARGET


        :param interval: One of the following: the sequence of the target region,
                         a tuple of the interval of <start>,<length>, or a list of
                         tuples of <start>,<length>
        :return self
        """
        return self.update({"SEQUENCE_TARGET": self._parse_interval(interval)})

    def excluded(
        self, interval: Union[str, Tuple[int, int], List[Tuple[int, int]]]
    ) -> DesignPresets:
        """
        Primers and oligos may not overlap any region specified in this tag.
        The associated value must be a space-separated list of <start>,<length> pairs
        where <start> is the index of the first base of the
        excluded region, and <length> is its length. This tag is useful for tasks such
        as excluding regions of low sequence quality or for excluding regions containing
        repetitive elements such as ALUs or LINEs.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_TEMPLATE
        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_EXCLUDED_REGION

        :param interval: One of the following: the sequence of the target region,
                         a tuple of the interval of <start>,<length>, or a list of
                         tuples of <start>,<length>
        :return: self
        """
        return self.update({"SEQUENCE_EXCLUDED_REGION": self._parse_interval(interval)})

    def pick_anyway(self, b=1) -> DesignPresets:
        """
        If true use primer provided in SEQUENCE_PRIMER,
        SEQUENCE_PRIMER_REVCOMP, or SEQUENCE_INTERNAL_OLIGO even if it violates
        specific constraints.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_ANYWAY

        :param b: default True
        :return self
        """
        return self.update({"PRIMER_PICK_ANYWAY": b})


def clip(x, mn, mx):
    return max(min(x, mx), mn)


class DesignBase:

    DEFAULT_PARAMS = default_boulderio()  #: default parameters
    DEFAULT_GRADIENT = dict(
        PRIMER_MAX_SIZE=(1, DEFAULT_PARAMS["PRIMER_MAX_SIZE"], 36),
        PRIMER_MIN_SIZE=(-1, 16, DEFAULT_PARAMS["PRIMER_MAX_SIZE"]),
        PRIMER_MAX_TM=(1, DEFAULT_PARAMS["PRIMER_MAX_SIZE"], 80),
        PRIMER_MIN_TM=(-1, 48, DEFAULT_PARAMS["PRIMER_MIN_TM"]),
        PRIMER_MAX_HAIRPIN_TH=(1, DEFAULT_PARAMS["PRIMER_MAX_HAIRPIN_TH"], 60),
    )  #: the default gradient to use for the :meth:`Design.run_and_optimize` method.
    _CHECK_PRIMERS = "check_primers"
    _GENERIC = "generic"
    _PICK_PRIMER_LIST = "pick_primer_list"
    _PICK_SEQUENCING_PRIMERS = "pick_sequencing_primers"
    _PICK_CLONING_PRIMERS = "pick_cloning_primers"
    _PICK_DISCRIMINATIVE_PRIMERS = "pick_discriminative_primers"

    def __init__(self):
        self.params = self.DEFAULT_PARAMS.copy()
        self.logger = logger(self)

    def run(self, params=None) -> Tuple[List[Dict], List[Dict]]:
        """Design primers. Optionally provide additional parameters.

        :param params:
        :return: results
        """
        if params is None:
            params = self.params
        res = primer3.bindings.designPrimers(params._sequence(), params._globals())
        pairs, explain = parse_primer3_results(res)
        return pairs, explain

    def run_and_optimize(
        self, max_iterations, params=None, gradient=None
    ) -> Tuple[List[dict], List[dict]]:
        """Design primers. If primer design is unsuccessful, relax parameters
        as defined in primer3plust.Design.DEFAULT_GRADIENT. Repeat for the specified
        number of max_iterations.

        :param max_iterations: the max number of iterations to perform relaxation
        :param params: optional parameters to provide
        :param gradient: optional gradient to provide. If not provided,
                            Design.DEFAULT_GRADIENT will be used. The gradient is a
                            dictionary off 3 tuples, the step the min and the max.
        :return: results
        """
        if gradient is None:
            gradient = self.DEFAULT_GRADIENT
        if params is None:
            params = self.params
        # n_return = params["PRIMER_NUM_RETURN"]
        pairs, explain = self.run(params)
        i = 0
        while i < max_iterations and len(pairs) == 0:
            i += 1
            update = self._update_dict(params, gradient=gradient)
            if update:
                self.logger.info("Updated: {}".format(update))
            else:
                break
            self.params.update(update)
            pairs, explain = self.run(params)
        return pairs, explain

    @staticmethod
    def _update_dict(params, gradient):
        update = {}
        for param_key, gradient_tuple in gradient.items():
            delta, mn, mx = gradient_tuple
            try:
                val = params[param_key] + delta
                val = clip(val, mn, mx)
                if params[param_key] != val:
                    update[param_key] = val
            except Exception as e:
                raise e
        return update

    @staticmethod
    def open_help():
        """Open the documentation help in a new browser tab."""
        webbrowser.open(DOCURL)

    def copy(self):
        """Copy this design and its parameters."""
        designer = self.__class__()
        designer.params = self.params.copy()

    def __copy__(self):
        return self.copy()


class Design(DesignBase, AllParameters):
    def __init__(self):
        """Initialize a new design. Set parameters using
        :meth:`Design.set`, which
        returns an instance of :class:`DesignPresets <primer3plus.design.DesignPresets>`

        .. code-block::

            design = Design()
            design.presets.template("AGGCTGTAGTGCTTGTAGCTGGTTGCGTTACTGTG")
            design.presets.left_sequence("GTAGTGCTTGTA")
            design.run()
        """
        super().__init__()
        self._set = DesignPresets(self)

    def set(self, key, value):
        self.params.defs[key].value = value

    def get(self, key):
        return self.params.defs[key]

    @property
    def presets(self) -> DesignPresets:
        """Return the :class:`DesignPresets <primer3plus.design.DesignPresets>`
        instance for this design."""
        return self._set

    def update(self, data: Dict[str, Any]):
        return self.params.update(data)


def new(params=None):
    """Start a new design"""
    design = Design()
    if params:
        design.params.update(params)
