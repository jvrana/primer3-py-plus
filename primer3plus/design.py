"""
To run designs, you copy the underlying parameters and change certain values.

There should be a number of presets available to run certain design tasks

There should also be a relaxation parameter (a subclass of Design?)

Async option

Design with overhangs
"""

from functools import wraps
import itertools
import webbrowser
from collections import Counter
import primer3
from .results import parse_primer3_results
from .params import default_boulderio
from typing import Tuple, List, Dict
import re


def _summarize_reasons(reasons):
    reason_dict = {}
    for reason in reasons:
        for k, v in reason.items():
            if "EXPLAIN" in k:
                for m in re.finditer("\s*([\w\s\-]+)\s+(\d+)", v):
                    reason_token = m.group(1)
                    num = int(m.group(2))
                    reason_dict.setdefault(k, Counter())[reason_token] += num
    return {k: dict(v) for k, v in reason_dict.items()}


def combine_results(results):
    """
    Combine and sort results. Combine all explainations.

    :param results:
    :type results:
    :return:
    :rtype:
    """
    all_pairs = []
    all_reasons = []
    for r in results:
        all_pairs += list(r[0].values())
        all_reasons.append(r[1])
    sorted_pairs = sorted(all_pairs, key=lambda x: x["PAIR"]["PENALTY"])
    explain = _summarize_reasons(all_reasons)
    return sorted_pairs, explain


def dispatch_iterable(params, max_results=10):
    """
    Dispatches the function to combine and sorts results if an argument is an iterable indicated in `params`.

    e.g. `params = [(1, (type, list)]` will dispatch the function if the first argument is a list or tuple.
    """

    def wrapped(f):
        @wraps(f)
        def _wrapped(*args, **kwargs):
            new_args = []
            is_iterable = False
            for pi, arg in enumerate(args):
                s = False
                for i, types in params:
                    if pi == i and type(arg) in types:
                        new_args.append(arg)
                        s = True
                        is_iterable = True
                if not s:
                    new_args.append([arg])
            if is_iterable:
                iterable_args = itertools.product(*new_args)
                results = []
                for _args in iterable_args:
                    if max_results > 0 and len(results) > max_results:
                        break
                    results.append(f(*_args, **kwargs))
                return combine_results(results)
            else:
                return f(*args, **kwargs)

        return _wrapped

    return wrapped


class DesignPresets(object):
    def __init__(self, design):
        self._design = design

    def _sequence_from_template(self, template, target):
        if isinstance(target, str):
            matches = self._get_index_of_match(template, target)
            if not matches:
                print("Target not in template")
                return None
            if len(matches) > 1:
                print("More than one target found")
                return None
            return matches[0]

    def _get_index_of_match(self, template, sequence):
        matches = []
        for m in re.finditer(sequence, template, re.IGNORECASE):
            matches.append((m.start(0), m.end(0)))
        return matches

    def _set(self, update: dict):
        self._design.params.update(update)
        return self

    def task(self, task: str):
        """
        This tag tells primer3 what task to perform. Legal values are:

        generic

        Design a primer pair (the classic primer3 task) if the PRIMER_PICK_LEFT_PRIMER=1, and
        PRIMER_PICK_RIGHT_PRIMER=1. In addition, pick an internal hybridization oligo if PRIMER_PICK_INTERNAL_OLIGO=1.

        NOTE: If PRIMER_PICK_LEFT_PRIMER=1, PRIMER_PICK_RIGHT_PRIMER=0 and PRIMER_PICK_INTERNAL_OLIGO=1, then behaves
        similarly to PRIMER_TASK=pick_primer_list.

        pick_detection_primers

        Deprecated alias for PRIMER_TASK=generic

        check_primers

        Primer3 only checks the primers provided in SEQUENCE_PRIMER, SEQUENCE_INTERNAL_OLIGO and
        SEQUENCE_PRIMER_REVCOMP. It is the only task that does not require a sequence. However, if SEQUENCE_TEMPLATE is
        provided, it is used.

        pick_primer_list

        Pick all primers in SEQUENCE_TEMPLATE (possibly limited by SEQUENCE_INCLUDED_REGION,
        SEQUENCE_EXCLUDED_REGION, SEQUENCE_PRIMER_PAIR_OK_REGION_LIST, etc.). Returns the primers sorted by quality
        starting with the best primers. If PRIMER_PICK_LEFT_PRIMER and PRIMER_PICK_RIGHT_PRIMER is selected primer3
        does not to pick primer pairs but generates independent lists of left primers, right primers, and, if requested,
         internal oligos.

        pick_sequencing_primers

        Pick primers suited to sequence a region. SEQUENCE_TARGET can be used to indicate several targets. The
        position of each primer is calculated for optimal sequencing results.

        pick_cloning_primers

        Pick primers suited to clone a gene were the start nucleotide and the end nucleotide of the PCR fragment must be
         fixed, for example to clone an ORF. SEQUENCE_INCLUDED_REGION must be used to indicate the first and the last
         nucleotide. Due to these limitations primer3 can only vary the length of the primers. Set PRIMER_PICK_ANYWAY=1
          to obtain primers even if they violate specific constraints.

        pick_discriminative_primers

        Pick primers suited to select primers which bind with their end at a specific position. This can be used to
        force the end of a primer to a polymorphic site, with the goal of discriminating between sequence variants.
        SEQUENCE_INCLUDED_REGION must be used to indicate the last nucleotide of the left (first nucleotide of included
        region) and the right primer (last nucleotide of included region). Due to these limitations primer3 can only
        vary the length of the primers. Set PRIMER_PICK_ANYWAY=1 to obtain primers even if they violate specific
        constraints.

        pick_pcr_primers

        Deprecated shortcut for the following settings:
        PRIMER_TASK=generic
        PRIMER_PICK_LEFT_PRIMER=1
        PRIMER_PICK_INTERNAL_OLIGO=0
        PRIMER_PICK_RIGHT_PRIMER=1

        WARNING: this task changes the values of PRIMER_PICK_LEFT_PRIMER, PRIMER_PICK_INTERNAL_OLIGO, and
        PRIMER_PICK_RIGHT_PRIMER in a way that is not obvious by looking at the input.

        pick_pcr_primers_and_hyb_probe

        Deprecated shortcut for the following settings:
        PRIMER_TASK=generic
        PRIMER_PICK_LEFT_PRIMER=1
        PRIMER_PICK_INTERNAL_OLIGO=1
        PRIMER_PICK_RIGHT_PRIMER=1

        WARNING: this task changes the values of PRIMER_PICK_LEFT_PRIMER, PRIMER_PICK_INTERNAL_OLIGO, and
        PRIMER_PICK_RIGHT_PRIMER in a way that is not obvious by looking at the input.

        pick_left_only

        Deprecated shortcut for the following settings:
        PRIMER_TASK=generic
        PRIMER_PICK_LEFT_PRIMER=1
        PRIMER_PICK_INTERNAL_OLIGO=0
        PRIMER_PICK_RIGHT_PRIMER=0

        WARNING: this task changes the values of PRIMER_PICK_LEFT_PRIMER, PRIMER_PICK_INTERNAL_OLIGO, and
        PRIMER_PICK_RIGHT_PRIMER in a way that is not obvious by looking at the input.

        pick_right_only

        Deprecated shortcut for the following settings:
        PRIMER_TASK=generic
        PRIMER_PICK_LEFT_PRIMER=0
        PRIMER_PICK_INTERNAL_OLIGO=0
        PRIMER_PICK_RIGHT_PRIMER=1

        WARNING: this task changes the values of PRIMER_PICK_LEFT_PRIMER, PRIMER_PICK_INTERNAL_OLIGO, and
        PRIMER_PICK_RIGHT_PRIMER in a way that is not obvious by looking at the input.

        pick_hyb_probe_only

        Deprecated shortcut for the following settings:
        PRIMER_TASK=generic
        PRIMER_PICK_LEFT_PRIMER=0
        PRIMER_PICK_INTERNAL_OLIGO=1
        PRIMER_PICK_RIGHT_PRIMER=0

        WARNING: this task changes the values of PRIMER_PICK_LEFT_PRIMER, PRIMER_PICK_INTERNAL_OLIGO, and
        PRIMER_PICK_RIGHT_PRIMER in a way that is not obvious by looking at the input.

        :param task:
        :type task:
        :return:
        :rtype:
        """
        self._set({"PRIMER_TASK": task})
        return self

    def as_cloning_task(self):
        return self.task("pick_cloning_primers")

    def as_generic_task(self):
        return self.task("generic")

    def template(self, template):
        self._set({"SEQUENCE_TEMPLATE": template})
        return self

    # TODO: set_iterations, set_num_return, set_force_return, set_gradient
    def primer_num_return(self, n):
        return self._set({"PRIMER_NUM_RETURN": n})

    def size(self, interval: tuple, opt=None):
        """
        Set the product size. Optionally include the optimal size.

        :param interval: a tuple of
        :type interval:
        :param opt:
        :type opt:
        :return:
        :rtype:
        """
        if opt is None:
            opt = int(sum(interval) / 2.0)
        return self._set(
            {"PRIMER_PRODUCT_SIZE_RANGE": interval, "PRIMER_PRODUCT_OPT_SIZE": opt}
        )

    def pair_region_list(self, region_list: List[Tuple[int, int, int, int]]):
        return self._set({"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": region_list})

    def left_sequence(self, primer: str):
        """
        The sequence of a left primer to check and around which to design right primers and optional internal
        oligos. Must be a substring of SEQUENCE_TEMPLATE.

        :param primer:
        :type primer:
        :return:
        :rtype:
        """
        return self._set({"SEQUENCE_PRIMER": primer, "PRIMER_PICK_RIGHT_PRIMER": 1})

    def right_sequence(self, primer: str):
        """
        The sequence of a right primer to check and around which to design left primers and optional internal
        oligos. Must be a substring of the reverse strand of SEQUENCE_TEMPLATE.

        :param primer:
        :type primer:
        :return:
        :rtype:
        """
        return self._set(
            {"SEQUENCE_PRIMER_REVCOMP": primer, "PRIMER_PICK_LEFT_PRIMER": 1}
        )

    def pick_left_only(self):
        return self._set({"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0})

    def pick_right_only(self):
        return self._set({"PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1})

    def internal_sequence(self, primer: str):
        """
        The sequence of an internal oligo to check and around which to design left and right primers. Must be a
        substring of SEQUENCE_TEMPLATE.

        :param primer:
        :type primer:
        :return:
        :rtype:
        """
        self._set({"SEQUENCE_INTERNAL_OLIGO": primer, "PRIMER_PICK_INTERNAL_OLIGO": 1})

    def primers(self, p1: str, p2: str):
        if p1:
            self.left_sequence(p1)
        if p2:
            self.right_sequence(p2)
        return self

    def included(self, interval: tuple):
        """
        A sub-region of the given sequence in which to pick primers. For example, often the first dozen or so bases
        of a sequence are vector, and should be excluded from consideration. The value for this parameter has the form

        <start>,<length>
        where <start> is the index of the first base to consider, and <length> is the number of subsequent bases
        in the primer-picking region.

        :param interval:
        :type interval:
        :return:
        :rtype:
        """
        return self._set({"SEQUENCE_INCLUDED_REGION": interval})

    def target(self, interval: tuple):
        """
        If one or more targets is specified then a legal primer pair must flank at least one of them. A target might
        be a simple sequence repeat site (for example a CA repeat) or a single-base-pair polymorphism, or an exon for
        resequencing. The value should be a space-separated list of

        <start>,<length>
        pairs where <start> is the index of the first base of a target, and <length> is its length.

        See also PRIMER_INSIDE_PENALTY, PRIMER_OUTSIDE_PENALTY. Has a different meaning when
        PRIMER_TASK=pick_sequencing_primers. See PRIMER_TASK for more information.

        :param interval:
        :type interval:
        :return:
        :rtype:
        """
        if isinstance(interval, str):
            interval = self._sequence_from_template(
                self._design.params["SEQUENCE_TEMPLATE"], interval
            )
        return self._set({"SEQUENCE_TARGET": interval})

    def excluded(self, interval):
        """
        Primers and oligos may not overlap any region specified in this tag. The associated value must be a
        space-separated list of

        <start>,<length>
        pairs where <start> is the index of the first base of the excluded region, and <length> is its length. This
        tag is useful for tasks such as excluding regions of low sequence quality or for excluding regions containing
        repetitive elements such as ALUs or LINEs.

        :param interval:
        :type interval:
        :return:
        :rtype:
        """
        if isinstance(interval, str):
            interval = self._sequence_from_template(
                self._design.params["SEQUENCE_TEMPLATE"], interval
            )
        return self._set({"SEQUENCE_EXCLUDED_REGION": interval})

    def pick_anyway(self, b=1):
        """
        If true use primer provided in SEQUENCE_PRIMER, SEQUENCE_PRIMER_REVCOMP, or SEQUENCE_INTERNAL_OLIGO even if
        it violates specific constraints.

        :param b:
        :type b:
        :return:
        :rtype:
        """
        return self._set({"PRIMER_PICK_ANYWAY": b})


class DesignBase(object):

    DEFAULT_PARAMS = default_boulderio()

    def __init__(self):
        self.params = self.DEFAULT_PARAMS.copy()

    def run(self, params=None) -> Tuple[List[Dict], List[Dict]]:
        if params is None:
            params = self.params
        results = primer3.bindings.designPrimers(params.sequence(), params.globals())
        pairs, explain = parse_primer3_results(results)
        return pairs, explain

    @staticmethod
    def open_help():
        webhelp = (
            "https://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/"
            + "primer3/src/libprimer3/primer3_manual.htm"
        )
        webbrowser.open(webhelp)

    def copy(self):
        designer = self.__class__()
        designer.params = self.params.copy()

    def __copy__(self):
        return self.copy()


class Design(DesignBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.set = DesignPresets(self)
