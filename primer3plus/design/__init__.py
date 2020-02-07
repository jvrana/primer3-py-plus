"""Core design modules.

For examples for how to use the design module, see the :doc:`Usage Docs <../usage>`
For a list of design parameters available, take a look at the
:ref:`BoulderIO Parameters <api_default_parameters>`

.. code-block:: python

    # a new task
    design = Design()

    # set template sequence
    design.settings.template("AGGTTGCGTGTGTATGGTCGTGTAGTGTGT")

    # set left primer sequence
    design.settings.left_sequence("GTTGCGTGTGT)

    # set as a cloning task
    design.settings.as_cloning_task()

    # run the design task
    design.run()
"""
import re
import webbrowser
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from warnings import warn

import primer3

from .interfaces import AllParameters
from .interfaces import ParameterAccessor
from .results import parse_primer3_results
from primer3plus.constants import DOCURL
from primer3plus.exceptions import Primer3PlusException
from primer3plus.exceptions import Primer3PlusRunTimeError
from primer3plus.exceptions import Primer3PlusWarning
from primer3plus.log import logger
from primer3plus.params import BoulderIO
from primer3plus.params import default_boulderio
from primer3plus.utils import anneal as anneal_primer
from primer3plus.utils import depreciated_warning


class DesignPresets:
    """Interface for setting design parameters. This is typically accessed from
    a :class:`Design <primer3plus.Design>` instance's.

    :meth:`Design <primer3plus.Design.set>` method. As in:

    .. code-block::

        design = Design()
        design.settings.left_sequence("AGGGAGATAGATA")
        design.run()
    """

    def __init__(self, design):
        """Initializes a new interface from a.

        :class:`~primer3plus.design.Design`.

        :param design: The design
        """
        self._design = design

    def _resolve(self):
        """Process any extra parameters and process BoulderIO so that it is
        digestable by primer3."""
        if self._design.PRIMER_USE_OVERHANGS.value:
            self._resolve_overhangs(self._design.PRIMER_MIN_ANNEAL_CHECK.value)
        if self._design.PRIMER_LONG_OK.value:
            self._resolve_max_lengths(lim=BoulderIO.PRIMER_MAX_SIZE_HARD_LIM)
            self._resolve_product_sizes()
        if not self._design.PRIMER_USE_OVERHANGS.value:
            if self._design.SEQUENCE_PRIMER_OVERHANG.value:
                warn(
                    Primer3PlusWarning(
                        "{} is non-empty (value={}) but {} was False. Overhang was"
                        " ignored".format(
                            self._design.SEQUENCE_PRIMER_OVERHANG.name,
                            self._design.SEQUENCE_PRIMER_OVERHANG.value,
                            self._design.PRIMER_USE_OVERHANGS.name,
                        )
                    )
                )
            if self._design.SEQUENCE_PRIMER_REVCOMP_OVERHANG.value:
                warn(
                    Primer3PlusWarning(
                        "{} is non-empty (value={}) but {} was False. Overhang was"
                        " ignored".format(
                            self._design.SEQUENCE_PRIMER_REVCOMP_OVERHANG.name,
                            self._design.SEQUENCE_PRIMER_REVCOMP_OVERHANG.value,
                            self._design.PRIMER_USE_OVERHANGS.name,
                        )
                    )
                )
        return self

    def _post_parse(self, pairs, explain) -> None:
        """Modify results from design parameters (e.g. overhangs)"""
        left_long_overhang = self._design._SEQUENCE_LONG_OVERHANG.value
        right_long_overhang = self._design._SEQUENCE_REVCOMP_LONG_OVERHANG.value

        for pair in pairs.values():
            for x in ["LEFT", "RIGHT"]:
                pair[x].setdefault("OVERHANG", "")
            if self._design.PRIMER_USE_OVERHANGS:
                pair["LEFT"]["OVERHANG"] = self._design.SEQUENCE_PRIMER_OVERHANG.value
                pair["RIGHT"][
                    "OVERHANG"
                ] = self._design.SEQUENCE_PRIMER_REVCOMP_OVERHANG.value
                if left_long_overhang:
                    pair["LEFT"]["SEQUENCE"] = (
                        left_long_overhang + pair["LEFT"]["SEQUENCE"]
                    )
                    pair["LEFT"]["OVERHANG"] = pair["LEFT"]["OVERHANG"][
                        : -len(left_long_overhang)
                    ]
                    pair["PAIR"]["PRODUCT_SIZE"] += len(left_long_overhang)
                    loc = pair["LEFT"]["location"]
                    pair["LEFT"]["location"] = [
                        loc[0] - len(left_long_overhang),
                        len(pair["LEFT"]["SEQUENCE"]),
                    ]

                if right_long_overhang:
                    pair["RIGHT"]["SEQUENCE"] = (
                        right_long_overhang + pair["RIGHT"]["SEQUENCE"]
                    )
                    pair["RIGHT"]["OVERHANG"] = pair["RIGHT"]["OVERHANG"][
                        : -len(right_long_overhang)
                    ]
                    pair["PAIR"]["PRODUCT_SIZE"] += len(right_long_overhang)
                    loc = pair["RIGHT"]["location"]
                    pair["RIGHT"]["location"] = [
                        loc[0] + len(right_long_overhang),
                        len(pair["RIGHT"]["SEQUENCE"]),
                    ]

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
        """Update an arbitrary parameter."""
        self._design.params.update(update)
        return self

    def task(self, task: str) -> "DesignPresets":
        """This tag tells primer3 what task to perform.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_TASK

        :param task: the task name
        :return self
        """
        self.update({"PRIMER_TASK": task})
        return self

    def as_cloning_task(self) -> "DesignPresets":
        """Set the design as a cloning task.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_TASK

        :return: self
        """
        return self.task("pick_cloning_primers")

    def as_generic_task(self) -> "DesignPresets":
        """Set the design as a generic task.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_TASK

        :return: self
        """
        return self.task("generic")

    def template(self, template: str) -> "DesignPresets":
        """Set the template sequence for the design. This sets the
        'SEQUENCE_TEMPLATE' parameter.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_TEMPLATE

        :param template: the template sequence
        :return: self
        """
        self.update({"SEQUENCE_TEMPLATE": template})
        return self

    # TODO: set_iterations, set_num_return, set_force_return, set_gradient
    def primer_num_return(self, n: int) -> "DesignPresets":
        """Set the number of primers to return for the design task.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_NUM_RETURN

        :param n: number of primers to return
        :return: self
        """
        return self.update({"PRIMER_NUM_RETURN": n})

    def product_size(
        self, interval: Union[Tuple[int, int], List[Tuple[int, int]]], opt=None
    ) -> "DesignPresets":
        """Set the product size. Optionally include the optimal size.

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
    ) -> "DesignPresets":
        """The list of regions from which to design primers.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER_PAIR_OK_REGION_LIST

        :param region_list: list of regions
        :return: self
        """
        return self.update({"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": region_list})

    def left_sequence(self, primer: str) -> "DesignPresets":
        """The sequence of a left primer to check and around which to design
        right primers and optional internal oligos. Must be a substring of
        SEQUENCE_TEMPLATE.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER

        :param primer: the primer sequence
        :return: self
        """

        return self.update({"SEQUENCE_PRIMER": primer, "PRIMER_PICK_RIGHT_PRIMER": 1})

    def right_sequence(self, primer: str) -> "DesignPresets":
        """The sequence of a right primer to check and around which to design
        left primers and optional internal oligos. Must be a substring of the
        reverse strand of SEQUENCE_TEMPLATE.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_PRIMER_REVCOMP
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER

        :param primer: the primer sequence
        :return: self
        """
        return self.update(
            {"SEQUENCE_PRIMER_REVCOMP": primer, "PRIMER_PICK_LEFT_PRIMER": 1}
        )

    @staticmethod
    def _trim_long(overhang: str, anneal: str, lim: int) -> Tuple[str, str, str]:
        """Fix the overhang and anneal from the hardcoded BoulderIO primer
        lim."""
        return overhang, anneal[:-lim], anneal[-lim:]

    def _get_left_overhang(self, min_primer_anneal: int):
        left = self._design.SEQUENCE_PRIMER.value
        if left:
            fwd, _ = anneal_primer(
                self._design.SEQUENCE_TEMPLATE.value, [left], n_bases=min_primer_anneal
            )
            if len(fwd) == 0:
                raise Primer3PlusRunTimeError("No annealing found for left sequence.")
            elif len(fwd) > 1:
                raise Primer3PlusRunTimeError(
                    "More than one annealing found for left sequence."
                )
            overhang = fwd[0]["overhang"]
            anneal = fwd[0]["anneal"]
            return overhang, anneal
        else:
            return "", left

    def _get_right_overhang(self, min_primer_anneal: int):
        right = self._design.SEQUENCE_PRIMER_REVCOMP.value
        if right:
            _, rev = anneal_primer(
                self._design.SEQUENCE_TEMPLATE.value, [right], n_bases=min_primer_anneal
            )
            if len(rev) == 0:
                raise Primer3PlusRunTimeError("No annealing found for right sequence.")
            elif len(rev) > 1:
                raise Primer3PlusRunTimeError(
                    "More than one annealing found for right "
                    "sequence {}.".format(self._design.SEQUENCE_PRIMER_REVCOMP)
                )
            overhang = rev[0]["overhang"]
            anneal = rev[0]["anneal"]
            return overhang, anneal
        else:
            return "", right

    def left_overhang(self, overhang: str) -> "DesignPresets":
        """Sets the left overhang sequence for the primer. This overhang will.

        *always* be in the overhang sequence regardless of other parameters.

        If using a primer that anneals with an overhang, this value will
        be appended to the 5' end of the overhang.

        :param overhang: overhang sequence
        :return: self
        """
        return self.update({"SEQUENCE_PRIMER_OVERHANG": overhang})

    def right_overhang(self, overhang: str) -> "DesignPresets":
        """Sets the right overhang sequence for the primer. This overhang will.

        *always* be in the overhang sequence regardless of other parameters.

        If using a primer that anneals with an overhang, this value will
        be appended to the 5' end of the overhang.

        :param overhang: overhang sequence
        :return: self
        """
        return self.update({"SEQUENCE_PRIMER_REVCOMP_OVERHANG": overhang})

    def use_overhangs(self, b: bool = True) -> "DesignPresets":
        """Set the BoulderIO to process overhangs.

        :param b: boolean to set
           :return: self
        """
        return self.update({"PRIMER_USE_OVERHANGS": b})

    def long_ok(self, b: bool = True) -> "DesignPresets":
        """Set the BoulderIO to process long primers.

        :param b: boolean to set
        :return: self
        """
        return self.update({"PRIMER_LONG_OK": b})

    def _resolve_product_sizes(self):
        """If there are long primers being used, the product_size is no longer
        valid as the trimmed sequence is no longer represented in the
        originally provided product size.

        This re-adjusts the product size to correspond the adjusted
        parameters.
        """
        # adjust product size range
        left_long_overhang = self._design._SEQUENCE_LONG_OVERHANG.value
        right_long_overhang = self._design._SEQUENCE_REVCOMP_LONG_OVERHANG.value
        product_sizes = self._design.PRIMER_PRODUCT_SIZE_RANGE.value

        x = len(left_long_overhang) + len(right_long_overhang)

        if isinstance(product_sizes[0], tuple):
            new_product_sizes = []
            for size in product_sizes:
                new_product_sizes.append((size[0] - x, size[1] - x))
            self._design.PRIMER_PRODUCT_SIZE_RANGE.value = new_product_sizes
        else:
            size = self._design.PRIMER_PRODUCT_SIZE_RANGE.value
            self._design.PRIMER_PRODUCT_SIZE_RANGE.value = [size[0] - x, size[1] - x]

    def _resolve_max_lengths(self, lim: int):
        """Fixes the annealing and overhang sequences for annealing sequences
        for primers over the :attr:`BoulderIO.

        <primer3plus.paramsBoulderIO.PRIMER_MAX_SIZE_HARD_LIM>`.
        Should always be run *after* :meth:`_resolve_overhangs`.
        """
        left_anneal = self._design.SEQUENCE_PRIMER.value
        right_anneal = self._design.SEQUENCE_PRIMER_REVCOMP.value
        left_over = self._design.SEQUENCE_PRIMER_OVERHANG.value
        right_over = self._design.SEQUENCE_PRIMER_REVCOMP_OVERHANG.value

        left_over, left_long_overhang, left_anneal = self._trim_long(
            left_over, left_anneal, lim=lim
        )
        right_over, right_long_overhang, right_anneal = self._trim_long(
            right_over, right_anneal, lim=lim
        )

        self.left_overhang(left_over + left_long_overhang)
        self.right_overhang(right_over + right_long_overhang)

        # save the sequences that were trimmed
        # TODO: will need to re-add these in the results in overhang, product size, and anneal
        # TODO: adjust the product_size
        # TODO: adjust any other regions
        # TODO: re-adjust tm and add any warnings

        # save values for long overhangs
        self._left_long_overhang(left_long_overhang)
        self._right_long_overhang(right_long_overhang)

        self.left_sequence(left_anneal)
        self.right_sequence(right_anneal)

    def _left_long_overhang(self, x):
        self.update({"_SEQUENCE_LONG_OVERHANG": x})

    def _right_long_overhang(self, x):
        self.update({"_SEQUENCE_REVCOMP_LONG_OVERHANG": x})

    def _resolve_overhangs(self, min_primer_anneal: int):
        """Sets the annealing and overhang sequences."""
        left_over, left_anneal = self._get_left_overhang(min_primer_anneal)
        _loverhang = self._design.SEQUENCE_PRIMER_OVERHANG.value
        if _loverhang:
            left_over = _loverhang + left_over
            # raise ValueError(
            #     "Left overhang already set to '{}'.".format(_loverhang)
            # )

        right_over, right_anneal = self._get_right_overhang(min_primer_anneal)
        _roverhang = self._design.SEQUENCE_PRIMER_REVCOMP_OVERHANG.value
        if _roverhang:
            right_over = _roverhang + right_over
            # raise ValueError(
            #     "Right overhang already set to '{}'.".format(_roverhang)
            # )

        self.left_overhang(left_over)
        self.right_overhang(right_over)
        self.left_sequence(left_anneal)
        self.right_sequence(right_anneal)

    def pick_left_only(self) -> "DesignPresets":
        """Design only the left primer.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER

        :return: self
        """
        return self.update(
            {"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0}
        )

    def pick_right_only(self) -> "DesignPresets":
        """Design only the right primer.

        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_LEFT_PRIMER
        http://primer3.ut.ee/primer3web_help.htm#PRIMER_PICK_RIGHT_PRIMER

        :return: self
        """
        return self.update(
            {"PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1}
        )

    def internal_sequence(self, primer: str) -> "DesignPresets":
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

    def primers(self, p1: str, p2: str) -> "DesignPresets":
        """Set the left and right primer sequences.

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

    def primers_with_overhangs(self, p1: str, p2: str) -> "DesignPresets":
        if p1:
            self.left_sequence_with_overhang(p1)
        if p2:
            self.right_sequence_with_overhang(p2)
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

    def included(self, interval: Union[str, Tuple[int, int]]) -> "DesignPresets":
        """Specify interval from which primers must be selected. A sub-region
        of the given sequence in which to pick primers. For example, often the
        first dozen or so bases of a sequence are vector, and should be
        excluded from consideration. The value for this parameter has the form.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_INCLUDED_REGION

        <start>,<length> where <start> is the index of the first base to consider, and
        <length> is the number of subsequent bases in the primer-picking region.

        :param interval: One of the following: the sequence of the target region,
                         a tuple of the interval of <start>,<length> or a str
        :return: self
        """
        if isinstance(interval, str):
            interval = self._interval_from_sequences(
                self._design.params["SEQUENCE_TEMPLATE"], interval
            )
        if not len(interval) == 2 or (
            not isinstance(interval, tuple) and not isinstance(interval, list)
        ):
            raise TypeError(
                "Expect an tuple or list of length 2 but found {}".format(interval)
            )
        interval = list(interval)
        return self.update({"SEQUENCE_INCLUDED_REGION": interval})

    def target(
        self, interval: Union[str, Tuple[int, int], List[Tuple[int, int]]]
    ) -> "DesignPresets":
        """Specify the interval that designed primers must flank. If one or
        more targets is specified then a legal primer pair must flank at least
        one of them. A target might be a simple sequence repeat site (for
        example a CA repeat) or a single-base-pair polymorphism, or an exon for
        resequencing. The value should be a space-separated list of.

        <start>,<length> pairs where <start> is the index of the first base of
        a target,and <length> is its length.  See also PRIMER_INSIDE_PENALTY,
        PRIMER_OUTSIDE_PENALTY. PRIMER_TASK=pick_sequencing_primers. See
        PRIMER_TASK for more information.

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
    ) -> "DesignPresets":
        """Primers and oligos may not overlap any region specified in this tag.
        The associated value must be a space-separated list of <start>,<length>
        pairs where <start> is the index of the first base of the excluded
        region, and <length> is its length. This tag is useful for tasks such
        as excluding regions of low sequence quality or for excluding regions
        containing repetitive elements such as ALUs or LINEs.

        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_TEMPLATE
        http://primer3.ut.ee/primer3web_help.htm#SEQUENCE_EXCLUDED_REGION

        :param interval: One of the following: the sequence of the target region,
                         a tuple of the interval of <start>,<length>, or a list of
                         tuples of <start>,<length>
        :return: self
        """
        return self.update({"SEQUENCE_EXCLUDED_REGION": self._parse_interval(interval)})

    def pick_anyway(self, b=1) -> "DesignPresets":
        """If true use primer provided in SEQUENCE_PRIMER,
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
    """Base design."""

    DEFAULT_PARAMS = default_boulderio  #: default parameters
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

    def __init__(
        self,
        gradient: Dict[
            str, Tuple[Union[float, int], Union[float, int], Union[float, int]]
        ] = None,
        params: BoulderIO = None,
        quiet_runtime: bool = False,
    ):
        """Initializes a new design.

        :param gradient: the design gradient.
        :param quiet_runtime: if True will siliently ignore any runtime errors.
        """
        if params is None:
            params = self.DEFAULT_PARAMS.copy()
        self.params = params
        self.logger = logger(self)
        self.gradient = gradient
        self.quiet_runtime = quiet_runtime

    def _raise_run_time_error(self, msg: str) -> Primer3PlusRunTimeError:
        """Raise a Primer3PlusRunTime exception. If parameters are named in the
        msg, print off some debugging information at the end of the message.

        :param msg: the error msg
        :return: the run time exception
        """
        parameter_explain = set()
        for name, value in self.params._params.items():
            if name in msg:
                parameter_explain.add("\t" + str(value))
        parameter_explain = sorted(parameter_explain)
        return Primer3PlusRunTimeError(msg + "\n" + "\n".join(parameter_explain))

    def _run(self, params: BoulderIO = None) -> Tuple[List[Dict], List[Dict]]:
        """Design primers. Optionally provide additional parameters.

        :param params:
        :return: results
        """
        if params is None:
            params = self.params
        try:
            res = primer3.bindings.designPrimers(params._sequence(), params._globals())
        except OSError as e:
            if not self.quiet_runtime:
                raise self._raise_run_time_error(str(e)) from e
            else:
                return {}, {"PRIMER_ERROR": str(e)}
        except Primer3PlusRunTimeError as e:
            if not self.quiet_runtime:
                raise self._raise_run_time_error(str(e)) from e
            else:
                return {}, {"PRIMER_ERROR": str(e)}
        except Primer3PlusException as e:
            raise self._raise_run_time_error(str(e)) from e

        pairs, explain = parse_primer3_results(res)
        self.settings._post_parse(pairs, explain)
        return pairs, explain

    def run(self) -> Tuple[List[Dict], List[Dict]]:
        """Design primers. Optionally provide additional parameters.

        :param params:
        :return: results
        """
        return self._run()

    def run_and_optimize(
        self,
        max_iterations,
        params: BoulderIO = None,
        gradient: Dict[
            str, Tuple[Union[float, int], Union[float, int], Union[float, int]]
        ] = None,
        run_kwargs: dict = None,
    ) -> Tuple[List[dict], List[dict]]:
        """Design primers and relax constraints. If primer design is
        unsuccessful, relax parameters as defined in
        primer3plust.Design.DEFAULT_GRADIENT. Repeat for the specified number
        of max_iterations.

        :param max_iterations: the max number of iterations to perform relaxation
        :param params: optional parameters to provide
        :param gradient: optional gradient to provide. If not provided,
                            Design.DEFAULT_GRADIENT will be used. The gradient is a
                            dictionary off 3 tuples, the step the min and the max.
        :return: results
        """
        if gradient is None:
            gradient = self.gradient or self.DEFAULT_GRADIENT
        if params is None:
            params = self.params
        pairs, explain = self._run(params)
        i = 0
        while i < max_iterations and len(pairs) == 0:
            i += 1
            update = self._update_dict(params, gradient=gradient)
            if update:
                self.logger.info("Updated: {}".format(update))
            else:
                self.logger.info("Reached end of gradient.")
                break
            self.params.update(update)
            pairs, explain = self._run(params)
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


class RestoreAfterRun:
    """Class to restore boulderio to its original parameters after a run."""

    def __init__(self, boulderio):
        self.params = boulderio

    def __enter__(self):
        for v in self.params._params.values():
            v.hold_restore()

    def __exit__(self, a, b, c):
        for v in self.params._params.values():
            v.restore()


class Design(DesignBase, AllParameters):
    def __init__(self):
        """Initialize a new design. Set parameters using.

        :attr:`Design.settings`, which
        returns an instance of
        :class:`DesignPresets <primer3plus.design.DesignPresets>`.

        Alternatively, parameters can be accessed more directly using
        the name of the parameter descriptor. For a list of parameters available, see
        :ref:`BoulderIO Parameters <api_default_parameters>`.

        .. code-block::

            design = Design()
            design.settings.template("AGGCTGTAGTGCTTGTAGCTGGTTGCGTTACTGTG")
            design.settings.left_sequence("GTAGTGCTTGTA")
            design.SEQUENCE_ID.value = "MY ID"
            design.run()
        """
        super().__init__()
        self._settings = DesignPresets(self)

    def set(self, key, value):
        self.params.defs[key].value = value

    def get(self, key):
        return self.params.defs[key]

    @property
    def settings(self) -> "DesignPresets":
        """Return the :class:`DesignPresets <primer3plus.design.DesignPresets>`
        instance for this design."""
        return self._settings

    @property
    def presets(self):
        depreciated_warning("'presets' has been renamed to 'settings'")
        return self.settings

    def update(self, data: Dict[str, Any]):
        """Update an arbitrary parameter."""
        return self.params.update(data)

    def run(self) -> Tuple[List[Dict], List[Dict]]:
        """Design primers. Optionally provide additional parameters.

        :param params:
        :return: results
        """
        with RestoreAfterRun(self.params):
            self.settings._resolve()
            return super()._run(None)

    def run_and_optimize(
        self,
        max_iterations,
        params: BoulderIO = None,
        gradient: Dict[
            str, Tuple[Union[float, int], Union[float, int], Union[float, int]]
        ] = None,
    ) -> Tuple[List[dict], List[dict]]:
        """Design primers. If primer design is unsuccessful, relax parameters
        as defined in primer3plust.Design.DEFAULT_GRADIENT. Repeat for the
        specified number of max_iterations.

        :param max_iterations: the max number of iterations to perform relaxation
        :param params: optional parameters to provide
        :param gradient: optional gradient to provide. If not provided,
                            Design.DEFAULT_GRADIENT will be used. The gradient is a
                            dictionary off 3 tuples, the step the min and the max.
        :return: results
        """
        with RestoreAfterRun(self.params):
            self.settings._resolve()
            return super().run_and_optimize(max_iterations)


def new(params=None):
    """Start a new design."""
    design = Design()
    if params:
        design.params.update(params)
