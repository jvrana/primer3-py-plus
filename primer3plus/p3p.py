import os
import re
import webbrowser
from collections import OrderedDict
from collections import namedtuple
from copy import deepcopy
from collections import Counter
import pandas as pd
import primer3
from typing import Dict, Any
import itertools
from functools import wraps
from primer3plus.log import logger

here = os.path.dirname(os.path.abspath(__file__))

# TODO: set relaxation ON or OFF
# TODO: group 'tasks' into an attribute
# TODO: better help documentation
# TODO: a CLI for quick design
# TODO: params should have attributes
# TODO: hardcode the attributes, with a helper method for determining primer version and boulderIO
# TODO: parse to boulderIO
class Primer3Params(object):
    """
    Reads the Primer3 documentation and creates the appropriate parameters.
    """

    POST_LOAD_DEFAULTS = {"PRIMER_EXPLAIN_FLAG": 1}
    default_file_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "primer3_params_raw.txt"
    )

    def __init__(self, params=None):
        if params is None:
            params = {}
        self.defaults = self._open_primer3_params()
        self.SEQUENCE = self._named_tuple_from_dict(
            "sequence", self.defaults["sequence"]
        )
        self.GLOBAL = self._named_tuple_from_dict("global", self.defaults["global"])
        self.PROGRAM = self._named_tuple_from_dict("program", self.defaults["program"])
        self.OTHER = self._named_tuple_from_dict("other", self.defaults["other"])
        self._all = ["SEQUENCE", "GLOBAL", "PROGRAM", "OTHER"]
        self.update(self.POST_LOAD_DEFAULTS)
        self.update(params)

    def reset(self):
        """
        Reset all parameters.

        :return: None
        """
        for key in self._all:
            setattr(self, key.upper(), self._named_tuple_from_dict(key.lower()))

    def update(self, value_dict: Dict[str, Any]):
        """
        Update parameters.

        :param value_dict:
        :return:
        """
        for key in self._all:
            param_tuple = getattr(self, key)
            param_dict = param_tuple._asdict()
            for k, v in value_dict.items():
                if k in param_dict:
                    param_dict[k]["value"] = v
            data = self.defaults[key.lower()]
            data.update(param_dict)
            new_tuple = self._named_tuple_from_dict(key, data)
            setattr(self, key.upper(), new_tuple)

    def __setitem__(self, key, value):
        self.update({key: value})

    def __getitem__(self, key):
        for param_type, param_dict in self.asdict().items():
            if key in param_dict:
                return param_dict[key]["value"]

    def __contains__(self, key):
        for param_type, param_dict in self.asdict().items():
            if key in param_dict:
                return True
        return False

    def _gettype(self, key):
        return getattr(self, key)

    def value(self, key):
        return self[key]["value"]

    def description(self, key):
        return self[key]["description"]

    def default(self, key):
        return self[key]["default"]

    def valuedict(self, clean=False):
        data = self.asdict(values_only=True)
        if clean:
            return {k: self._clean_params(v) for k, v in data.items()}
        return data

    def asdict(self, values_only=False):
        d = {}
        for ptype in self._all:
            paramdict = self._gettype(ptype)._asdict()
            if values_only:
                d[ptype] = {k: v["value"] for k, v in paramdict.items()}
            else:
                d[ptype] = paramdict
        return d

    def _named_tuple_from_dict(self, name, data):
        return namedtuple(name.upper(), data.keys())(*data.values())

    @staticmethod
    def _parse_primer3_docs(docstr):
        """
        Parse the docs of the primer3 website
        (https://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/primer3/src
        /libprimer3/primer3_manual.htm#globalTags)

        :param docstr: doc string
        :type docstr: basestring
        :return: params
        :rtype: dict
        """

        params = {}
        pattern = "(?P<name>\w+)\s+\((?P<type>[\w\s\"]+)\;\s+default\s+(?P<default>.+)\)\n(?P<description>.+)\n\n"

        type_dict = {
            "size range list": list,
            "string": str,
            "interval list": list,
            "nucleotide sequence": str,
            "int": int,
            "space separated integers": list,
            "float": float,
            "ambiguous nucleotide sequence": str,
            "boolean": bool,
            "semicolon separated list of integer \"quadruples\"": list,
            "semicolon separated list of integer quadruples": list
        }

        for m in re.finditer(pattern, docstr):
            data = m.groupdict()
            data["type_raw"] = data["type"]
            data["type"] = type_dict[data["type"]]

            default = data["default"]
            data["default_raw"] = default
            if data["type"] is str:
                if default == "empty":
                    default = ""
                default = str(default)
            elif data["type"] is list:
                if default == "empty":
                    default = []
                elif re.match("(\d+)-(\d+)", default):
                    list_match = re.match("(\d+)-(\d+)", default)
                    default = [int(list_match.group(1)), int(list_match.group(2))]
                default = list(default)
            elif data["type"] is bool:
                default = bool(int(default))
            else:
                default = data["type"](default)

            data["default"] = default
            data["value"] = default
            name = data["name"]
            params[name] = data

        return params

    @classmethod
    def _open_primer3_params(cls, filepath=None):
        if filepath is None:
            filepath = cls.default_file_path
        with open(filepath, "r") as f:
            params_txt = f.read()
            params = cls._parse_primer3_docs(params_txt)
        sequence_args = {}
        global_args = {}
        program_args = {}
        other_args = {}

        for k, v in params.items():
            if k.startswith("SEQUENCE"):
                sequence_args[k] = v
            elif k.startswith("PRIMER"):
                global_args[k] = v
            elif k.startswith("P3"):
                program_args[k] = v
            else:
                raise Exception("Parameter {} not recognized".format(k))

        return {
            "sequence": sequence_args,
            "global": global_args,
            "program": program_args,
            "other": other_args,
        }

    @staticmethod
    def _clean_params(params: dict) -> dict:
        """
        Removes empty lists and empty strings from params
        :return:
        :rtype:
        """
        cleaned = dict(params)
        ignore = ["SEQUENCE_ID"]
        for k in params:
            if k not in ignore:
                v = params[k]
                if hasattr(v, "__len__") and len(v) == 0:
                    cleaned.pop(k)
        return cleaned

    @staticmethod
    def param_values(param_dict, key="value"):
        return {k: v[key] for k, v in param_dict.items()}

    def df(self) -> pd.DataFrame:
        rows = []
        for param_type, params in self.asdict().items():
            for name, param in params.items():
                row = OrderedDict()
                row["name"] = name
                row["param_type"] = param_type
                row["value"] = param["value"]
                row["default"] = param["default"]
                row["type"] = str(param["type"])
                row["type_raw"] = param["type_raw"]
                row["description"] = param["description"]
                rows.append(row)
        return pd.DataFrame(rows, columns=rows[0].keys())

    def __add__(self, data: dict) -> dict:
        copied = self.copy()
        for k, v in data.items():
            if k in self:
                copied[k] = self[k] + v

    def copy(self):
        return deepcopy(self)


def dict_diff(d1, d2):
    diff = {}
    for k1, v1 in d1.items():
        if k1 in d2:
            v2 = d2[k1]
            if v1 != v2:
                diff[k1] = (v1, v2)
    return diff


PARAMS = Primer3Params()

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


class Primer3Design(object):
    # {param: (delta, min, max)
    DEFAULT_GRADIENT = dict(
        PRIMER_MAX_SIZE=(1, PARAMS["PRIMER_MAX_SIZE"], 36),
        PRIMER_MIN_SIZE=(-1, 16, PARAMS["PRIMER_MAX_SIZE"]),
        PRIMER_MAX_TM=(1, PARAMS["PRIMER_MAX_SIZE"], 80),
        PRIMER_MIN_TM=(-1, 48, PARAMS["PRIMER_MIN_TM"]),
        PRIMER_MAX_HAIRPIN_TH=(1, PARAMS["PRIMER_MAX_HAIRPIN_TH"], 60),
    )

    def __init__(self):
        self.params = Primer3Params()
        self.iterations = 5
        self.logger = logger(self)

    @staticmethod
    def _parse_primer3_results(results_dict):
        """
        Parse the primer3 results.

        :param results_dict:
        :type results_dict:
        :return:
        :rtype:
        """
        num_pairs = results_dict["PRIMER_PAIR_NUM_RETURNED"]
        num_left = results_dict["PRIMER_LEFT_NUM_RETURNED"]
        num_right = results_dict["PRIMER_RIGHT_NUM_RETURNED"]

        pairs = {}
        other = {}
        for i in range(max([num_pairs, num_left, num_right])):
            pairs.setdefault(i, {})
        key_pattern = "PRIMER_(?P<label>[a-zA-Z]+)_(?P<pair_id>\d+)_(?P<key>.+)"
        location_pattern = "PRIMER_(?P<label>[a-zA-Z]+)_(?P<pair_id>\d+)\s*$"
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

    @classmethod
    def _design_from_params(cls, params: Primer3Params, parse=True) -> list:
        param_dict = params.valuedict(clean=True)
        results = primer3.bindings.designPrimers(
            param_dict["SEQUENCE"], param_dict["GLOBAL"]
        )
        if parse:
            pairs, other = cls._parse_primer3_results(results)
            return pairs, other
        return results

    def design_from_params(
        self, params: Primer3Params, parse=True, max_iterations=None, gradient=None
    ) -> list:
        if gradient is None:
            gradient = {}
        gradient_dict = self.DEFAULT_GRADIENT
        gradient_dict.update(gradient)

        for i in self.logger.tqdm(
            range(max_iterations), "INFO", desc="Design iteration"
        ):
            pairs, other = self._design_from_params(params, parse)
            if pairs:
                return pairs, other
            else:
                self.apply_gradient(params, gradient_dict=gradient_dict)
        return pairs, other

    def run(self, parse=True, max_iterations=None, gradient=None) -> list:
        if max_iterations is None:
            max_iterations = self.iterations
        return self.design_from_params(self.params, parse, max_iterations, gradient)

    def design(
        self,
        new_seq_args: dict,
        new_global_args: dict,
        parse=True,
        return_with_params=False,
    ) -> list:
        params = self.params.copy()
        params.update(new_seq_args)
        params.update(new_global_args)
        results = self._design_from_params(params, parse)
        if return_with_params:
            return results, params
        return results

    def apply_gradient(cls, params, gradient_dict=None):
        if gradient_dict is None:
            gradient_dict = cls.DEFAULT_GRADIENT

        update = {}
        for param_key, gradient_tuple in gradient_dict.items():
            delta, mn, mx = gradient_tuple
            val = params[param_key] + delta
            # perform clip
            update[param_key] = max(min(val, mx), mn)
        params.update(update)
        return params

    #
    # def design_with_relaxation(self, params, max_iterations=15, gradient=None):
    #     if gradient is None:
    #         gradient = {}
    #     gradient_dict = self.DEFAULT_GRADIENT
    #     gradient_dict.update(gradient)
    #
    #     for i in range(max_iterations):
    #         pairs, other = self._design_from_params(params)
    #         if pairs:
    #             return pairs, other
    #         else:
    #             self.apply_gradient(params, gradient_dict=gradient_dict)
    #     return pairs, other
    #
    # def design_with_relaxation(self, new_seq_args: dict, new_global_args: dict, max_iterations=15, gradient=None):
    #     params = self.params.copy()
    #     params.update(new_seq_args)
    #     params.update(new_global_args)
    #
    #     if gradient is None:
    #         gradient = {}
    #     gradient_dict = self.DEFAULT_GRADIENT
    #     gradient_dict.update(gradient)
    #
    #     for i in range(max_iterations):
    #         pairs, other = self._design_from_params(params)
    #         if pairs:
    #             return pairs, other
    #         else:
    #             self.apply_gradient(params, gradient_dict=gradient_dict)
    #     return pairs, other

    def set_iterations(self, n):
        self.iterations = n

    def set_template(self, template):
        """
        The sequence from which to choose primers. The sequence must be presented 5' -> 3' (i.e, in the normal way).
        In general, the bases may be upper or lower case, but lower case letters are treated specially if PRIMER_
        LOWERCASE_MASKING is set. The entire sequence MUST be all on a single line. (In other words, the sequence cannot
        span several lines.)

        :param template:
        :type template:
        :return:
        :rtype:
        """
        self.params.update({"SEQUENCE_TEMPLATE": template})
        return self

    def set_task(self, task: str):
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
        self.params.update({"PRIMER_TASK": task})
        return self

    def set_cloning_task(self):
        return self.set_task("pick_cloning_primers")

    # TODO: set_iterations, set_num_return, set_force_return, set_gradient

    def set(self, update: dict):
        self.params.update(update)
        return self

    def set_num_return(self, n):
        return self.set({"PRIMER_NUM_RETURN": n})

    def set_size(self, interval: tuple, opt=None):
        if opt is None:
            opt = int(sum(interval) / 2.0)
        return self.set(
            {"PRIMER_PRODUCT_SIZE_RANGE": interval, "PRIMER_PRODUCT_OPT_SIZE": opt}
        )

    def set_left(self, primer: str):
        """
        The sequence of a left primer to check and around which to design right primers and optional internal
        oligos. Must be a substring of SEQUENCE_TEMPLATE.

        :param primer:
        :type primer:
        :return:
        :rtype:
        """
        return self.set({"SEQUENCE_PRIMER": primer, "PRIMER_PICK_RIGHT_PRIMER": 1})

    def set_left_only(self):
        return self.set({"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0})

    def set_right_only(self):
        return self.set({"PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1})

    def set_right(self, primer: str):
        """
        The sequence of a right primer to check and around which to design left primers and optional internal
        oligos. Must be a substring of the reverse strand of SEQUENCE_TEMPLATE.

        :param primer:
        :type primer:
        :return:
        :rtype:
        """
        return self.set(
            {"SEQUENCE_PRIMER_REVCOMP": primer, "PRIMER_PICK_LEFT_PRIMER": 1}
        )

    def set_internal(self, primer: str):
        """
        The sequence of an internal oligo to check and around which to design left and right primers. Must be a
        substring of SEQUENCE_TEMPLATE.

        :param primer:
        :type primer:
        :return:
        :rtype:
        """

    def set_primers(self, p1: str, p2: str):
        if p1:
            self.set({"SEQUENCE_PRIMER": p1, "PRIMER_PICK_LEFT_PRIMER": 1})
        if p2:
            self.set({"SEQUENCE_PRIMER_REVCOMP": p2, "PRIMER_PICK_RIGHT_PRIMER": 1})
        return self

    def set_included(self, interval: tuple):
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
        return self.set({"SEQUENCE_INCLUDED_REGION": interval})

    def set_target(self, interval: tuple):
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
        return self.set({"SEQUENCE_TARGET": interval})

    def set_target_from_template(self, template, target):
        if isinstance(target, str):
            matches = self._get_index_of_match(template, target)
            if not matches:
                print("Target not in template")
                return None
            if len(matches) > 1:
                print("More than one target found")
                return None
            return matches[0]

    def set_excluded(self, interval):
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
        return self.set({"SEQUENCE_EXCLUDED_REGION": interval})

    def set_pick_anyway(self, b=1):
        """
        If true use primer provided in SEQUENCE_PRIMER, SEQUENCE_PRIMER_REVCOMP, or SEQUENCE_INTERNAL_OLIGO even if
        it violates specific constraints.

        :param b:
        :type b:
        :return:
        :rtype:
        """
        return self.set({"PRIMER_PICK_ANYWAY": b})

    def pick_cloning_primers(
        self, seq, addnl_params={}, max_iterations=None, gradient=None
    ):
        """
        Design a left primer given the template and right primer sequences.

        :param seq: seq to pick primers from. Primers will be picked on the ends of the sequence.
        :type seq: basestring
        :param addnl_params: additional design parameters
        :type addnl_params: dict
        :param max_iterations: maximum number of iterations to use for the relaxation procedure
        :type max_iterations: int
        :param gradient: gradient dictionary to use for the relaxation procedure
        :type gradient: dict
        :return: tuple of designed primers (dict) and other parameters (dict)
        :rtype: tuple
        """
        return (
            self.copy()
            .set_template(seq)
            .set_included([0, len(seq)])
            .set_size([len(seq), len(seq)])
            .set(addnl_params)
            .set_task("pick_cloning_primers")
            .run(max_iterations=max_iterations, gradient=gradient)
        )

    def pick_pcr_primers(
        self,
        template,
        include_region,
        size_range,
        opt_size=None,
        addnl_params={},
        max_iterations=None,
        gradient=None,
    ):
        """
        Design a left and right primers given the template sequence.

        :param template: template to pick primers from
        :type template: basestring
        :param include_region: A sub-region of the given sequence in which to pick primers.
        :type include_region: list | tuple
        :param size_range: product size range
        :type size_range: list | tuple
        :param opt_size: optimal produce size
        :type opt_size: int
        :param addnl_params: additional design parameters
        :type addnl_params: dict
        :param max_iterations: maximum number of iterations to use for the relaxation procedure
        :type max_iterations: int
        :param gradient: gradient dictionary to use for the relaxation procedure
        :type gradient: dict
        :return: tuple of designed primers (dict) and other parameters (dict)
        :rtype: tuple
        """
        return (
            self.copy()
            .set_template(template)
            .set_included(include_region)
            .set_size(size_range, opt_size)
            .set(addnl_params)
            .set_task("generic")
            .run(max_iterations=max_iterations, gradient=gradient)
        )

    @dispatch_iterable([(2, (list, tuple))])
    def pick_left_pcr_primer(
        self,
        template,
        primer,
        include_region,
        size_range,
        opt_size=None,
        addnl_params={},
        max_iterations=None,
        gradient=None,
    ) -> tuple:
        """
        Design a left primer given the template and right primer sequences.

        :param template: template to pick primers from
        :type template: basestring
        :param primer: sequence of the left primer
        :type primer: basestring
        :param include_region: A sub-region of the given sequence in which to pick primers.
        :type include_region: list | tuple
        :param size_range: product size range
        :type size_range: list | tuple
        :param opt_size: optimal produce size
        :type opt_size: int
        :param addnl_params: additional design parameters
        :type addnl_params: dict
        :param max_iterations: maximum number of iterations to use for the relaxation procedure
        :type max_iterations: int
        :param gradient: gradient dictionary to use for the relaxation procedure
        :type gradient: dict
        :return: tuple of designed primers (dict) and other parameters (dict)
        :rtype: tuple
        """
        return (
            self.copy()
            .set_template(template)
            .set_included(include_region)
            .set_size(size_range, opt_size)
            .set_right(primer)
            .set(addnl_params)
            .set_task("generic")
            .run(max_iterations=max_iterations, gradient=gradient)
        )

    # TODO: add target to pick pcr primers

    @dispatch_iterable([(2, (list, tuple))])
    def pick_right_pcr_primer(
        self,
        template,
        primer,
        include_region=(),
        size_range=(),
        opt_size=None,
        addnl_params={},
        max_iterations=None,
        gradient=None,
    ):
        """
        Design a left primer given the template and right primer sequences.

        :param template: template to pick primers from
        :type template: basestring
        :param primer: sequence of the right primer
        :type primer: basestring
        :param include_region: A sub-region of the given sequence in which to pick primers.
        :type include_region: list | tuple
        :param size_range: product size range
        :type size_range: list | tuple
        :param opt_size: optimal produce size
        :type opt_size: int
        :param addnl_params: additional design parameters
        :type addnl_params: dict
        :param max_iterations: maximum number of iterations to use for the relaxation procedure
        :type max_iterations: int
        :param gradient: gradient dictionary to use for the relaxation procedure
        :type gradient: dict
        :return: tuple of designed primers (dict) and other parameters (dict)
        :rtype: tuple
        """
        return (
            self.copy()
            .set_template(template)
            .set_included(include_region)
            .set_size(size_range, opt_size)
            .set_left(primer)
            .set(addnl_params)
            .set_task("generic")
            .run(max_iterations=max_iterations, gradient=gradient)
        )

    @dispatch_iterable([(2, (list, tuple)), (3, (list, tuple))])
    def check_pcr_primers(
        self,
        template,
        p1,
        p2,
        include_region=(),
        size_range=(),
        opt_size=None,
        target=None,
        addnl_params={},
        max_iterations=None,
        gradient=None,
    ):
        design = (
            self.copy()
            .set_template(template)
            .set_included(include_region)
            .set_size(size_range, opt_size)
            .set_left(p1)
            .set_right(p2)
        )

        if target:
            if isinstance(target, str):
                design.set_target_from_template(template, target)
            else:
                design.set_target(target)
        return (
            design.set(addnl_params)
            .set_task("generic")
            .run(max_iterations=max_iterations, gradient=gradient)
        )

    def pick_sequencing_primers(
        self,
        template,
        target,
        left_only=False,
        right_only=False,
        addnl_params={},
        max_iterations=None,
        gradient=None,
    ):
        """
        Pick sequencing primers for a template given a target

        :param template: template sequence
        :type template: basestring
        :param target: target for the sequencing primers. Either a substring of the template sequence or a interval list
        :type target: basestring | list
        :param addnl_params: additional design parameters
        :type addnl_params: dict
        :param max_iterations: maximum number of iterations to use for the relaxation procedure
        :type max_iterations: int
        :param gradient: gradient dictionary to use for the relaxation procedure
        :type gradient: dict
        :return: tuple of designed primers (dict) and other parameters (dict)
        :rtype: tuple
        """
        d = self.copy().set_task("pick_sequencing_primers")
        d.set({"PRIMER_SEQUENCING_LEAD": 50})
        d.set_template(template)
        if isinstance(target, str):
            d.set_target_from_template(template, target)
        else:
            d.set_target(target)
        if right_only:
            d.set_right_only()
        elif left_only:
            d.set_left_only()
        d.set(addnl_params)
        return d.run(max_iterations=max_iterations, gradient=gradient)

    @dispatch_iterable([(1, (tuple, list)), (2, (tuple, list))])
    def check_primers(self, p1, p2, template="", addnl_params={}):
        """
        Check the penalty score for a primer pair against an optional template.

        :param p1: primer sequence
        :type p1: basestring
        :param p2: primer sequence
        :type p2: basestring
        :param template: optional template sequence
        :type template: basestring
        :return: tuple of designed primers (dict) and other parameters (dict)
        :rtype: tuple
        """
        d = self.copy().set_task("check_primers")
        d.set_primers(p1, p2)
        d.set_pick_anyway()
        if template:
            d.set_template(template)
        d.set(addnl_params)
        return d.run()

    def _get_index_of_match(self, template, sequence):
        matches = []
        for m in re.finditer(sequence, template, re.IGNORECASE):
            matches.append((m.start(0), m.end(0)))
        return matches

    # def pick_cloning_primers(self, seq, addnl_params={}, max_iterations=1, gradient=None):
    #     params = self.params.copy()
    #     update_dict = {
    #         'PRIMER_TASK': 'pick_cloning_primers',
    #         'SEQUENCE_TEMPLATE': seq,
    #         'SEQUENCE_INCLUDED_REGION': [0, len(seq)],
    #         'PRIMER_PRODUCT_SIZE_RANGE': [len(seq), len(seq)]
    #     }
    #     update_dict.update(addnl_params)
    #     params.update(update_dict)
    #     return self.design_from_params(params, max_iterations=max_iterations, gradient=gradient)
    #
    # def pick_pcr_primers(self, template, include_region, size_range, opt_size=None,
    #                      addnl_params={}, max_iterations=1, gradient=None):
    #     params = self.params.copy()
    #     size_range = list(size_range)
    #     if opt_size is None:
    #         opt_size = int(sum(size_range)/2.0)
    #     update_dict = {
    #         'PRIMER_TASK': 'pick_pcr_primers',
    #         'SEQUENCE_TEMPLATE': template,
    #         'SEQUENCE_INCLUDED_REGION': list(include_region),
    #         'PRIMER_PRODUCT_SIZE_RANGE': size_range,
    #         'PRIMER_PRODUCT_OPT_SIZE': opt_size
    #     }
    #     update_dict.update(addnl_params)
    #     params.update(update_dict)
    #     return self.design_from_params(params, max_iterations=max_iterations, gradient=gradient)
    #
    #
    # def pick_right_pcr_primer(self, template, left_primer, include_region, size_range, opt_size=None,
    #                      addnl_params={}, max_iterations=1, gradient=None):
    #     params = self.params.copy()
    #     size_range = list(size_range)
    #     if opt_size is None:
    #         opt_size = int(sum(size_range)/2.0)
    #     update_dict = {
    #         'PRIMER_TASK': 'generic',
    #         'PRIMER_PICK_LEFT_PRIMER': 0,
    #         'PRIMER_PICK_RIGHT_PRIMER': 1,
    #         'PRIMER_PICK_INTERNAL_OLIGO': 0,
    #         'SEQUENCE_TEMPLATE': template,
    #         'SEQUENCE_PRIMER': left_primer,
    #         'SEQUENCE_INCLUDED_REGION': list(include_region),
    #         'PRIMER_PRODUCT_SIZE_RANGE': size_range,
    #         'PRIMER_PRODUCT_OPT_SIZE': opt_size
    #     }
    #     update_dict.update(addnl_params)
    #     params.update(update_dict)
    #     return self.design_from_params(params, max_iterations=max_iterations, gradient=gradient)
    #
    # def pick_left_pcr_primer(self, template, right_primer, include_region, size_range, opt_size=None,
    #                           addnl_params={}, max_iterations=1, gradient=None):
    #     params = self.params.copy()
    #     size_range = list(size_range)
    #     if opt_size is None:
    #         opt_size = int(sum(size_range)/2.0)
    #     update_dict = {
    #         'PRIMER_TASK': 'generic',
    #         'PRIMER_PICK_LEFT_PRIMER': 1,
    #         'PRIMER_PICK_RIGHT_PRIMER': 0,
    #         'PRIMER_PICK_INTERNAL_OLIGO': 0,
    #         'SEQUENCE_TEMPLATE': template,
    #         'SEQUENCE_PRIMER_REVCOMP': right_primer,
    #         'SEQUENCE_INCLUDED_REGION': list(include_region),
    #         'PRIMER_PRODUCT_SIZE_RANGE': size_range,
    #         'PRIMER_PRODUCT_OPT_SIZE': opt_size
    #     }
    #     update_dict.update(addnl_params)
    #     params.update(update_dict)
    #     return self.design_from_params(params, max_iterations=max_iterations, gradient=gradient)

    # TODO: bayesian optimization for PENALTY vs results

    @staticmethod
    def open_help():
        webhelp = (
            "https://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/"
            + "primer3/src/libprimer3/primer3_manual.htm"
        )
        webbrowser.open(webhelp)

    def copy(self):
        return deepcopy(self)
