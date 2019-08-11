import os
import re
from collections import MutableMapping
from copy import deepcopy
import webbrowser
from .expected_params import expected_opts


class ParamTypes(object):
    GLOBAL = "GLOBAL"
    PROGRAM = "PROGRAM"
    OTHER = "OTHER"
    SEQUENCE = "SEQUENCE"
    CATEGORY = "category"


class Constants(object):
    DOCURL = "https://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_manual.htm"
    REL_PARAM_PATH = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "primer3_params_raw.txt"
    )


class ParameterType(object):
    def __init__(self, name, description, type, default, category):
        self.name = name
        self.description = description
        self.type = type
        self.default = default
        self.category = category

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "<{cls} {name} {type} default={default}>".format(
            cls=self.__class__.__name__,
            name=self.name,
            type=self.type,
            default=self.default,
        )


class Parameter(ParameterType):
    def __init__(self, ptype, value=None):
        self.ptype = ptype
        self._value = None
        if value is None:
            self.value = self.ptype.default
        else:
            self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, v):
        if not issubclass(self.ptype.type, type(v)):
            raise TypeError(
                "Paramater {ptype} must be of type {type}".format(
                    ptype=self.ptype, type=self.ptype.type
                )
            )
        self._value = self.ptype.type(v)

    @property
    def name(self):
        return self.ptype.name

    def set_default(self):
        self.value = self.ptype.default

    def copy(self):
        p = self.__class__(self.ptype, deepcopy(self.value))
        return p

    def __str__(self):
        return "<{cls} {name} {type} value={value} default={default}>".format(
            cls=self.__class__.__name__,
            name=self.name,
            type=self.ptype.type,
            value=self.value,
            default=self.ptype.default,
        )

    def __repr__(self):
        return str(self)


class BoulderIO(MutableMapping):
    POST_LOAD_DEFAULTS = {"PRIMER_EXPLAIN_FLAG": 1}
    EXPECTED = expected_opts[:]

    def __init__(self):
        self.params = {}

    def update(self, data_dict):
        for k, v in data_dict.items():
            self[k] = v

    def _post_load(self) -> dict:
        self.update(self.POST_LOAD_DEFAULTS)

    def load(self, param_dict):
        for k, v in param_dict.items():
            ptype = ParameterType(
                name=v["name"],
                type=v["type"],
                default=v["default"],
                description=v["description"],
                category=v[ParamTypes.CATEGORY],
            )
            p = Parameter(ptype, ptype.default)
            self.params[p.name] = p
        self._post_load()
        missing = self._check_missing()
        if missing:
            raise ValueError(
                "The following keys are missing: {}".format(" ".join(missing))
            )

    def _check_missing(self):
        missing = []
        for key in self.EXPECTED:
            if key not in self:
                missing.append(key)
        return missing

    def _raise_no_key(self, key):
        return KeyError(
            "{key} not in params. See docs for help: {url}".format(
                key=key, url=Constants.DOCURL
            )
        )

    def __contains__(self, key):
        return key in self.params

    def __setitem__(self, key, value):
        try:
            self.params[key].value = value
        except KeyError:
            raise self._raise_no_key(key)

    def __getitem__(self, key):
        try:
            return self.params[key].value
        except KeyError:
            raise self._raise_no_key(key)

    def __delitem__(self, key):
        try:
            self.params[key].set_default()
        except KeyError:
            raise self._raise_no_key(key)

    def __len__(self):
        return len(self.params)

    def __iter__(self):
        for k, v in self.params.items():
            yield k

    def online_help(self, open=False, key=None):
        if key:
            url = "{url}#{key}".format(Constants.DOCURL, key)
        else:
            url = Constants.DOCURL
        if open:
            webbrowser.open(url)
        return url

    def _by_category(self, category):
        return {
            k: v.value for k, v in self.params.items() if v.ptype.category == category
        }

    def globals(self, clean=True):
        data = self._by_category(ParamTypes.GLOBAL)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def program(self, clean=True):
        data = self._by_category(ParamTypes.PROGRAM)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def sequence(self, clean=True):
        data = self._by_category(ParamTypes.SEQUENCE)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def other(self, clean=True):
        data = self._by_category(ParamTypes.OTHER)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def set_defaults(self):
        for v in self.params.values():
            v.set_default()

    def values(self):
        for v in self.params.values():
            yield v.value

    def items(self):
        for k, v in self.params.items():
            yield k, v.value

    def copy(self):
        copied = BoulderIO()
        for k, v in self.params.items():
            copied.params[k] = v.copy()
        return copied

    @staticmethod
    def _clean_dictionary(params: dict) -> dict:
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


class ParamParser(object):
    """
    Reads the Primer3 documentation and creates the appropriate parameters.
    """

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

        expected_names = "(?P<name>{})".format("|".join(expected_opts))

        catch_all = "^{name_pattern}(?P<rest>\s+\(.+)\n".format(
            name_pattern=expected_names
        )

        type_dict = {
            "size range list": list,
            "string": str,
            "interval list": list,
            "nucleotide sequence": str,
            "int": int,
            "space separated integers": list,
            "float": float,
            "decimal": float,
            "ambiguous nucleotide sequence": str,
            "boolean": bool,
            'semicolon separated list of integer "quadruples"': list,
            "semicolon separated list of integer quadruples": list,
        }

        for caught in re.finditer(catch_all, docstr, flags=re.MULTILINE):
            caught_dict = caught.groupdict()
            name = caught_dict["name"]
            rest = caught_dict["rest"]

            pattern = "\s+\((?P<type>.+?);\s+default\s+(?P<default>.+)\)"
            m = re.match(pattern, rest)
            if not m:
                raise Exception("Did not catch:\n{}".format(caught.group(0)))

            type_str = m.groupdict()["type"]
            type_str = re.match("([a-zA-z\s]+)", type_str).group(1).strip()
            default = m.groupdict()["default"]
            ptype = type_dict[type_str]

            if ptype is str:
                if default == "empty":
                    default = ""
                default = str(default)
            elif ptype is list:
                if default == "empty":
                    default = []
                elif re.match("(\d+)-(\d+)", default):
                    list_match = re.match("(\d+)-(\d+)", default)
                    default = [int(list_match.group(1)), int(list_match.group(2))]
                default = list(default)
            elif ptype is bool:
                default = bool(int(default))
            elif ptype is int:
                default = int(default)
            elif ptype is float:
                default = float(default)
            else:
                raise ValueError(str(caught.group(0)))

            params[name] = {
                "name": name,
                "type": ptype,
                "type_raw": type_str,
                "default": default,
                "description": "",
            }
        if not params:
            raise Exception("Could not load parameters")
        return params

    @classmethod
    def _open_primer3_params(cls, filepath=None):
        if filepath is None:
            filepath = Constants.REL_PARAM_PATH
        with open(filepath, "r") as f:
            params_txt = f.read()
            params = cls._parse_primer3_docs(params_txt)

        for k, v in params.items():
            if k.startswith("SEQUENCE"):
                param_type = ParamTypes.SEQUENCE
            elif k.startswith("PRIMER"):
                param_type = ParamTypes.GLOBAL
            elif k.startswith("P3"):
                param_type = ParamTypes.PROGRAM
            else:
                raise Exception("Parameter {} not recognized".format(k))
            v[ParamTypes.CATEGORY] = param_type

        return params

    @classmethod
    def open(cls, filepath=None):
        return cls._open_primer3_params(filepath)


def default_boulderio():
    param_dict = ParamParser.open()
    boulderio = BoulderIO()
    boulderio.load(param_dict)
    return boulderio
