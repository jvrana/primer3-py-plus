from __future__ import annotations

import os
import re
import webbrowser
from collections import Mapping
from copy import deepcopy
from typing import Any
from typing import Dict
from typing import Iterator

from .expected_opts import _expected_opts
from primer3plus.constants import DOCURL


class ParamTypes:
    GLOBAL = "GLOBAL"
    PROGRAM = "PROGRAM"
    OTHER = "OTHER"
    SEQUENCE = "SEQUENCE"
    CATEGORY = "category"


REL_PARAM_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "primer3_params_raw.txt"
)


class ParameterType:
    def __init__(self, name, description, type, default, category):
        """
        Initialize a BoulderIO parameter type.

        :param name: name of the parameter
        :param description: description of the parameter
        :param type: expected type of the parameter
        :param default: default value of the parameter
        :param category: parameter category (from 'sequence', 'program', 'globals', or
                'other')
        """
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


class Parameter:
    def __init__(self, ptype: ParameterType, value=None):
        """
        Initialize a parameter from a
        :class:`Parameter <primer3plus.params.ParameterType>`

        :param ptype: parameter type
        :param value: value to set parameter
        """
        self.ptype = ptype
        self._value = None
        if value is None:
            self.value = self.ptype.default
        else:
            self.value = value

    @property
    def value(self) -> Any:
        """
        Return the value of the parameter.

        :return: Any
        """
        return self._value

    @value.setter
    def value(self, v: Any):
        """
        Set the value of the parameter

        :param v: the value
        :return: None
        """
        if not issubclass(self.ptype.type, type(v)):
            raise TypeError(
                "Paramater {ptype} must be of type {type} not {nottype}".format(
                    ptype=self.ptype, type=self.ptype.type, nottype=type(v)
                )
            )
        self._value = self.ptype.type(v)

    @property
    def name(self):
        """
        Return the name of the parameter

        :return: parameter name
        """
        return self.ptype.name

    def set_default(self):
        """
        Set the parameter to its default.

        :return: None
        """
        self.value = self.ptype.default

    def copy(self) -> Parameter:
        """
        Make a copy of this parameter.

        :return: None
        """
        p = self.__class__(self.ptype, deepcopy(self.value))
        return p

    def help(self, open=False):
        return "{}#{}".format(DOCURL, self.name)

    def __str__(self) -> str:
        return "<{cls} {name} {type} value={value} default={default}>".format(
            cls=self.__class__.__name__,
            name=self.name,
            type=self.ptype.type,
            value=self.value,
            default=self.ptype.default,
        )

    def __repr__(self) -> str:
        return str(self)


class BoulderIO(Mapping):
    POST_LOAD_DEFAULTS = {"PRIMER_EXPLAIN_FLAG": 1}
    EXPECTED = _expected_opts[:]

    def __init__(self):
        self._params = {}

    def update(self, data_dict: Dict[str, Any]):
        """
        Update the parameters from a dictionary of key:values.

        :param data_dict: update dictionary
        :return: None
        """
        for k, v in data_dict.items():
            self[k] = v

    def _post_load(self):
        self.update(self.POST_LOAD_DEFAULTS)

    def load(self, param_dict):
        """
        Load parameters from a dictionary.

        :param param_dict:
        :return:
        """
        for k, v in param_dict.items():
            ptype = ParameterType(
                name=v["name"],
                type=v["type"],
                default=v["default"],
                description=v["description"],
                category=v[ParamTypes.CATEGORY],
            )
            p = Parameter(ptype, ptype.default)
            self._params[p.name] = p
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

    @staticmethod
    def _raise_no_key(key):
        return KeyError(
            "{key} not in params. See docs for help: {url}".format(key=key, url=DOCURL)
        )

    @staticmethod
    def online_help(open=False, key=None) -> str:
        """
        Display online help in a browser tab.

        :param open: if True, open browser tab. Else return url.
        :param key: optional parameter key
        :return:
        """
        if key:
            url = "{url}#{key}".format(url=DOCURL, key=key)
        else:
            url = DOCURL
        if open:
            webbrowser.open(url)
        return url

    def _by_category(self, category):
        return {
            k: v.value for k, v in self._params.items() if v.ptype.category == category
        }

    def _globals(self, clean=True) -> Dict[str, Any]:
        """
        Return global parameters.

        :param clean: if True, will remove empty lists and empty strings from params.
        :return: parameter values as a dict.
        """
        data = self._by_category(ParamTypes.GLOBAL)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def _program(self, clean=True):
        """
        Return program parameters.

        :param clean: if True, will remove empty lists and empty strings from params.
        :return: parameter values as a dict.
        """
        data = self._by_category(ParamTypes.PROGRAM)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def _sequence(self, clean=True):
        """
        Return sequence parameters.

        :param clean: if True, will remove empty lists and empty strings from params.
        :return: parameter values as a dict.
        """
        data = self._by_category(ParamTypes.SEQUENCE)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def _other(self, clean=True):
        """
        Return other parameters.

        :param clean: if True, will remove empty lists and empty strings from params.
        :return: parameter values as a dict.
        """
        data = self._by_category(ParamTypes.OTHER)
        if clean:
            data = self._clean_dictionary(data)
        return data

    def all(self):
        return {k: v for k, v in self.items()}

    def set_defaults(self):
        """
        Set all parameters to their defaults.

        :return:
        """
        for v in self._params.values():
            v.set_default()

    def values(self) -> Iterator[Any]:
        """
        Iterator for parameter values.

        :return: iterator over parameter values
        """
        for v in self._params.values():
            yield v.value

    def items(self):
        for k, v in self._params.items():
            yield k, v.value

    def copy(self):
        copied = BoulderIO()
        for k, v in self._params.items():
            copied._params[k] = v.copy()
        return copied

    @staticmethod
    def _clean_dictionary(params: dict) -> dict:
        """Removes empty lists and empty strings from params.

        :return: :rtype:
        """
        cleaned = dict(params)
        ignore = ["SEQUENCE_ID"]
        for k in params:
            if k not in ignore:
                v = params[k]
                if hasattr(v, "__len__") and len(v) == 0:
                    cleaned.pop(k)
        return cleaned

    @property
    def defs(self):
        return dict(self._params)

    def __contains__(self, key: str) -> bool:
        return key in self._params

    def __setitem__(self, key: str, value: Any):
        try:
            self._params[key].value = value
        except KeyError:
            raise self._raise_no_key(key)

    def __getitem__(self, key: str) -> Any:
        try:
            return self._params[key].value
        except KeyError:
            raise self._raise_no_key(key)

    # def __delitem__(self, key: str):
    #     try:
    #         self._params[key].set_default()
    #     except KeyError:
    #         raise self._raise_no_key(key)

    def __len__(self) -> int:
        return len(self._params)

    def __iter__(self) -> Iterator[str]:
        for k, v in self._params.items():
            yield k

    def __str__(self):
        return str(self._params)


class ParamParser:
    """Reads the Primer3 documentation and creates the appropriate
    parameters."""

    @staticmethod
    def _parse_primer3_docs(docstr):
        """Parse the docs of the primer3 website
        (https://htmlpreview.github.io.

        /?https://github.com/libnano/primer3-py/master/primer3/src.
        /libprimer3/primer3_manual.htm#globalTags)  :param docstr: doc string :type
        docstr: basestring :return: params :rtype: dict
        """

        params = {}

        expected_names = "(?P<name>{})".format("|".join(_expected_opts))

        catch_all = "^{name_pattern}(?P<rest>\\s+\\(.+)\n".format(
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

            pattern = r"\s+\((?P<type>.+?);\s+default\s+(?P<default>.+)\)"
            m = re.match(pattern, rest)
            if not m:
                raise Exception("Did not catch:\n{}".format(caught.group(0)))

            type_str = m.groupdict()["type"]
            type_str = re.match(r"([a-zA-z\s]+)", type_str).group(1).strip()
            default = m.groupdict()["default"]
            ptype = type_dict[type_str]

            if ptype is str:
                if default == "empty":
                    default = ""
                default = str(default)
            elif ptype is list:
                if default == "empty":
                    default = []
                elif re.match(r"(\d+)-(\d+)", default):
                    list_match = re.match(r"(\d+)-(\d+)", default)
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
    def _open_primer3_params(cls, filepath=None) -> BoulderIO:
        if filepath is None:
            filepath = REL_PARAM_PATH
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
    def open(cls, filepath=None) -> BoulderIO:
        """
        Open the parameters from the filepath.

        :param filepath: filepath
        :return:
        """
        return cls._open_primer3_params(filepath)


def default_boulderio() -> BoulderIO:
    """
    Open the default parameters as a :class:`BoulderIO <primer3plus.params.params>`

    :return: the BoulderIO instance.
    """
    param_dict = ParamParser.open()
    boulderio = BoulderIO()
    boulderio.load(param_dict)
    return boulderio
