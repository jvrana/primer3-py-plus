from primer3plus.design.interfaces import ParameterDescriptor, ParameterAccessor
from primer3plus.params import BoulderIO, default_boulderio


def test_descriptor():
    class Foo(object):

        SEQUENCE_ID = ParameterDescriptor("SEQUENCE_ID")

        def __init__(self):
            self.params = default_boulderio()

    foo = Foo()
    foo.params["SEQUENCE_ID"] = "bar"
    assert foo.SEQUENCE_ID == "bar"

    foo.SEQUENCE_ID = "baz"
    assert foo.SEQUENCE_ID == "baz"
    assert foo.params["SEQUENCE_ID"] == "baz"


def test_accessor():
    class Foo(object):

        P = ParameterAccessor()

        def __init__(self):
            self.params = default_boulderio()

    foo = Foo()
    foo2 = Foo()
    foo.params["SEQUENCE_ID"] = "bar"
    assert foo.P.SEQUENCE_ID == "bar"

    foo.P.SEQUENCE_ID = "baz"
    assert foo.P.SEQUENCE_ID == "baz"
    assert foo.params["SEQUENCE_ID"] == "baz"
    assert not foo2.params["SEQUENCE_ID"] == "baz"
    assert not foo2.P.SEQUENCE_ID == "baz"
