from primer3plus import Design
from primer3plus.design.interfaces import ParameterDescriptor
from primer3plus.params import default_boulderio


def test_descriptor():
    class Foo:

        SEQUENCE_ID = ParameterDescriptor("SEQUENCE_ID")

        def __init__(self):
            self.params = default_boulderio()

    foo = Foo()
    foo.params["SEQUENCE_ID"] = "bar"
    assert foo.SEQUENCE_ID.value == "bar"

    foo.SEQUENCE_ID.value = "baz"
    assert foo.SEQUENCE_ID.value == "baz"
    assert foo.params["SEQUENCE_ID"] == "baz"


def test_design_accessor():

    design = Design()
    design.SEQUENCE_ID = "baz"
    assert design.params["SEQUENCE_ID"] == "baz"
    # assert design.P.SEQUENCE_ID == 'baz'
