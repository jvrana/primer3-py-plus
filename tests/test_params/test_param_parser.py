from primer3plus.params import ParamParser, BoulderIO
import pytest


@pytest.fixture(scope="function")
def params():
    param_data = ParamParser.open()

    params = BoulderIO()
    params.load(param_data)
    return params


def test_open(params):
    assert params


def test_open_and_load():
    param_data = ParamParser.open()

    params = BoulderIO()
    params.load(param_data)


class TestIter(object):
    def test_iter(self, params):
        x = list(params)
        assert len(x) == len(params)
        print(x)

    def test_items(self, params):
        x = list(params.items())
        assert len(x) == len(params)
        print(x)

    def test_value(self, params):
        x = list(params.values())
        assert len(x) == len(params)
        print(x)


class TestSettersAndGetters(object):
    def test_set_value_raises_type_error(self, params):
        with pytest.raises(TypeError):
            params["SEQUENCE_TARGET"] = 7
        params["SEQUENCE_TARGET"] = []

    def test_set_value(self, params):
        params["SEQUENCE_TARGET"] = [1, 2, 3, 4]
        assert params["SEQUENCE_TARGET"] == [1, 2, 3, 4]

    def test_del_value(self, params):
        params["SEQUENCE_TARGET"] = [1, 2, 3, 4]
        assert params["SEQUENCE_TARGET"] == [1, 2, 3, 4]
        del params["SEQUENCE_TARGET"]
        assert params["SEQUENCE_TARGET"] == []

    def test_set_default(self, params):
        params["SEQUENCE_TARGET"] = [1, 2, 3, 4]
        assert params["SEQUENCE_TARGET"] == [1, 2, 3, 4]
        params.set_defaults()
        assert params["SEQUENCE_TARGET"] == []


class TestCopy(object):
    def test_copy(self, params):

        params_copy = params.copy()
        params["SEQUENCE_TARGET"].append(4)
        assert params["SEQUENCE_TARGET"] == [4]
        assert not params_copy["SEQUENCE_TARGET"] == [4]

        assert (
            params_copy.params["SEQUENCE_TARGET"].ptype
            is params.params["SEQUENCE_TARGET"].ptype
        )


class TestCategory(object):
    def test_global(self, params):
        assert params.globals

    def test_sequence(self, params):
        assert params.sequence

    ## not implemented
    # def test_program(self, params):
    #     assert params.program
