from .__version__ import __authors__
from .__version__ import __homepage__
from .__version__ import __repo__
from .__version__ import __title__
from .__version__ import __version__
from .design import Design
from .log import logger


def new(params=None):
    """Start a new design."""
    design = Design()
    if params:
        design.params.update(params)
    return design
