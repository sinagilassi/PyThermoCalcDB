# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, AnnotatedValue
from pythermodb_settings.utils import (
    to_annotated_value,
)
from pythermodb_settings.decorators import calculation_info
from .core.fractions import (
    _calc_fractions,
    _calc_fractions_from_mapping,
    _calc_component_fractions,
)

# NOTE: logger set
logger = logging.getLogger(__name__)


# ======================================================================
# *** Public annotated API
# ======================================================================

# ::: annotated for numpy array

# ::: annotated for mapping

# ::: annotated for component mapping


# all
_all_ = []
