# import libs
import logging
from typing import List, Optional, Dict, Literal, cast
from pythermodb_settings.models import Component

# NOTE: logger setup
logger = logging.getLogger(__name__)


def map_component_state(
    component: Component,
    options: Dict[str, str] = {
        "g": "IG",
        "l": "LIQ",
        "s": "SOL",
    },
    **kwargs
) -> str:
    """
    Map the phase state for a single component.

    Parameters
    ----------
    component : Component
        A Component object.
    options : Dict[str, str], optional
        A dictionary mapping component states to phase states, by default { "g": "IG", "l": "LIQ", "s": "SOL" }.
    kwargs : additional keyword arguments
        Additional keyword arguments (not used).
        - phase_state_options : Dict[str, str], optional
            A dictionary to override the default options.

    Returns
    -------
    str
        The phase state of the component.
    """
    try:
        # SECTION: set default options
        options = {
            **options, **kwargs['phase_state_options']
        } if 'phase_state_options' in kwargs else options

        # SECTION: map component state
        # NOTE: get phase state from options
        phase_state = options.get(component.state.lower(), None)

        # NOTE: phase state found
        if phase_state is None:
            raise ValueError(
                f"Phase state not found for component state: {component.state}"
            )

        return phase_state
    except Exception as e:
        logger.error(f"Error setting component state: {e}")
        raise
