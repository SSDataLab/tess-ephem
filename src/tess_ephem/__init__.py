import logging

# Configure logging
log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler())

from .ephem import ephem, TessEphem  # noqa: E402

__version__ = "0.1.0"
__all__ = ["ephem", "TessEphem"]
