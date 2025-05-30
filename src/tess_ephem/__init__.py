import logging

# Configure logging
log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler())

from .ephem import TessEphem, ephem  # noqa: E402

__version__ = "0.5.2"
__all__ = ["ephem", "TessEphem"]
