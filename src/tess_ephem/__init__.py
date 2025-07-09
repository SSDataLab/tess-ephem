import logging

# Configure logging
log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler())

from .ephem import TessEphem, ephem  # noqa: E402

__version__ = "0.6.1"
__all__ = ["ephem", "TessEphem"]
