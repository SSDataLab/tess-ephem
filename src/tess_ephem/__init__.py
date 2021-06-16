import logging

# Configure logging
log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler())

from .ephem import ephem, TessEphem  # noqa: E402

__version__ = "0.4.dev"
__all__ = ["ephem", "TessEphem"]
