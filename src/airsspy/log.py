"""Logging configuration for airsspy.

Provides a shared ``setup_logging()`` function used by the CLI and available
to library / workflow users who want timestamped log output.
"""

import logging
import sys


def setup_logging(level: int = logging.INFO) -> None:
    """Configure the root logger with timestamps.

    Args:
        level: Logging level (default ``INFO``).  Use ``logging.DEBUG`` for
            verbose output or ``logging.WARNING`` for quiet mode.
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        force=True,
    )
