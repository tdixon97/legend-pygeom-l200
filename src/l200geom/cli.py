from __future__ import annotations

import argparse
import logging
import sys

from pyg4ometry import gdml

from . import _version
from .core import construct

log = logging.getLogger(__name__)


def dump_gdml_cli() -> None:
    parser = argparse.ArgumentParser(
        prog="legend-pygeom-l200",
        description="legend-pygeom-l200 command line interface",
    )

    # global options
    parser.add_argument(
        "--version", action="store_true", help="""Print pygama version and exit"""
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="""Increase the program verbosity""",
    )
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        help="""Increase the program verbosity to maximum""",
    )

    parser.add_argument(
        "filename",
        default="l200.gdml",
        help="""File name for the output GDML geometry.""",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger("l200geom").setLevel(logging.DEBUG)
    elif args.debug:
        logging.root.setLevel(logging.DEBUG)

    if args.version:
        print(version.__version__)  # noqa: T201
        sys.exit()

    log.info(f"exporting GDML geometry to {args.filename}")
    w = gdml.Writer()
    w.addDetector(construct())
    w.write(args.filename)
