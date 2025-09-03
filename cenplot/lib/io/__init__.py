"""
Module for reading track data.
"""

from .bed9 import read_bed9
from .bed_hor import read_bed_hor
from .bed_label import read_bed_label
from .bed_identity import read_bed_identity
from .tracks import read_tracks, read_track

__all__ = [
    "read_bed9",
    "read_bed_hor",
    "read_bed_label",
    "read_bed_identity",
    "read_track",
    "read_tracks",
]
