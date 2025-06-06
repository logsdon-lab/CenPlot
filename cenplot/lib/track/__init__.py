"""
Module for track settings and types.
"""

from .types import Track, TrackType, TrackPosition, TrackList, LegendPosition
from .settings import (
    SelfIdentTrackSettings,
    BarTrackSettings,
    LabelTrackSettings,
    HORTrackSettings,
    HOROrtTrackSettings,
    PositionTrackSettings,
    LegendTrackSettings,
    SpacerTrackSettings,
    LocalSelfIdentTrackSettings,
    StrandTrackSettings,
    TrackSettings,
)

__all__ = [
    "Track",
    "TrackType",
    "TrackPosition",
    "TrackList",
    "LegendPosition",
    "SelfIdentTrackSettings",
    "LocalSelfIdentTrackSettings",
    "BarTrackSettings",
    "LabelTrackSettings",
    "HORTrackSettings",
    "HOROrtTrackSettings",
    "PositionTrackSettings",
    "StrandTrackSettings",
    "LegendTrackSettings",
    "SpacerTrackSettings",
    "TrackSettings",
]
