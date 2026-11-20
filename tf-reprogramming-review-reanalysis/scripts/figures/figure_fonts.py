"""Locate optional system fonts without redistributing font files."""

import os
from pathlib import Path


def font_paths():
    override = os.environ.get("FIGURE_FONT_DIR")
    if override:
        return (Path(override),)
    if os.name == "nt":
        return (Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts",)
    return (Path("/mnt/c/Windows/Fonts"),)


def system_font(name):
    for folder in font_paths():
        path = folder / name
        if path.is_file():
            return path
    return None
