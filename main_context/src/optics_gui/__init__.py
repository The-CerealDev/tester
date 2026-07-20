"""
Packaged backend for the ISIS RCS optics GUI.
"""

try:
    from importlib.metadata import version
except ImportError:  # pragma: no cover
    version = None

if version is None:
    __version__ = "0.1.0"
else:
    try:
        __version__ = version("optics-gui")
    except Exception:  # pragma: no cover
        __version__ = "0.1.0"

__all__ = ["__version__"]
