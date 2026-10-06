from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("eefinder")
except PackageNotFoundError:
    __version__ = "unknown"
