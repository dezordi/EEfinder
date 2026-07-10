from importlib.metadata import PackageNotFoundError, version
import sys

try:
    __version__ = version("eefinder")
except PackageNotFoundError:
    sys.exit(1)
