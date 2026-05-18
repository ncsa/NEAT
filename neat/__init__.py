from importlib.metadata import PackageNotFoundError, version

try:
    __version__: str = version("neat-genreads")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
