from importlib.metadata import PackageNotFoundError, version as _version

def neat_version() -> str:
    """
    Return NEAT's package version.
    """
    try:
        return _version("neat")
    except PackageNotFoundError:
        return "unknown"
