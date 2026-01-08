try:
    from importlib.metadata import version
    __version__ = version('genbank_to')
except Exception:
    __version__ = 'unknown'
