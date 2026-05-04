import logging
import pytest


@pytest.fixture(autouse=True)
def _isolate_neat_logging():
    """
    Close and remove any FileHandlers attached to NEAT loggers before each test.
    Prevents 'ValueError: I/O operation on closed file' errors when a FileHandler
    from a previous test is still attached after its underlying file is closed.
    Propagation is left intact so caplog can capture NEAT log output.
    """
    def _close_file_handlers(logger):
        for h in list(logger.handlers):
            if isinstance(h, logging.FileHandler):
                logger.removeHandler(h)
                try:
                    h.close()
                except Exception:
                    pass

    for name, logger in list(logging.Logger.manager.loggerDict.items()):
        if (name == "neat" or name.startswith("neat.")) and isinstance(logger, logging.Logger):
            _close_file_handlers(logger)

    yield

    for name, logger in list(logging.Logger.manager.loggerDict.items()):
        if (name == "neat" or name.startswith("neat.")) and isinstance(logger, logging.Logger):
            _close_file_handlers(logger)