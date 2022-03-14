import logging
import sys

neat_log = logging.getLogger(__name__)


def premature_exit(err_code: int):
    print("Quitting NEAT...")
    print('-------------------------------------------------------------------------\n')
    sys.exit(err_code)


def print_lines(message: list, error_type: str):
    """
    Helper function for print_and_log
    """
    for line in message:
        print(f'{error_type.upper()} - {line}')


def log_lines(message: list, error_type: str):
    """
    Helper function for print_and_log7
    """

    if error_type.lower() == 'critical':
        for line in message:
            neat_log.critical(line)
    elif error_type.lower() == 'error':
        for line in message:
            neat_log.error(line)
    elif error_type.lower() == 'warning':
        for line in message:
            neat_log.warning(line)
    elif error_type.lower() == 'info':
        for line in message:
            neat_log.info(line)
    elif error_type.lower() == 'debug':
        for line in message:
            neat_log.debug(line)
    else:
        log_mssg(f'BUG: Unknown error type', 'critical')
        premature_exit(1)


def log_mssg(mssg: str, error_type: str):
    """
    Prints and logs any message. Only for info level messages

    :param mssg: a string with the message.
    :param error_type: A string that indicates the level of error. Has to be a declared type for the logging module
    """
    message = mssg.split('\n')
    message[1:] = [f'\t - {x}' for x in message[1:]]
    log_lines(message, error_type)
