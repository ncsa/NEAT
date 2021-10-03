import logging
import sys


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
            logging.critical(line)
    elif error_type.lower() == 'error':
        for line in message:
            logging.error(line)
    elif error_type.lower() == 'warning':
        for line in message:
            logging.warning(line)
    elif error_type.lower() == 'info':
        for line in message:
            logging.info(line)
    elif error_type.lower() == 'debug':
        for line in message:
            logging.debug(line)
    else:
        print_and_log(f'BUG: Unknown error type', 'critical')
        premature_exit(1)


def print_and_log(mssg: str, error_type: str):
    """
    Prints and logs any message. Only for info level messages

    :param mssg: a string with the message.
    :param error_type: A string that indicates the level of error. Has to be a declared type for the logging module
    """
    message = mssg.split('\n')
    message[1:] = [f'\t - {x}' for x in message[1:]]
    print_lines(message, error_type.lower())
    log_lines(message, error_type)
