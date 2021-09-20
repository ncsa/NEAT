import logging
import sys


def will_exit(err_code: int):
    print("Quitting NEAT...")
    print('-------------------------------------------------------------------------\n')
    sys.exit(err_code)


def print_and_log(message: str, error_type: str):
    """
    Prints and logs any message. Only for info level messages

    :param message: a string with the message.
    :param error_type: A string that indicates the level of error. Has to be a declared type for the logging module
    """
    print(f'{error_type.upper()} - {message}')
    if error_type.lower() == 'critical':
        logging.critical(message)
    elif error_type.lower() == 'error':
        logging.error(message)
    elif error_type.lower() == 'warning':
        logging.warning(message)
    elif error_type.lower() == 'info':
        logging.info(message)
    elif error_type.lower() == 'debug':
        logging.debug(message)
    else:
        print_and_log(f'BUG: Unknown error type', 'critical')
        will_exit(1)
