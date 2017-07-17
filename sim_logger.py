import logging
from datetime import datetime

log_file = datetime.today()
log_file = str(log_file).split('.')[0].replace(
    '-', '_').replace(' ', '__').replace(':', '_')
log_file = "simulator_{}.log".format(log_file)

logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    filemode='w',
    format='%(asctime)s %(message)s',
    datefmt='%m/%d - %I:%M:%S')

print('logging to {}'.format(log_file))
