import os
import logging
from datetime import datetime

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = os.path.join(dir_path, '..', 'logs')

if not os.path.isdir(dir_path):
    os.makedirs(dir_path)

"""
log_file_detail = datetime.today()
log_file_detail = str(log_file_detail).split('.')[0].replace(
    '-', '_').replace(' ', '__').replace(':', '_')
log_file = os.path.join(dir_path, "simulator_{}.log".format(log_file_detail))
"""

log_file = os.path.join(dir_path, "simulator.log")
print('logging to {}'.format(log_file))

logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    filemode='w',
    format='%(asctime)s %(message)s',
    datefmt='%m/%d - %I:%M:%S')

logging.getLogger().addHandler(logging.StreamHandler())


