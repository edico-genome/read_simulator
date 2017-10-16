import sys
import json
import yaml
import logging

logger = logging.getLogger(__name__)


class Settings(object):
    def __init__(self, run_settings_file):
        self.runs = self.load_config(run_settings_file)

    @staticmethod
    def load_config(f_path):
        try:
            with open(f_path) as stream:
                j = json.load(stream)
            return j
        except IOError:
            print("Failed to open file: {}".format(f_path))
            sys.exit(1)
        except Exception as e:
            print("Failed to import file: {} as json. Exception: {}".format(f_path, e))
            sys.exit(1)

    def log_yaml(self, y_dict):
        y = yaml.safe_dump(
            y_dict, indent=4,
            allow_unicode=True,
            default_flow_style=False)
        logger.info(y)

    def print_settings(self):
        for idx, run in enumerate(self.runs):
            logger.info("RUN {}".format(idx+1))
            self.log_yaml(run)
