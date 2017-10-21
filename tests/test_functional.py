import pytest, os
from dSim import Settings, main

this_dir_path = os.path.dirname(os.path.realpath(__file__))
test_json = os.path.join(this_dir_path, "run_test.json")

def test_basic_end2end():
    settings = Settings(test_json)
    rv = main(settings)
    assert (rv == 0)
