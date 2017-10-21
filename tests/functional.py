import pytest
from dSim.dSim import Settings, main

def test_basic_end2end():
    settings = Settings("run_test.json")
    rv = main()
    assert (rv == 0)
