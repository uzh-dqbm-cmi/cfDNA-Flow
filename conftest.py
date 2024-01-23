#! /usr/bin/env python3

import pytest


def pytest_addoption(parser):
    parser.addoption("--inst_folder", action="store", default="./",
                     help="Installation folder. Default: './' , used for testing from sources, "
                          "'/' would be for testing from a singularity image")


@pytest.fixture
def inst_folder(request):
    inst_folder = request.config.getoption("--inst_folder")
    inst_folder = inst_folder if inst_folder.endswith("/") else inst_folder + "/"
    return inst_folder

# for local run:
# pytest -q --inst_folder="./" test_cfDNA.py

# for singularity image run:
# pytest -q --inst_folder="/" test_cfDNA.py
