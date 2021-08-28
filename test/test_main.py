import easydev
import os
import tempfile
import subprocess
import sys
import sequana_pipelines.lora.main as m

from . import test_dir


def test_standalone_subprocess(tmpdir):
    wkdir = tmpdir.mkdir("wkdir")

    cmd = "sequana_lora  --input-directory {test_dir}/data "
    cmd += f"--working-directory {wkdir} --run-mode local --force"
    subprocess.call(cmd.split())


def test_standalone_script(tmpdir):
    wkdir = tmpdir.mkdir("wkdir")
    sys.argv = ["test", "--input-directory", f"{test_dir}/data", "--working-directory",
        str(wkdir), "--run-mode", "local", "--force"]
    m.main()



def test_version():
    cmd = "sequana_lora --version"
    subprocess.call(cmd.split())
