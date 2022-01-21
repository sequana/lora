import os
import subprocess
import sys

import sequana_pipelines.lora.main as m

from . import test_dir


def test_standalone_subprocess(tmpdir):
    input_dir = os.sep.join((test_dir, 'resources'))
    cmd = ["test", "--input-directory", input_dir, "--working-directory", str(tmpdir), "--force"]
    subprocess.call(cmd)


def test_standalone_script(tmpdir):
    input_dir = os.sep.join((test_dir, 'resources'))
    sys.argv = ["test", "--input-directory", input_dir, "--working-directory", str(tmpdir), "--force"]
    m.main()


def test_version():
    cmd = ["sequana_lora", "--version"]
    subprocess.call(cmd)
