import subprocess
import sys

import sequana_pipelines.lora.main as main

from . import test_dir

from click.testing import CliRunner

def test_standalone_subprocess(tmpdir):



    input_dir = test_dir / "resources"
    cmd = ["test", "--input-directory", str(input_dir), "--working-directory", str(tmpdir), 
            "--force", "--assembler", "flye"]
    subprocess.call(cmd)


def test_standalone_script(tmpdir):
    input_dir = test_dir / "resources"
    runner = CliRunner()

    results = runner.invoke(main.main, ["--input-directory", str(input_dir), "--working-directory", 
        str(tmpdir), "--force", "--pacbio", "--assembler", "flye"])

    assert results.exit_code == 0


def test_standalone_script_nanopore(tmpdir):
    input_dir = test_dir / "resources"
    runner = CliRunner()

    args = [
        "--input-directory",
        str(input_dir),
        "--working-directory",
        str(tmpdir),
        "--force",
        "--nanopore",
        "--mode",
        "eukaryotes",
        "--assembler",
        "flye"
    ]
    results = runner.invoke(main.main, args)
    assert results.exit_code == 0


def test_version():
    cmd = ["sequana_lora", "--version"]
    subprocess.call(cmd)

def test_help():
    cmd = ["sequana_lora", "--help"]
    subprocess.call(cmd)
