import asyncio
import json
import shutil
from pathlib import Path
from typing import Dict, Tuple

from sequana_pipelines import lora

LORA_PATH = Path(lora.__file__).parent


def base_version_parser(stdout: str):
    return stdout.split("\n", 1)[0].split()[-1]


def strip_version_parser(stdout: str):
    return stdout.strip()


def seqtk_version_parser(stderr: str):
    return stderr.strip().split("\n")[1].split()[-1]


def unicycler_version_parser(stdout: str):
    return stdout.split("\n", 1)[0].split()[1]


def ccs_version_parser(stdout: str):
    return stdout.split("\n", 1)[0].split()[1]


def checkm_version_parser(stdout: str):
    return stdout.split("\n", 1)[1].split()[2]


# dict to select the good parser with the first key indicate the error code
VERSION_PARSER = {
    0: {
        "ccs": ccs_version_parser,
        "checkm": checkm_version_parser,
        "circlator": strip_version_parser,
        "flye": strip_version_parser,
        "unicycler": unicycler_version_parser,
        "hifiasm": strip_version_parser,
        "minimap2": strip_version_parser,
    },
    1: {
        "seqtk": seqtk_version_parser,
    },
}
APPTAINER_CORRESPONDANCE = {
    "canu_correction": "canu",
    "sequana_coverage": "sequana",
}


def get_tools_versions(config):
    """Run used tools to get version for dependencies."""
    # get used tools
    with open(LORA_PATH / "requirements.json") as fhin:
        tools = json.load(fhin)
    assembler = config["assembler"]
    used_tools = {
        assembler: tools["assembler"][assembler],
        **tools["utils"],
        **dict(_iter_optional_tools(tools["optional"], config)),
    }
    # just to avoid issue with deprecated lora
    apptainers = config.get("apptainers", {})
    return asyncio.run(run_tools_versions(used_tools, apptainers))


def _iter_optional_tools(tools, config):
    for tool, subtools in tools.items():
        if config[tool]["do"]:
            try:
                for subtool, cmd in subtools.items():
                    yield subtool, cmd
            except AttributeError:
                # whitout subtools, subtools is the version cmd
                yield tool, subtools


def get_apptainer_version(tool: str, apptainers: Dict[str, str]) -> str:
    try:
        apptainer = apptainers[tool]
    except KeyError:
        # some rule name differ from apptainer name
        apptainer = apptainers.get(APPTAINER_CORRESPONDANCE.get(tool, ""), "No version found")
    # sequana container finish with "_version.img"
    return apptainer.split("_")[-1].strip(".img")


async def run_version(tool: str, cmd: str, has_apptainer: bool, apptainers: Dict[str, str]) -> Tuple[str, str]:
    proc = await asyncio.create_subprocess_shell(cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate()

    # Code error 127: command not found
    if proc.returncode is None or proc.returncode == 127:
        # user may use apptainer
        if has_apptainer:
            return tool, get_apptainer_version(tool, apptainers)
        return tool, "Tool not found"

    # Some software return error code 1 with version information
    if proc.returncode != 0 and tool not in VERSION_PARSER[proc.returncode]:
        return tool, "No version found"

    # get the good parser for the tool
    version_parser = VERSION_PARSER.get(proc.returncode, {}).get(tool, base_version_parser)
    output = stdout.decode() if proc.returncode == 0 else stderr.decode()
    try:
        return tool, version_parser(output)
    except IndexError:
        return tool, "No version found"


async def run_tools_versions(used_tools: Dict["str", "str"], apptainers: Dict["str", "str"]):
    # check if user has apptainer/singularity
    has_apptainer = (shutil.which("apptainer") or shutil.which("singularity")) is not None
    result = await asyncio.gather(
        *[run_version(tool, cmd, has_apptainer, apptainers) for tool, cmd in used_tools.items()]
    )
    return result
