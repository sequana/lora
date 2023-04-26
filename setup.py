from setuptools import setup


# handle sequana git link
with open("requirements.txt") as fh:
    requirements = [req.rstrip() if not req.startswith("git+") else req.rstrip().split("egg=")[-1] for req in fh]

_MAJOR = 0
_MINOR = 2
_MICRO = 0
version = f"{_MAJOR}.{_MINOR}.{_MICRO}"
release = f"{_MAJOR}.{_MINOR}"

metainfo = {
    "authors": {"main": ("thomas cokelaer", "thomas.cokelaer@pasteur.fr")},
    "version": version,
    "license": "new BSD",
    "url": "https://github.com/sequana/",
    "description": "long read assembly pipeline",
    "platforms": ["Linux", "Unix", "MacOsX", "Windows"],
    "keywords": ["pacbio,assembly,snakemake,sequana,canu"],
    "classifiers": [
        'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
}

NAME = "lora"


setup(
    name="sequana_{}".format(NAME),
    version=version,
    maintainer=metainfo["authors"]["main"][0],
    maintainer_email=metainfo["authors"]["main"][1],
    author=metainfo["authors"]["main"][0],
    author_email=metainfo["authors"]["main"][1],
    long_description=open("README.rst").read(),
    keywords=metainfo["keywords"],
    description=metainfo["description"],
    license=metainfo["license"],
    platforms=metainfo["platforms"],
    url=metainfo["url"],
    classifiers=metainfo["classifiers"],
    # package installation
    packages=[
        "sequana_pipelines.lora",
        "sequana_pipelines.lora.presets",
        "sequana_pipelines.lora.rules",
        "sequana_pipelines.lora.src",
        "sequana_pipelines.lora.src.templates"
    ],
    install_requires=requirements,
    # This is recursive include of data files
    exclude_package_data={"": ["__pycache__"]},
    package_data={
        "": [
            "*.yaml",
            "*.json",
            "requirements.txt",
            "*png",
            "*yml",
            "*smk",
            "*html",
            "*js",
            "*css"
        ],
    },
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "sequana_lora=sequana_pipelines.lora.main:main",
        ]
    },
)
