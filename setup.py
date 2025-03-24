from setuptools import setup, find_packages

#VERSION = '0.0.1' 
DESCRIPTION = 'SCOOTI'
LONG_DESCRIPTION = 'SCOOTI: Single Cell Optimization OBjective and Tradeoff Inference'

# Setting up
setup(
    name="SCOOTI",
    use_scm_version={"write_to": "SCOOTI/_version.py"},# Enables setuptools_scm to derive the version from SCM metadata
    setup_requires=["setuptools_scm"],  # Ensures setuptools_scm is installed during setup
    packages=find_packages(),           # Automatically finds packages (directories with __init__.py)
    author="Chandrasekaran Lab, University of Michigan",
    author_email="csriram@umich.edu",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    classifiers= [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Research",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Linux :: Ubuntu",
    ],
    python_requires=">=3.6",
)



"""
SCOOTI/
├── SCOOTI/
│   ├── __init__.py
│   ├── every_function.py
│   └── main.py
├── tests/
│   └── test_every_function.py
├── setup.py
├── README.md
├── LICENSE
└── MANIFEST.in   # (if you need to include non-code files)
"""
