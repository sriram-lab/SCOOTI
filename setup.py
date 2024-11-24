from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'SCOOTI'
LONG_DESCRIPTION = 'SCOOTI: Single Cell Optimization OBjective and Tradeoff Inference'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="SCOOTI", 
        version=VERSION,
        author="Chandrasekaran Lab, University of Michigan",
        author_email="csriram@umich.edu",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'SCOOTI'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Research",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: Linux :: Ubuntu",
        ]
)
