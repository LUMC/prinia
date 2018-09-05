from setuptools import setup, find_packages

setup(
    name='prinia',
    version='0.1',
    packages=find_packages(),
    install_requires=["biopython", "pysam", "pyfaidx",
                      "lxml", "pyvcf"],
    url='',
    license='',
    author='sander bollen',
    author_email='a.h.b.bollen@lumc.nl',
    description='primer design',
    entry_points={
        "console_scripts": [
            "primerdesign = prinia.primerdesign:main"
        ]
    }
)
