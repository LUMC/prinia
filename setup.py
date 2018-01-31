from setuptools import setup

setup(
    name='prinia',
    version='0.1',
    packages=['prinia'],
    install_requires=["biopython", "pysam", "pyfaidx",
                      "lxml", "pyvcf"],
    url='',
    license='',
    author='sander bollen',
    author_email='a.h.b.bollen@lumc.nl',
    description='primer design',
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "primerdesign = prinia.primerdesign:main"
        ]
    }
)
