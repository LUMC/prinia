from setuptools import setup

setup(
    name='prinia',
    version='0.1',
    packages=['prinia'],
    install_required=["pyfaidx", "fastools", "lxml"],
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
