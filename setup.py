
from setuptools import setup, find_namespace_packages

setup(
    name='bifrost_reporter',
    version='0.1.0',
    packages=find_namespace_packages(include=['bifrost_reporter',
                                              'bifrost_reporter.*'
                                             ]),
    scripts=['bin/bifrost_reporter'],
    install_requires=['pandas',
                      'pymongo',
                      'envyaml',
                      'openpyxl',
		      'seaborn'],
    include_package_data=True,
    zip_safe=False,
    author='Simone Scrima',
    author_email='sscr@dksund.dk',
    description='Bifrost Reporter',
    url='https://github.com/ssi-dk/bifrost_reporter',
)
