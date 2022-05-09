from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='kasearch',
    version='0.0.1',
    description='Search OAS for similar sequences',
    license='BSD 3-clause license',
    maintainer='Tobias Olsen; Brennan Abanades',
    long_description=long_description,
    long_description_content_type='text/markdown',
    include_package_data=True,
    packages=find_packages(include=('kasearch', 'kasearch.*')),
    install_requires=[
        'numpy',
        'numba',
    ],
)
