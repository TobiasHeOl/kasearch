from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='kasearch',
    version='0.0.15',
    description='KA-Search: Rapid and exhaustive sequence identity search of known antibodies',
    license='BSD 3-clause license',
    maintainer='Tobias Hegelund Olsen; Brennan Abanades kenyon',
    long_description=long_description,
    long_description_content_type='text/markdown',
    include_package_data=True,
    packages=find_packages(include=('kasearch', 'kasearch.*')),
    package_data={
        '': ['*.txt']
    },
    install_requires=[
        'pandas',
        'joblib',
        'requests',
        'numpy',
        'jax',
        'jaxlib'
    ],
)
