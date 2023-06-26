from setuptools import setup, find_packages

setup(
    name='big_library_design',
    version='0.1.0',
    packages=find_packages(include=['big_library_design','big_library_design/*']),
)
