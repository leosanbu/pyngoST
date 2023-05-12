from setuptools import setup

setup(
    name='pyngoST',
    version='1.0.0',
    author='Leonor Sanchez-Buso',
    author_email='sanchez_leobus@gva.es',
    description='pyngoST: multiple sequence typing of Neisseria gonorrhoeae large assembly collections',
    packages=['pyngoST'],
    install_requires=[
        'urllib3==1.26.6',
        'requests',
        'biopython',
        'pandas',
        'pyfaidx',
        'openpyxl',
        'pyahocorasick'
    ],
)
