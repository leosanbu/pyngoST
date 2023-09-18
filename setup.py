from setuptools import setup

setup(
    name='pyngoST',
    version='1.0.0',
    author='Leonor Sanchez-Buso',
    author_email='leonor.sanchez@fisabio.es',
    description='pyngoST: multiple sequence typing of Neisseria gonorrhoeae large assembly collections',
    url='https://github.com/leosanbu/pyngoST',
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
