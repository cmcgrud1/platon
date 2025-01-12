import setuptools
from platon import __version__, name

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = name,
    version = __version__,
    author = "Michael Zhang, Yayaati Chachan",
    author_email = "zmzhang@caltech.edu",
    description = "A package to compute transmission spectra and retrieve atmospheric parameters from transmission spectra",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/ideasrule/platon",
    packages = setuptools.find_packages(),
    classifiers = (
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
    include_package_data = True,
    zip_safe = False,
    install_requires = [
        "numpy",
        'scipy < 1.3;python_version < "3.5"',
        'scipy;python_version >= "3.5"',
        "emcee", "dynesty",
        "future", "nose", "setuptools", "configparser",
        'matplotlib < 3;python_version < "3.5"',
        'matplotlib;python_version >= "3.5"']
)
