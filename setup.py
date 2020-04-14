import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = {}
with open("pygmes/version.py") as fp:
    exec(fp.read(), version)
setuptools.setup(
    name="pygmes",
    version=version["__version__"],
    author="Paul Saary",
    author_email="saary@ebi.ac.uk",
    description="Run GeneMark-ES using pretrained models",
    url="https://github.com/openpaul/pygmes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    py_modules=["api"],
    entry_points={"console_scripts": ["pygmes = pygmes.api:main"]},
    install_requires=["ete3", "pyfaidx>=0.5.8"],
    packages=setuptools.find_packages(),
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
    ],
)
