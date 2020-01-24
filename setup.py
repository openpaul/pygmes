import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pygmes",
    version="0.1",
    author="Paul Saary",
    author_email="saary@ebi.ac.uk",
    description="Run GeneMark-ES with some bells",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    py_modules=[],
    entry_points={"console_scripts": ["pygmes = modules.__main__:main"]},
    install_requires=["ete3", "pyfaidx"],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
