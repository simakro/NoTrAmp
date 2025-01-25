from setuptools import setup


__version__ = ""
exec(open('notramp/version.py').read())

with open("README.md", "r", encoding="utf-8") as desc:
    long_description = desc.read()

setup(
    name="notramp",
    version=__version__,
    description="Super-fast Normalization and Trimming for Amplicon Sequencing Data (Long- and Short-read)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/simakro/NoTrAmp.git",
    author="Simon Magin",
    author_email="simon.magin@uk-essen.de",
    license="BSD-2",
    packages=[
        "notramp",
        "notramp.resources",
        "notramp.test",
        "notramp.test_self",
    ],
    include_package_data=True,
    entry_points={"console_scripts": ["notramp = notramp.notramp_main:run_notramp"]},
    install_requires=[
        "psutil",
    ],
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
