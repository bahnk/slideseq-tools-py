"""
Packaging for slideseq-tools.
"""

from setuptools import setup, find_packages

setup(
    name="slideseq_tools",
    version="0.0.0",
    description="Tools to process Slide-seq data.",
    author="Nourdine Bah",
    author_email="nourdinebah@gmail.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["click"],
    entry_points={
        "console_scripts": [
            "synthetic_coordinates = slideseq_tools.scripts.synthetic_coordinates:main",
            "synthetic_data = slideseq_tools.scripts.synthetic_data:main",
            "check_slideseq_samplesheet = slideseq_tools.scripts.check_slideseq_samplesheet:main"
        ]
    },
)
