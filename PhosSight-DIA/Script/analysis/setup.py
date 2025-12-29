from setuptools import setup, find_packages

setup(
    name="phossight-analysis",
    version="0.1.0",
    packages=find_packages(),
    py_modules=["protein_species_mapper", "PhosSight_paper_style"],
    python_requires=">=3.9"
)