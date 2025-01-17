from setuptools import setup, find_packages

setup(
    name="PDBminer",
    version="1.0.0",
    description="Tool for mining PDB data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Kristine Degn",
    author_email="krde@dtu.com",
    url="https://github.com/ELELAB/PDBminer",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "biopython==1.83",
        "matplotlib==3.8.3",
        "networkx==3.2.1",
        "numpy==1.26.4",
        "pandas==2.2.1",
        "PyYAML==6.0.1",
        "requests==2.31.0",
        "seaborn==0.13.2",
    ],
    entry_points={
        "console_scripts": [
            "PDBminer=PDBminer.PDBminer:main",
            "PDBminer2coverage=PDBminer.PDBminer2coverage:main",
            "PDBminer2network=PDBminer.PDBminer2network:main",
        ],
    },
)
