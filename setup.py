from setuptools import setup, find_packages

setup(
    name="Metagenomics_pipeline4_V3",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "plotly", "kaleido",'distinctipy','numpy',
    ],
    entry_points={
    "console_scripts": [
        
"run_metagenomics_pl2=Metagenomics_pipeline4_V2.scripts.run_metagenomics_pl2:main",
    ],
},

    author="Harouna",
    author_email="harounasoum17@gmail.com",
    description="A bioinformatics pipeline for trimming, host depletion, and taxonomic classification",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Harounas/Metagenomics_pipeline4_V3.git",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
