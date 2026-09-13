from setuptools import setup, find_packages

setup(
    name="theraflow",
    version="1.0.0",
    description="Multi-Tier In Silico Evaluation & Discovery Framework for Transcriptome-Driven Generative Chemistry",
    author="SydneyBioX",
    packages=find_packages(),
    py_modules=["run_pipeline"],
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.20.0",
        "scipy>=1.9.0",
        "scikit-learn>=1.0.0",
        "joblib>=1.2.0",
        "tqdm>=4.64.0",
        "chembl-webresource-client>=0.10.8",
        "posebusters>=0.2.0",
        "rdkit",
        "selfies>=2.1.0",
        "requests>=2.28.0",
    ],
    entry_points={
        "console_scripts": [
            "theraflow=run_pipeline:main",
        ],
    },
    python_requires=">=3.8",
)
