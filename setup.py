"""
setup.py — nustar_uplim
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="nustar_uplim",
    version="1.0.0",
    author="nustar_uplim contributors",
    description="NuSTAR non-detection upper limit analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-org/nustar_uplim",
    packages=find_packages(exclude=["tests*", "docs*"]),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21",
        "scipy>=1.7",
        "astropy>=5.0",
        "matplotlib>=3.5",
        "regions>=0.6",
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov", "flake8"],
    },
    entry_points={
        "console_scripts": [
            "nustar-uplim=nustar_uplim.pipeline:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="nustar x-ray astronomy upper-limits non-detection",
)
