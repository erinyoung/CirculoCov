from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='circulocov',
    version='0.1.20230104',
    author='Erin Young',
    author_email='eriny@utah.gov',
    url="https://github.com/erinyoung/CirculoCov",
    description='Circular-Aware Coverage for Draft Genomes',
    long_description=long_description,
    long_description_content_type="text/markdown", 
    packages=find_packages(include=['circulocov', 'circulocov.*']),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta ',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GPL-3.0',
        'Topic :: Scientific/Engineering :: Bio-Informatics '
        ],
    keywords = [
        "bioinformatics",
        "coverage",
        "visualization",
    ],
    python_requires=">=3.8, <4",
    install_requires=[
        "pandas",
        "pysam",
        "matplotlib",
        "numpy",
    ],
    entry_points={
        'console_scripts': [
            'circulocov=circulocov.circulocov:main'
            ],
    },
)