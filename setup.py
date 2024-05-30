from setuptools import setup,find_packages

setup(
    name="pattools",
    packages=find_packages(),
    version="0.0.10",
    author="gomics",
    author_email="it@gomicsgene.com",
    description="pattools is a BS-seq analysis tool suite based on pat format",
    classifiers=[
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        "Operating System :: OS Independent"
    ],
    package_data={
        'pattools.deconv.moss': ['*.csv'],
        'pattools.deconv.sun': ['*.csv'],
        'pattools.deconv.loyfer': ['*.tsv']
    },
    install_requires=['pysam', 'cvxpy', 'scipy'],
    entry_points={
        'console_scripts': [
            'pattools=pattools:main',
        ],
    },
    py_modules=[],
    python_requires='>=3.8, <4',
)
