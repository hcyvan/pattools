from setuptools import setup

setup(
    name="pattools",
    packages=['pattools'],
    version="0.0.1",
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
    install_requires=['pysam', 'cvxpy', 'scipy'],
    entry_points={
        'console_scripts': [
            'pattools=pattools:main',
        ],
    },
    py_modules=[],
    python_requires='>=3.8, <4',
)
