from setuptools import setup, find_packages

setup(
    name="pattools-methy",
    packages=find_packages(),
    version="0.1.0",
    author="Department of research and development, Zhejiang Gaomei Genomics",
    author_email="it@gomicsgene.com",
    description="pattools is a BS-seq analysis tool suite based on pat format",
    url='https://hcyvan.github.io/pattools/index.html',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
    ],
    package_data={
        'pattools.deconv.moss': ['*.csv'],
        'pattools.deconv.sun': ['*.csv'],
        'pattools.deconv.loyfer': ['*.tsv']
    },
    install_requires=['pysam', 'cvxpy', 'scipy', 'pandas', 'matplotlib', 'scikit-learn', 'mpi4py', 'hdbscan'],
    entry_points={
        'console_scripts': [
            'pattools=pattools:main',
        ],
    },
    py_modules=[],
    python_requires='>=3.10, <4',
)
