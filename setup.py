from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="metamodel",
    version="0.0.1",
    author="Christopher C Powers",
    author_email="c-11060@uri.edu",
    description='''Constructs a draft metabolic model from a metatransctriptomic
    metagenomic dataset''',
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="",
    packages = find_packages(),
    entry_points = {'console_scripts': ['metamodel = metamodel.__main__:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires = ['argh',
        'psamm >= 1.0',
        'biopython >= 1.79',
        'wget >= 3.2',
        'docopt >= 0.6.2'
        ]
)
