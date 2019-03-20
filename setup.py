import os
import sys
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open("README.md") as f:
    readme = f.read()

classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2.7',
    ]

# Requirements.
req2and3 = 'scipy,pandas,nose,sympy,scikit-image,subprocess32'.split(',')
if sys.version_info[0] == 3:
    requirements = req2and3 + ['matplotlib']
else:
    requirements = req2and3 + ['matplotlib==2.2.4']

setup(
    name = "FAST",
    version = "1.0.0",
    description = "Fast Automated Spud Trekker",
    long_description = readme,
    packages = ["FAST"],
    scripts = ['bin/fast', 'bin/lima', 'bin/ppss', 'bin/stack2tifs'],
    package_data = {},
    install_requires = requirements,
    author = "Tural Aksel et. al.",
    author_email = "",
    url = "http://github.com/turalaksel/FAST",
    license='GPL',
    classifiers=classifiers,
)
