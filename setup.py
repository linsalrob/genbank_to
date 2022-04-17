import setuptools
from distutils.core import Extension
import os

# Read the markdown files for the long description
def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        if os.path.exists(filename):
            with open(filename, "r", encoding=encoding) as f:
                buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

def get_version():
    with open("VERSION", 'r') as f:
        v = f.readline().strip()
        return v

def get_requirements():
    reqs = []
    with open('requirements.txt', 'r') as f:
        for l in f:
            reqs.append(l.strip())
    return reqs

def main():
    setuptools.setup(
        name="genbank_to",
        version=get_version(),
        description="Convert GenBank format files to a swath of other formats",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Rob Edwards",
        platforms='any',
        keywords="genbanke bioinformatics microbiology genome genomics",
        author_email="raedwards@gmail.com",
        url='https://github.com/linsalrob/genbank_to',
        license='The MIT License (MIT)',
        packages=setuptools.find_packages(),
        include_package_data=True,
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Unix',
            'Programming Language :: Python :: 3.9',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires = [
            'biopython >= 1.74',
            'numpy >= 1.16.0',
            'pandas',
            'bcbio-gff >= 0.6.6'
        ],
        entry_points={
            "console_scripts": ["genbank_to = GenBankToLib.main:run"]
        }

    )

if __name__ == "__main__":
    main()
