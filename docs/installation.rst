Installation
============

genbank_to can be installed in several ways depending on your needs and preferences.

Using Conda (Recommended)
--------------------------

This is the easiest and recommended method for most users:

.. code-block:: bash

   mamba create -n genbank_to genbank_to
   conda activate genbank_to
   genbank_to --help

If you prefer using conda instead of mamba:

.. code-block:: bash

   conda create -n genbank_to genbank_to
   conda activate genbank_to
   genbank_to --help

Using pip
---------

You can install genbank_to from PyPI using pip. We recommend using a virtual environment:

.. code-block:: bash

   # Create a virtual environment
   python -m venv venv
   
   # Activate the virtual environment
   # On Linux/macOS:
   source venv/bin/activate
   # On Windows:
   # venv\Scripts\activate
   
   # Install genbank_to
   pip install genbank_to
   
   # Verify installation
   genbank_to --help

From Source (Development)
--------------------------

If you want to contribute to genbank_to or use the latest development version:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/linsalrob/genbank_to.git
   cd genbank_to
   
   # Create a virtual environment
   python -m venv venv
   source venv/bin/activate
   
   # Install in editable mode
   pip install -e .
   
   # Verify installation
   genbank_to --help

Dependencies
------------

genbank_to requires Python 3.9 or later and the following packages:

- **biopython** (>=1.74): For parsing GenBank files
- **numpy** (>=1.16.0): For numerical operations
- **pandas**: For data manipulation
- **bcbio-gff** (>=0.6.6): For GFF format support

These dependencies will be automatically installed when you install genbank_to using any of the methods above.

Verifying Installation
----------------------

After installation, verify that genbank_to is correctly installed:

.. code-block:: bash

   # Check version
   genbank_to --version
   
   # View help
   genbank_to --help

You should see the version information and command-line options.

Testing Your Installation
--------------------------

You can test your installation with a sample GenBank file:

.. code-block:: bash

   # Download test data
   wget https://raw.githubusercontent.com/linsalrob/genbank_to/main/test/NC_001417.gbk
   
   # Convert to FASTA
   genbank_to -g NC_001417.gbk -n NC_001417.fna
   
   # Check the output
   head NC_001417.fna

Troubleshooting
---------------

ImportError: No module named 'GenBankToLib'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error occurs when genbank_to is not properly installed. Make sure you've activated your virtual environment and installed the package:

.. code-block:: bash

   source venv/bin/activate  # or conda activate genbank_to
   pip install genbank_to

Command not found: genbank_to
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This happens when the installation directory is not in your PATH. If using a virtual environment, make sure it's activated. If using conda, activate the environment:

.. code-block:: bash

   conda activate genbank_to

BioPython version issues
~~~~~~~~~~~~~~~~~~~~~~~~

If you encounter BioPython-related errors, ensure you have a compatible version:

.. code-block:: bash

   pip install --upgrade biopython>=1.74

Upgrading
---------

To upgrade to the latest version:

Using pip:

.. code-block:: bash

   pip install --upgrade genbank_to

Using conda:

.. code-block:: bash

   conda update genbank_to

Uninstalling
------------

To uninstall genbank_to:

Using pip:

.. code-block:: bash

   pip uninstall genbank_to

Using conda:

.. code-block:: bash

   conda remove genbank_to
