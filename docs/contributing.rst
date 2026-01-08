Contributing
============

We welcome contributions to genbank_to! This document provides guidelines for contributing to the project.

Getting Started
---------------

1. **Fork the Repository**

   Fork the repository on GitHub: https://github.com/linsalrob/genbank_to

2. **Clone Your Fork**

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/genbank_to.git
      cd genbank_to

3. **Set Up Development Environment**

   .. code-block:: bash

      # Create a virtual environment
      python -m venv venv
      source venv/bin/activate
      
      # Install in editable mode with test dependencies
      pip install -e ".[test]"

4. **Create a Branch**

   .. code-block:: bash

      git checkout -b feature/your-feature-name

How to Contribute
-----------------

Reporting Bugs
~~~~~~~~~~~~~~

If you find a bug, please create an issue on GitHub with:

- A clear, descriptive title
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Your environment (OS, Python version, genbank_to version)
- Sample data if possible (or a link to it)

Use our `bug report template <https://github.com/linsalrob/genbank_to/issues/new>`_.

Suggesting Enhancements
~~~~~~~~~~~~~~~~~~~~~~~

We welcome feature requests! Please create an issue with:

- A clear description of the feature
- Why you think it would be useful
- Example use cases
- Any relevant examples from other tools

Use our `feature request template <https://github.com/linsalrob/genbank_to/issues/new>`_.

Adding Output Formats
~~~~~~~~~~~~~~~~~~~~~

We're always interested in supporting new output formats! If you'd like to add a new format:

1. Create an issue describing the format
2. Provide links to format specifications
3. Include example input and expected output
4. If you've already implemented it, submit a pull request!

Code Contributions
------------------

Making Changes
~~~~~~~~~~~~~~

1. **Write Code**

   - Follow the existing code style (see Style Guidelines below)
   - Keep changes focused and minimal
   - Add docstrings to new functions
   - Update type hints where applicable

2. **Add Tests**

   Add tests for your changes in the ``test/`` directory:

   .. code-block:: python

      def test_new_feature():
          """Test description."""
          # Your test code
          assert expected == actual

3. **Run Tests**

   .. code-block:: bash

      # Run all tests
      pytest
      
      # Run specific test file
      pytest test/test_genbank_to.py
      
      # Run with coverage
      pytest --cov=GenBankToLib

4. **Update Documentation**

   If you've added a new feature or changed behavior:
   
   - Update relevant ``.rst`` files in ``docs/``
   - Add examples to ``docs/examples.rst``
   - Update API documentation in ``docs/api.rst``
   - Update the README.md if appropriate

5. **Commit Changes**

   Write clear, descriptive commit messages:

   .. code-block:: bash

      git add .
      git commit -m "Add support for XYZ format
      
      - Implement xyz_converter function
      - Add tests for XYZ format
      - Update documentation with examples"

6. **Push and Create Pull Request**

   .. code-block:: bash

      git push origin feature/your-feature-name

   Then create a pull request on GitHub with:
   
   - Clear description of changes
   - Link to related issues
   - Screenshots if applicable
   - Confirmation that tests pass

Style Guidelines
----------------

Python Code Style
~~~~~~~~~~~~~~~~~

- Follow PEP 8 for Python code
- Use 4 spaces for indentation
- Maximum line length: 120 characters
- Use descriptive variable names

.. code-block:: python

   # Good
   def genbank_to_format(genbank_file, output_format='fasta'):
       """Convert GenBank file to specified format."""
       pass
   
   # Avoid
   def gb2fmt(gbf, fmt='fasta'):
       pass

Docstring Format
~~~~~~~~~~~~~~~~

Use Google-style docstrings:

.. code-block:: python

   def genbank_to_format(gbkf, output_format='fasta', skip_pseudo=True):
       """
       Convert a GenBank file to the specified format.
       
       Args:
           gbkf (str): Path to the GenBank file.
           output_format (str, optional): Desired output format. 
               Default is 'fasta'.
           skip_pseudo (bool, optional): Skip pseudogenes. 
               Default is True.
       
       Returns:
           str: Path to the output file.
       
       Raises:
           FileNotFoundError: If the GenBank file doesn't exist.
           ValueError: If the output format is not supported.
       
       Example:
           >>> genbank_to_format('genome.gbk', 'gff3')
           'genome.gff3'
       """
       pass

Testing Guidelines
------------------

Writing Tests
~~~~~~~~~~~~~

- Write tests for all new functionality
- Use descriptive test names
- Test edge cases and error conditions
- Keep tests independent and isolated

.. code-block:: python

   import pytest
   from GenBankToLib import genbank_to_faa
   
   def test_genbank_to_faa_basic():
       """Test basic protein extraction."""
       results = list(genbank_to_faa('test/NC_001417.gbk'))
       assert len(results) > 0
       assert all(len(r) == 3 for r in results)
   
   def test_genbank_to_faa_skip_pseudo():
       """Test that pseudogenes are skipped by default."""
       results_with = list(genbank_to_faa('test/pseudo.gbk', skip_pseudo=False))
       results_without = list(genbank_to_faa('test/pseudo.gbk', skip_pseudo=True))
       assert len(results_with) > len(results_without)
   
   def test_genbank_to_faa_file_not_found():
       """Test error handling for missing file."""
       with pytest.raises(FileNotFoundError):
           list(genbank_to_faa('nonexistent.gbk'))

Test Data
~~~~~~~~~

- Place test files in the ``test/`` directory
- Use small, representative test files
- Document the source of test data
- Include expected output files when appropriate

Documentation Guidelines
------------------------

Documentation Structure
~~~~~~~~~~~~~~~~~~~~~~~

- Use reStructuredText (.rst) format
- Keep documentation up-to-date with code changes
- Include practical examples
- Link to relevant sections

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

Build the documentation locally to check your changes:

.. code-block:: bash

   cd docs
   pip install -r requirements-docs.txt
   make html
   
   # View the built documentation
   open _build/html/index.html

Documentation Style
~~~~~~~~~~~~~~~~~~~

- Use clear, concise language
- Provide code examples
- Explain both what and why
- Include expected output when helpful

Pull Request Process
--------------------

1. **Before Submitting**

   - Ensure all tests pass
   - Update documentation
   - Check code style
   - Rebase on latest main branch

2. **Pull Request Description**

   Include:
   
   - Summary of changes
   - Motivation and context
   - Type of change (bug fix, new feature, etc.)
   - How to test
   - Related issues

3. **Review Process**

   - Address reviewer feedback
   - Keep discussion constructive
   - Update PR as needed
   - Be patient - maintainers review in their spare time

4. **After Merge**

   - Delete your branch
   - Update your fork
   - Celebrate! 🎉

Code of Conduct
---------------

We are committed to providing a welcoming and inclusive environment. Please be respectful and professional in all interactions.

- Be kind and courteous
- Respect differing viewpoints
- Accept constructive criticism gracefully
- Focus on what's best for the community
- Show empathy towards others

Communication
-------------

- **GitHub Issues**: Bug reports and feature requests
- **Pull Requests**: Code contributions and discussions
- **Email**: raedwards@gmail.com for security issues or private concerns

Recognition
-----------

Contributors are recognized in:

- The repository's contributor list on GitHub
- The CITATION.cff file for academic citations
- Release notes for significant contributions

Thank You!
----------

Thank you for contributing to genbank_to! Your contributions help make bioinformatics tools more accessible and useful for the research community.

Questions?
----------

If you have questions about contributing, please:

1. Check existing issues and documentation
2. Create a new issue with your question
3. Reach out to the maintainers

We're here to help!
