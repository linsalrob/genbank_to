# genbank_to Documentation

This directory contains the Sphinx documentation for genbank_to.

## Building the Documentation Locally

### Prerequisites

Install the required packages:

```bash
pip install -r requirements.txt
```

### Build HTML Documentation

```bash
make html
```

The built documentation will be in `_build/html/`. Open `_build/html/index.html` in your browser to view it.

### Build Other Formats

```bash
# PDF (requires LaTeX)
make latexpdf

# Plain text
make text

# ePub
make epub
```

### Clean Build Files

```bash
make clean
```

## Documentation Structure

- `index.rst` - Main landing page
- `installation.rst` - Installation instructions
- `quickstart.rst` - Quick start guide
- `usage.rst` - Command-line usage reference
- `output_formats.rst` - Output format descriptions
- `examples.rst` - Comprehensive examples
- `api.rst` - Python library API reference
- `contributing.rst` - Contribution guidelines
- `changelog.rst` - Version history

## ReadTheDocs

The documentation is automatically built and published on ReadTheDocs at:
https://genbank-to.readthedocs.io

Configuration is in `.readthedocs.yaml` in the repository root.

## Contributing to Documentation

When contributing to the documentation:

1. Use reStructuredText (.rst) format
2. Test builds locally before committing
3. Keep language clear and concise
4. Include code examples where helpful
5. Update the table of contents in `index.rst` if adding new pages

## Style Guide

- Use descriptive section headers
- Include code blocks with syntax highlighting
- Provide both command-line and library examples
- Link to related sections using cross-references
- Use admonitions (note, warning, tip) appropriately
