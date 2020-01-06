This directory contains different pieces of documentation.

The main part is the OpenMolcas documentation in
[reStructuredText](http://docutils.sourceforge.net/rst.html) markup format.
HTML and PDF versions can be generated with
[Sphinx](http://www.sphinx-doc.org). Note that most of this precedes the
creation of OpenMolcas and it is probably outdated in several points. It may
also mention features not available in OpenMolcas.

The `.rst` files contain embedded blocks that are used to generate command-line
help (accessed with `pymolcas help_doc`) and an XML input description for
MolGUI. These blocks use the `.. xmldoc::` directive and follow a format
similar to that described
[here](https://gitlab.com/Molcas/OpenMolcas/wikis/Programming%20guide/Documentation).

To build the documentation, after configuring OpenMolcas with `cmake`, run
`make doc_html` or `make doc_pdf`.

In addition, some subdirectories have other content:

`doxygen`

Configuration and auxiliary files to generate source code documentation with
[Doxygen](http://www.doxygen.org/). The documentation itself, where it exists,
is included in the source files.

To build the doxygen documentation just run `doxygen` inside this directory.
