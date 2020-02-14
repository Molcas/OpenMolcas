OpenMolcas
==========

OpenMolcas is a quantum chemistry software package developed by scientists
and intended to be used by scientists. It includes programs to apply many
different electronic structure methods to chemical systems, but its key
feature is the multiconfigurational approach, with methods like CASSCF and
CASPT2.

OpenMolcas is not a fork or reimplementation of
[Molcas](http://www.molcas.org), it *is* a large part of the Molcas codebase
that has been released as free and open-source software (FOSS) under the Lesser
General Public License (LGPL). Some parts of Molcas remain under a different
license by decision of their authors (or impossibility to reach them), and are
therefore not included in OpenMolcas.

**Latest reference**:
"OpenMolcas: From Source Code to Insight."
*J. Chem. Theory Comput.* **15** (2019) 5925-5964.
DOI: [10.1021/acs.jctc.9b00532](https://doi.org/10.1021/acs.jctc.9b00532)

Installation
------------

For more detailed information, please refer to the [wiki
pages](/../wikis/home).

OpenMolcas is configured with [CMake](https://cmake.org). A quick way to get it
up and running is the following:

1.  Clone the repository:

    ```
    git clone https://gitlab.com/Molcas/OpenMolcas.git
    ```

2.  Get the `lapack` submodule (only needed if you don't use another linear
    algebra library like MKL or OpenBLAS):

    ```
    cd OpenMolcas
    git submodule update --init External/lapack
    cd ..
    ```

3.  Create a new directory and run `cmake` from it:

    ```
    mkdir build
    cd build
    cmake ../OpenMolcas
    ```

4.  Compile with `make`:

    ```
    make
    ```

5.  Run the verification suite (failures in "grayzone" tests are expected):

    ```
    pymolcas verify
    ```

For running other calculations you should define the `MOLCAS` environment
variable to point to the `build` directory. Run `pymolcas --help` to see the
available options of the script. In particular it is recommended to run:
```
pymolcas -setup
```
for your first installation.

Documentation
-------------

The documentation can be found in the
[`doc`](https://gitlab.com/Molcas/OpenMolcas/tree/master/doc) directory, you
can read it in [HTML format](https://molcas.gitlab.io/OpenMolcas/sphinx/) or
[PDF format](https://molcas.gitlab.io/OpenMolcas/Manual.pdf). Note that most
of it precedes the creation of OpenMolcas and it is probably outdated in
several points. It may also mention features not available in OpenMolcas.

Help
----

OpenMolcas is a community-supported software and as such it doesn't have an
official technical support. If you have any problems or questions, you can use
the [Issues](/../issues) page or the [Molcas
forum](https://cobalt.itc.univie.ac.at/molcasforum/index.php), and hopefully
some other user or developer will be able to help you.

If you need technical support, you can acquire a [Molcas
license](http://www.molcas.org/order.html).

Contributing
------------

Since OpenMolcas is FOSS, you can download it, modify it and distribute it
freely (according to the terms of the LGPL). If you would like your
contributions to be included in the main repository, please contact one of the
developers, write a message in the
[forum](https://cobalt.itc.univie.ac.at/molcasforum/index.php) or submit a
[merge request](https://docs.gitlab.com/ee/gitlab-basics/README.html). We are
currently only accepting [Molcas
developers](http://www.molcas.org/cgi-bin/dev.plx) as members of the `Molcas`
group, but everyone is welcome to send patches, suggestions and bug reports. 
