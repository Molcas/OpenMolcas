************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2022, Roland Lindh                                     *
************************************************************************

Notes on the dft_util.

OpenMolcas is since January 2022 completely interfaced with LibXC.

Functionals are defined in data/functionals.txt (or .src) through their
unique name in LibXC and their associated coefficients.

The routine "driver" contains logic for the combinations of functionals that are
used in the orbital-free embedding implementation by Aquilante et al.

Note that "driver" should not be cloned ever. If new functionalities are to be
added these should be incoporated in this routine.

In the case of functionals not present in LibXC the corresponding subroutines
are placed in this directory. "driver" does support the use of external
functions in combinations with or without functionals of the LibXC library.
