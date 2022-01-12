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

Functionals are defined in driver.f90 through their uniqe func_id in LibXC and
their associated coefficients are set there too.
The same routine contains logic for the combinations of functionals that are
used in the orbital-free embedding implementation by Aquilante et al.

Note that driver.f90 should not be clone ever. If new functionalities are to be
added these should be incoporated in this routine.

In the case of functionals not present in LibXC the corresponding subroutines
are places in src/raw_functionals. Driver.f90 do support the use of external
function in combinations with or not functionals of the LibXC library.
