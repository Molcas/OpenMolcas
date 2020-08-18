************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************

      INTEGER FUNCTION mn2K(iab,lSym)
      use pso_stuff
      Integer  iab, lSym

      lab = iOff_ij2K(lSym) + iab
      mn2K = ij2K(lab)

      End
