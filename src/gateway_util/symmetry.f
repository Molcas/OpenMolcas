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
      integer function symmetry(isym,jsym)
CSVC: returns the product of two irreps (1-based indexing of irreps!)
      implicit none
      integer isym, jsym
      symmetry = IEOR(isym-1,jsym-1)+1
      return
      end
