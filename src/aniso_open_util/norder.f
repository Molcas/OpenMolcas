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
      Integer function norder(icoord,int_code,lmax)

      Implicit None
      Integer lmax
      Integer icoord(lmax),int_code(lmax)
      Integer nb,isite

      nb=0
      Do isite=1,lmax
         nb=nb+icoord(isite)*int_code(isite)
      End Do

      norder=nb+1
      Return
      End
