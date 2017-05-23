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
      Subroutine vClr(a,lda,nItem)
      Implicit Real*8 (a-h,o-z)
      Dimension a(*)
      Inda=1
      Do 100 i=1,nItem
         A(Inda)=0.0D0
         Inda=Inda+lda
100   Continue
      Return
      End
