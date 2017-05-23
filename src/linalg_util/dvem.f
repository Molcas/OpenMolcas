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
      Subroutine DVEM(n,rin1,inc1,rin2,inc2,rout,inc3)
      Real*8 rin1(inc1,*),rin2(inc2,*),rout(inc3,*)
      Integer n,i
      Do i=1,n
       rout(1,i)=rin1(1,i)*rin2(1,i)
      End Do
      Return
      End
