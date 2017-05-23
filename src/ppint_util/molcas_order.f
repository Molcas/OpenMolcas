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
      Subroutine Molcas_Order(GOut,na,nb)
      Implicit Real*8 (a-h,o-z)
      Real*8 GOut(na*nb,2)
*                                                                      *
************************************************************************
*                                                                      *
C     Call RecPrt('Gout(Argos)',' ',Gout,nb,na)
      call dcopy_(na*nb,GOut(1,1),1,GOut(1,2),1)
*
      Do ia = 1, na
         Do ib = 1, nb
            iba=(ia-1)*nb+ib
            iab=(ib-1)*na+ia
            GOut(iab,1)=GOut(iba,2)
         End Do
      End Do
C     Call RecPrt('Gout(Molcas)',' ',Gout,na,nb)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
