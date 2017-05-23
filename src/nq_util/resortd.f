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
      Subroutine ResortD(D_Old,D_New,iBas,iCmp,jBas,jCmp)
      Implicit Real*8 (a-h,o-z)
*
      Real*8 D_Old(iBas,jBas,iCmp,jCmp), D_New(iBas,iCmp,jBas,jCmp)
*
      Do jC = 1, jCmp
         Do jB = 1, jBas
            Do iC = 1, iCmp
               Do iB = 1, iBas
                  D_New(iB,iC,jB,jC)=D_Old(iB,jB,iC,jC)
               End Do
            End Do
         End Do
      End Do
*
      Return
      End
