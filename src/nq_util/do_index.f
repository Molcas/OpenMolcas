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
      Subroutine Do_Index(Index,NrBas,NrBas_Eff,iCmp)
      Implicit Real*8 (a-h,o-z)
      Integer Index(NrBas_Eff*iCmp)
*
      iAdd=NrBas-NrBas_Eff
      Do iB_Eff = 1, NrBas_Eff
         iB = iB_Eff + iAdd
         Do iC = 1, iCmp
            iCB = (iC-1)*NrBas + iB
            iCB_Eff = (iC-1)*NrBas_Eff + iB_Eff
            Index(iCB_Eff)=iCB
         End Do
      End Do
*
      Return
      End
