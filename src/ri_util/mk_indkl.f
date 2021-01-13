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
      Subroutine Mk_Indkl(Indkl_OnOff,Indkl,nkl)
      Implicit Real*8 (a-h,o-z)
      Integer Indkl_OnOff(nkl), Indkl(nkl)
*
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call iVcPrt('Mk_Indkl: Indkl_OnOff',' ',Indkl_OnOff,nkl)
#endif
      ikl = 0
      Do jkl = 1, nkl
         If (Indkl_OnOff(jkl).eq.1) Then
            ikl = ikl + 1
            Indkl(jkl)=ikl
         Else
            Indkl(jkl)=0
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call iVcPrt('Mk_Indkl: Indkl',' ',Indkl,nkl)
#endif
*
      Return
      End
