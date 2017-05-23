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
      Function ExtNuc (ipExt,natom)
      Implicit Real*8 (A-H,O-Z)
*
* Compute Z - ExtPot interactions
*
#include "espf.fh"
*
      Real*8 E,ExtNuc
*
      Call QEnter('extnuc')
      iPL = iPL_espf()
*
      Call Get_dArray('Effective nuclear Charge',Charge,nAtom)
*
      E = Zero
      Do INuc = 1, natom
         E = E + Charge(iNuc)*Work(ipExt+(INuc-1)*10)
      End Do
      If (E.ne.zero .and. iPL.ge.3) Then
         Write(6,*) ' '
         Write(6,1000) E
1000     Format(' Ext Pot/(QM nuclei and MM charges) energy =',
     &          F16.10,' hartrees')
      End If
      ExtNuc = E
*
      Call QExit('extnuc')
      Return
      End
