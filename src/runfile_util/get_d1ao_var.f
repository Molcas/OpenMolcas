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
      Subroutine Get_D1ao_Var(ipD1ao,nDens)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"

      Character*24 Label
#ifdef _DEBUG_
#include "run_common.fh"
#endif
      Logical      Found

      Call Get_iScalar('System BitSwitch',iOption)
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*                                                                      *
************************************************************************
*                                                                      *
      Label='D1aoVar'
      Call qpg_dArray(Label,Found,nDens)
      If(.not.Found .or. nDens.eq.0) Then
         Call Get_D1ao(ipD1ao,nDens)
      Else
         Call GetMem('D1ao_var','Allo','Real',ipD1ao,nDens)
         Call get_dArray(Label,Work(ipD1ao),nDens)
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      if(is_nSym.eq.0) then
       Call Get_iScalar('nSym',nSym)
       is_nSym=1
      endif
      if(is_nBas.eq.0) then
       Call Get_iArray('nBas',nBas,nSym)
       is_nBas=1
      endif
      Write(6,*) 'variational 1st order density matrix'
      ii=ipD1ao
      Do iIrrep = 0, nSym - 1
         If (nBas(iIrrep).gt.0) Then
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',Work(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End If
      End Do
#endif
*
      Return
      End
