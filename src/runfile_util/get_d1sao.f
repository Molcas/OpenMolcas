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
      Subroutine Get_D1sao(D1sao,nD1sao)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"
      Character*24 Label
#ifdef _DEBUGPRINT_
#include "run_common.fh"
#endif
      Logical      Found
      Real*8 D1sao(nD1sao)

      Call Get_iScalar('System BitSwitch',iOption)
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*                                                                      *
************************************************************************
*                                                                      *
      Label='D1sao'
      Call qpg_dArray(Label,Found,nDens)
      If(.not.Found .or. nDens.eq.0) Then
         Call SysAbendMsg('get_d1sao','Did not find',Label)
      End If
      If (nDens/=nD1sao) Then
         Write (6,*) 'Get_D1sao: nDens/=nD1sao'
         Call Abend()
      End If
      Call Get_dArray(Label,D1sao,nD1sao)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      if(is_nSym.eq.0) then
       Call Get_iScalar('nSym',nSym)
       is_nSym=1
      endif
      if(is_nBas.eq.0) then
       Call Get_iArray('nBas',nBas,nSym)
       is_nBas=1
      endif
      Write(6,*) 'variational 1st order density matrix'
      ii=1
      Do iIrrep = 0, nSym - 1
         If (nBas(iIrrep).gt.0) Then
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',D1sao(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End If
      End Do
#endif
*
      Return
      End
