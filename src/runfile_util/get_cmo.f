************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Roland Lindh                                           *
************************************************************************
*  Get_Cmo
*
*> @brief
*>   Return the symmetry blocked MO coefficients as read from the runfile
*> @author R. Lindh
*>
*> @details
*> The utility will read the symmetry blocked MO coefficients from the runfile.
*>
*> @param[out] CMO array of symmetry blocked MO coefficients
*> @param[out] nCMO  Number of elements in the array of symmetry blocked MO coefficients
************************************************************************
      Subroutine Get_Cmo(CMO,nCMO)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      Character Label*24
      Integer is_nSym, nSym
      Integer is_nBas, nBas(0:7)
      save is_nSym, is_nBas
      data is_nSym/0/, is_nBas/0/
      save nSym, nBas
      Logical   IfTest, Found
      Data      IfTest/.False./
      Real*8 CMO(nCMO)
      Call Get_iScalar('System BitSwitch',iOption)
#ifdef _DEBUG_
      IfTest=.True.
#endif

      Label='Last orbitals'
      Call qpg_dArray(Label,Found,mCmo)
      If(.not.Found) Then
         Call SysAbendMsg('get_CMO','Could not find',Label)
      End If
      If (mCMO/=nCMO) Then
         Write (6,*) 'Get_CMO: mCMO/=nCMO'
         Write (6,*) 'nCMO=',nCMO
         Write (6,*) 'mCMO=',mCMO
         Call Abend()
      End If
      Call Get_dArray(Label,CMO,nCMO)
*                                                                      *
************************************************************************
*                                                                      *
      If (IfTest) Then
         if(is_nSym.eq.0) then
          Call Get_iScalar('nSym',nSym)
          is_nSym=1
         endif
         if(is_nBas.eq.0) then
          Call Get_iArray('nBas',nBas,nSym)
          is_nBas=1
         endif
         Write (6,*) ' Input Orbitals from RUNFILE'
         Write (6,*)
         ii=1
         Do iIrrep = 0, nSym-1
            If (nBas(iIrrep).gt.0) Then
               Write (6,*) ' Symmetry Block',iIrrep
               Call RecPrt(' ',' ',CMO(ii),nBas(iIrrep),nBas(iIrrep))
               Write (6,*)
            End If
            ii = ii + nBas(iIrrep)**2
         End Do
      End If
*
      Return
      End
