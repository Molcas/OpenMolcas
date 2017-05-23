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
      Subroutine Get_Cmo(ipCMO,nCMO)
************************************************************
*
*   <DOC>
*     <Name>Get\_Cmo</Name>
*     <Syntax>Call Get\_Cmo(ipCMO,nCMO)</Syntax>
*     <Arguments>
*       \Argument{ipCMO}{pointer to array of symmetry blocked MO coefficients.}{Integer}{out}
*       \Argument{nCMO}{Number of elements in the array of symmetry blocked MO coefficients.}{Integer}{out}
*     </Arguments>
*     <Purpose>To return a pointer in Work of the location of the symmetry blocked MO coefficients as read from the run file.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>The utility will read the symmetry blocked MO coefficients from the run file.
*     </Description>
*    </DOC>
*
************************************************************
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
      Call Get_iScalar('System BitSwitch',iOption)
#ifdef _DEBUG_
      IfTest=.True.
#endif

      Label='Last orbitals'
      Call qpg_dArray(Label,Found,nCmo)
      If(.not.Found) Then
         Call SysAbendMsg('get_CMO','Could not find',Label)
      End If
      Call GetMem('CMO','Allo','Real',ipCMO,nCMO)
      Call Get_dArray(Label,Work(ipCMO),nCMO)
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
         ii=ipCMO
         Do iIrrep = 0, nSym-1
            If (nBas(iIrrep).gt.0) Then
               Write (6,*) ' Symmetry Block',iIrrep
               Call RecPrt(' ',' ',Work(ii),nBas(iIrrep),nBas(iIrrep))
               Write (6,*)
            End If
            ii = ii + nBas(iIrrep)**2
         End Do
      End If
*
      Return
      End
