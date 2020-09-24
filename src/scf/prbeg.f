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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine PrBeg(Meth)
************************************************************************
*                                                                      *
*     purpose: Print the header to iterations                          *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*

#include "mxdm.fh"
#include "infscf.fh"
      Character*(*) Meth
      Character*10 Label
      Character *4 cUHF
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUGPRINT_
      Call qEnter('PrBeg')
#endif
*
      If (jPrint.ge.2) Then
*
      Write(6,*)
      call CollapseOutput(1,'Convergence information')
      iDummy_run=0
      cUHF='    '
      if(iUHF.eq.1) cUHF='UHF '
      Label=Meth(1:10)
      If (nIter(nIterP).gt.0) Then
        Write(6,'(31x,A,A,A)') cUHF,Label,' iterations: Energy and '//
     &                     'convergence statistics'
        Write(6,*)
        Write(6,'(A,A,A)')'Iter     Tot. ',Label,' One-elec.   '//
     &          '    Two-elec.     Energy      Max Dij or  Max Fij '//
     &          '     DNorm      TNorm      AccCon     Time'
        Write(6,'(A)')'         Energy          Energy          '//
     &          'Energy        Change      Delta Norm           '//
     &          '                                     in Sec.'
      Else
         iDummy_run=1
         Write(6,'(45x,A)')'No optimization is performed'
         If (InVec.eq.1) Then
            Write(6,'(29x,A)') 'Results refer to orbitals obtained '
     &                       //'from core diagonalization'
         Else If (InVec.eq.2) Then
            Write(6,'(34x,A,A)') 'Results refer to input orbitals '
     &            //'read from ',Trim(SCF_FileOrb)
         Else If (InVec.eq.3) Then
            Write(6,'(34x,A)') 'Results refer to density matrix '
     &                       //'read from COMOLD'
         End If
      End If
*
      End If
*
#ifdef _DEBUGPRINT_
      Call qExit('PrBeg')
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
