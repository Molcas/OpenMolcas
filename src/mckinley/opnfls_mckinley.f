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
      SubRoutine OpnFls_McKinley()
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep, lIrrep
      Implicit Real*8(a-h,o-z)
#include "Molcas.fh"
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "disp2.fh"
#include "etwas.fh"
      Character*8 Method, MckLbl
      Character*288 Header
      iOpt = 1
      iRC = -1
      MckLbl='Title'
      Call cWrMck(iRC,iOpt,MckLbl,1,Header,iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'OpnFls: Error writing to MCKINT'
         Write (6,'(A,A)') 'MckLbl=',MckLbl
         Call QTrace
         Call Abend()
      End If
      iOpt = 1
      iRC = -1
      MckLbl='nSym'
      Call iWrMck(iRC,iOpt,MckLbl,1,[nIrrep],iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'OpnFls: Error writing to MCKINT'
         Write (6,'(A,A)') 'MckLbl=',MckLbl
         Call QTrace
         Call Abend()
      End If
      iOpt = 0
      iRC = -1
      MckLbl='nBas'
      Call iWrMck(iRC,iOpt,MckLbl,1,nBas,iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'OpnFls: Error writing to MCKINT'
         Write (6,'(A,A)') 'MckLbl=',MckLbl
         Call QTrace
         Call Abend()
      End If
      iOpt = 0
      iRC = -1
      MckLbl='SymOp'
      Call cWrMck(iRC,iOpt,MckLbl,1,lirrep(0),iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'OpnFls: Error writing to MCKINT'
         Write (6,'(A,A)') 'MckLbl=',MckLbl
         Call QTrace
         Call Abend()
      End If
      iOpt = 0
      iRC = -1
      MckLbl='ldisp'
      Call iWrMck(iRC,iOpt,MckLbl,1,ldisp,iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'OpnFls: Error writing to MCKINT'
         Write (6,'(A,A)') 'MckLbl=',MckLbl
         Call QTrace
         Call Abend()
      End If
      ngrad=0
      Do i=0,nIrrep-1
       nGrad=nGrad+ldisp(i)
      End Do
      iOpt = 0
      iRC = -1
      MckLbl='chdisp'
      Call cWrMck(iRC,iOpt,MckLbl,1,chdisp(1),iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'OpnFls: Error writing to MCKINT'
         Write (6,'(A,A)') 'MckLbl=',MckLbl
         Call QTrace
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*...  Get the method label
*
      Call Get_cArray('Relax Method',Method,8)
      If ( Method.eq.'RHF-SCF ' ) then
         nMethod=SCF
      Else if ( Method.eq.'CASSCF  ' ) then
         nMethod=RASSCF
      Else if ( Method.eq.'CASSCFSA' ) then
         nMethod=RASSCF
         Call Get_iScalar('SA ready',iGo)
         If (lHss.and.iGo.ne.2) Then
            Write (6,*)
            Write (6,*) ' Wavefunction type: RASSCF-SA'
            Write (6,*)
            Write (6,*) ' This option is not allowed when computing'//
     &                  ' the Hessian. Use the RHS option!'
            Call Quit_OnUserError()
         End If
      Else
         Write (6,*) ' OpnFls: Wavefunction type:',Method
         Write (6,*) '         Illegal type of wave function!'
         Write (6,*) '         McKinley can not continue'
         Write (6,*)
         Call Quit_OnUserError()
      End If
*
      Return
      End
