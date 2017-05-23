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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine GenerateB(CMO,nBas,nOrb2Loc,ipLbl_AO,ipLbl,nComp,Debug)
C
C     Author: T.B. Pedersen
C
C     Purpose: generate the dipole matrices for Boys localisation, i.e.
C              transform from AO to MO basis.
C
      Implicit Real*8 (a-h,o-z)
      Real*8  CMO(*)
      Integer ipLbl_AO(nComp), ipLbl(nComp)
      Logical Debug
#include "WrkSpc.fh"

      If (nBas.lt.1 .or. nOrb2Loc.lt.1) Return

      lDbar = nBas*nOrb2Loc
      Call GetMem('Dbar','Allo','Real',ipDbar,lDbar)
      Do iComp = 1,nComp
         Call DGEMM_('N','N',nBas,nOrb2Loc,nBas,
     &              1.0d0,Work(ipLbl_AO(iComp)),nBas,CMO,nBas,
     &              0.0d0,Work(ipDbar),nBas)
         Call DGEMM_('T','N',nOrb2Loc,nOrb2Loc,nBas,
     &              1.0d0,CMO,nBas,Work(ipDbar),nBas,
     &              0.0d0,Work(ipLbl(iComp)),nOrb2Loc)
      End Do
      Call GetMem('Dbar','Free','Real',ipDbar,lDbar)

      If (Debug) Then
         Write(6,*)
         Write(6,*) 'In GenerateB'
         Write(6,*) '------------'
         Write(6,*) '[Assuming doubly occupied orbitals]'
         Do iComp = 1,nComp
            ip0 = ipLbl(iComp) - 1
            Cmp = 0.0d0
            Do iMO = 1,nOrb2Loc
               Cmp = Cmp + Work(ip0+nOrb2Loc*(iMO-1)+iMO)
            End Do
            Cmp = 2.0d0*Cmp
            Write(6,'(A,I5,1X,F15.8)')
     &      'Component, Exp. Val.:',iComp,Cmp
            Do j = 1,nOrb2Loc-1
               Do i = j+1,nOrb2Loc
                  kij = ip0 + nOrb2Loc*(j-1) + i
                  kji = ip0 + nOrb2Loc*(i-1) + j
                  Tst = Work(kij) - Work(kji)
                  If (abs(Tst) .gt. 1.0d-14) Then
                     Write(6,*) 'GenerateB: broken symmetry!'
                     Write(6,*) '  Component: ',iComp
                     Write(6,*) '  i and j  : ',i,j
                     Write(6,*) '  Dij      : ',Work(kij)
                     Write(6,*) '  Dji      : ',Work(kji)
                     Write(6,*) '  Diff.    : ',Tst
                     Call SysAbendMsg('GenerateB','Broken symmetry!',
     &                                ' ')
                  End If
               End Do
            End Do
         End Do
      End If

      End
