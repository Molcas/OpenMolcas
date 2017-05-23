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
      SubRoutine ComputeFuncB2(nOrb2Loc,ipLbl,nComp,Functional,Debug)
C
C     Author: T.B. Pedersen
C
C     Purpose: compute Boys localisation functional B2.
C
      Implicit Real*8 (a-h,o-z)
      Integer ipLbl(nComp)
      Logical Debug
#include "WrkSpc.fh"

      Functional = 0.0d0
      Do iComp = 1,nComp
         ip0 = ipLbl(iComp) - 1
         Do i = 1,nOrb2Loc
            Functional = Functional + (Work(ip0+nOrb2Loc*(i-1)+i))**2
         End Do
      End Do

      If (Debug) Then
         Write(6,*)
         Write(6,*) 'In ComputeFuncB2'
         Write(6,*) '----------------'
         Write(6,*) 'Functional B2 = ',Functional
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
                     Write(6,*) 'ComputeFuncB2: broken symmetry!'
                     Write(6,*) '  Component: ',iComp
                     Write(6,*) '  i and j  : ',i,j
                     Write(6,*) '  Dij      : ',Work(kij)
                     Write(6,*) '  Dji      : ',Work(kji)
                     Write(6,*) '  Diff.    : ',Tst
                     Call SysAbendMsg('ComputeFuncB2',
     &                                'Broken symmetry!',' ')
                  End If
               End Do
            End Do
         End Do
      End If

      End
