!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************
      SubRoutine PrMtrx(Label,lOper,nComp,ip,Matrix)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, January '91                  *
!***********************************************************************
      Use Basis_Info, only: nBas
      use Gateway_global, only: PrPrt
      use Symmetry_Info, only: nIrrep
      Implicit None
      Integer nComp
      Character(LEN=*) Label
      Real*8 Matrix(*)
      Integer ip(nComp), lOper(nComp)

!     Local variables
      Character(LEN=80) Line
      Logical Type
      Integer iComp, ip1, iSmLbl, iIrrep, jIrrep
!
!
      Do iComp = 1, nComp
         ip1 = ip(iComp)
         iSmLbl = lOper(iComp)
         If (Prprt) iSmLbl = iAnd(1,iSmLbl)
         Type = .True.
         Do iIrrep = 0, nIrrep - 1
            If (nBas(iIrrep).le.0) Cycle
            Do jIrrep = 0, iIrrep
               If (nBas(jIrrep).le.0) Cycle
               If (iAnd(iSmLbl,2**iEor(iIrrep,jIrrep)).eq.0) Cycle
               If (Type) Then
                  Type = .False.
                  Write (6,*)
                  Write (6,*)
                  Write (6,'(A,A,A,I2)')
     &                  ' SO Integrals of type ', Label,' Component ',
     &                     iComp
               End If
               Line=''
               If (iIrrep.eq.jIrrep) Then
                  Write (Line,'(1X,A,I1)')
     &            ' Diagonal Symmetry Block ', iIrrep+1
                  Call TriPrt(Line,' ',
     &                 Matrix(ip1),nBas(iIrrep))
                  ip1 = ip1 + nBas(iIrrep)*(nBas(iIrrep)+1)/2
               Else
                  Write (Line,'(1X,A,I1,A,I1)')
     &            ' Off-diagonal Symmetry Block ',
     &            iIrrep+1, ',' , jIrrep+1
                  Call RecPrt(Line,' ',
     &                        Matrix(ip1),nBas(iIrrep),nBas(jIrrep))
                  ip1 = ip1 + nBas(iIrrep)*nBas(jIrrep)
               End If
            End Do
         End Do
      End Do
!
      Return
      End SubRoutine PrMtrx
