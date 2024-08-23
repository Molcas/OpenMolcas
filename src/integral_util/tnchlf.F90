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
! Copyright (C) 1990,1992,1994,1996, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************
      Subroutine Tnchlf(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,        &
     &                  lZeta,nVec,IncVec,A1,A2,A3,Indij)
!***********************************************************************
!                                                                      *
! Object: to do a half transformation. The loop over the two matrix-   *
!         matrix multiplications is segmented such that the end of the *
!         intermediate matrix will not push the start of the same out  *
!         from the cache.                                              *
!                                                                      *
! Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!                                                                      *
!             Modified to decontraction May 1996, by R. Lindh          *
!***********************************************************************
      use Constants, only: Zero
      Implicit None
      Integer nCntr1,nPrm1,nCntr2,nPrm2,lZeta,nVec,IncVec
      Real*8 Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2),                &
     &       A1(nVec,nCntr1,nCntr2), A2(nPrm2,IncVec,nCntr1),           &
     &       A3(lZeta,nVec)
      Integer Indij(lZeta)

      Logical Seg1, Seg2
      Integer iPrm1, iCntr1, iPrm2, iCntr2, iiVec, mVec, iVec, iZeta
!
!-----Check if the basis set is segmented
!
      Seg1=.False.
      Do iPrm1 = nPrm1, 1, -1
         Do iCntr1 = nCntr1, 1, -1
            If (Coeff1(iPrm1,iCntr1).eq.Zero) Then
               Seg1=.True.
               Go To 10
            End If
         End Do
      End Do
 10   Continue
!
      Seg2=.False.
      Do iPrm2 = nPrm2, 1, -1
         Do iCntr2 = nCntr2, 1, -1
            If (Coeff2(iPrm2,iCntr2).eq.Zero) Then
               Seg2=.True.
               Go To 20
            End If
         End Do
      End Do
 20   Continue
!
!-----Set output matrix to zero
!
      Call FZero(A3,nVec*lZeta)
!
!-----Loop sectioning
!
      Do iiVec = 1, nVec, IncVec
         mVec = Min(IncVec,nVec-iiVec+1)
!--------Set intermediate matrix to zero
         Call FZero(A2,nPrm2*IncVec*nCntr1)
!
         If (Seg2) Then
!
!-----First quarter transformation, (x,AB) -> (b,x,a)
!
      Do iPrm2 = 1, nPrm2
         Do iCntr2 = 1, nCntr2
!-----------Check for zero due to segmented basis
            If (Abs(Coeff2(iPrm2,iCntr2)).gt.Zero) Then
               Do iCntr1 = 1, nCntr1
                  Do iVec = 1, mVec
                     A2(iPrm2,iVec,iCntr1) = A2(iPrm2,iVec,iCntr1) +    &
     &                 Coeff2(iPrm2,iCntr2)                             &
     &                 *A1(iVec+iiVec-1,iCntr1,iCntr2)
                  End Do
               End Do
            End If
         End Do
      End Do
!
         Else    ! Seg2
!
      Do iPrm2 = 1, nPrm2
         Do iCntr2 = 1, nCntr2
            Do iCntr1 = 1, nCntr1
               Do iVec = 1, mVec
                  A2(iPrm2,iVec,iCntr1) = A2(iPrm2,iVec,iCntr1) +       &
     &              Coeff2(iPrm2,iCntr2)                                &
     &              *A1(iVec+iiVec-1,iCntr1,iCntr2)
               End Do
            End Do
         End Do
      End Do
!
         End If   ! Seg2
!
         If (Seg1) Then
!
!-----Second quarter transformation
!
         Do iCntr1 = 1, nCntr1
            Do iZeta = 1, lZeta
               iPrm2=(Indij(iZeta)-1)/nPrm1 + 1
               iPrm1=Indij(iZeta)-(iPrm2-1)*nPrm1
!--------------Check for zero due to segmented basis
               If (Abs(Coeff1(iPrm1,iCntr1)).gt.Zero) Then
                  Do iVec = iiVec, iiVec+mVec-1
                     A3(iZeta,iVec) = A3(iZeta,iVec) +                  &
     &                 Coeff1(iPrm1,iCntr1)*                            &
     &                 A2(iPrm2,iVec-iiVec+1,iCntr1)
                  End Do ! iVec
               End If
            End Do       ! iZeta
         End Do          ! iCntr1
!
         Else
!
!-----Second quarter transformation
!
         Do iCntr1 = 1, nCntr1
            Do iZeta = 1, lZeta
               iPrm2=(Indij(iZeta)-1)/nPrm1 + 1
               iPrm1=Indij(iZeta)-(iPrm2-1)*nPrm1
               Do iVec = iiVec, iiVec+mVec-1
                  A3(iZeta,iVec) = A3(iZeta,iVec) +                     &
     &              Coeff1(iPrm1,iCntr1)*                               &
     &              A2(iPrm2,iVec-iiVec+1,iCntr1)
               End Do ! iVec
            End Do    ! iZeta
         End Do       ! iCntr1
!
      End If
!
!-----End of loop sectioning
!
      End Do
!
      Return
      End Subroutine Tnchlf
