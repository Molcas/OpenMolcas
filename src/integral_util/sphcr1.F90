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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************
      SubRoutine SphCr1(Win,ijkla,                                      &
     &                  Scrt,nScrt,                                     &
     &                  Coeff3,kCar,kSph,Tr3,Pr3,                       &
     &                  Coeff4,lCar,lSph,Tr4,Pr4,Wout,mcd)
!***********************************************************************
!                                                                      *
! Object : to transform the two-electron integrals from cartesian      *
!          gaussians to real spherical harmonic gaussians.             *
!                                                                      *
!          Observe that most of the time Win and Wout will overlap.    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to back projection to cartesian gaussians,      *
!             January '92.                                             *
!***********************************************************************
      Implicit None
      Integer ijkla, nScrt, kCar,kSph,lCar,lSph,mcd
      Real*8 Win(ijkla*kSph*lSph), Scrt(nScrt),                         &
     &       Coeff3(kCar,kCar), Coeff4(lCar,lCar),                      &
     &       Wout(mcd*ijkla)
      Logical Tr3, Pr3, Tr4, Pr4
!
!     Call RecPrt(' In SphCr1: P(AB|CD) ',' ',Win,ijkla,kSph*lSph)
      If (Tr3.and.Tr4) Then
!        Call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
!--------Starting with IJKL,AB,CD transforming to d,IJKL,AB,C
!        Call xxDGeMul(Coeff4,lCar,'N',
!    &               Win,ijkla*kSph,'T',
!    &               Scrt,lCar,
!    &               lCar,lSph,ijkla*kSph)
         Call NTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kSph)
!
!        Call RecPrt(' In SphCr: P(AB|Cd) ',' ',Scrt,lCar*ijkla,kSph)
!        Call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
!--------Transform d,IJKL,AB,C to cd,IJKL,AB
!        Call xxDGeMul(Coeff3,kCar,'N',
!    &               Scrt,lCar*ijkla,'T',
!    &               Wout,kCar,
!    &               kCar,kSph,lCar*ijkla)
         Call NTMul(Coeff3,Scrt,Wout,kCar,kSph,lCar*ijkla)
      Else If (Tr4) Then
!        Call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
!--------Starting with IJKL,AB,cD transforming to d,IJKL,AB,c
!        Call xxDGeMul(Coeff4,lCar,'N',
!    &               Win,ijkla*kCar,'T',
!    &               Scrt,lCar,
!    &               lCar,lSph,ijkla*kCar)
         Call NTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)
!--------Transpose d,IJKL,AB,c to cd,IJKL,AB
         Call DGeTMO(Scrt,lCar*ijkla,lCar*ijkla,kCar,Wout,kCar)
      Else If (Tr3) Then
!--------Transpose IJKL,AB,C,d to d,IJKL,AB,C
         Call DGeTMO(Win,ijkla*kSph,ijkla*kSph,lCar,Scrt,lCar)
!
!        Call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
!        Transform d,IJKL,AB,c to cd,IJKL,AB
!        Call xxDGeMul(Coeff3,kCar,'N',
!    &               Scrt,lCar*ijkla,'T',
!    &               Wout,kCar,
!    &               kCar,kSph,lCar*ijkla)
         Call NTMul(Coeff3,Scrt,Wout,kCar,kSph,lCar*ijkla)
      Else
!---------Transpose IJKL,AB,cd to cd,IJKL,AB
          If (kCar*lCar.ne.1) Then
             call dcopy_(ijkla*kCar*lCar,Win,1,Scrt,1)
             Call DGeTMO(Scrt,ijkla,ijkla,kCar*lCar,Wout,kCar*lCar)
          Else
             call dcopy_(ijkla*kCar*lCar,Win,1,Scrt,1)
             call dcopy_(ijkla*kCar*lCar,Scrt,1,Wout,1)
          End If
      End If
!
!     Call RecPrt(' In SphCr1: P(AB|cd)  ',' ',Wout,mcd,ijkla)
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Pr3)
         Call Unused_logical(Pr4)
      End If
      End SubRoutine SphCr1
