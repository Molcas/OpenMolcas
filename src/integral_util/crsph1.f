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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine CrSph1(Win,ijkla,
     &                  Scrt,nScrt,
     &                  Coeff3,kCar,kSph,Tr3,Pr3,
     &                  Coeff4,lCar,lSph,Tr4,Pr4,Wout,mcd)
************************************************************************
*                                                                      *
* Object : to transform the two-electron integrals from cartesian      *
*          gaussians to spherical gaussians.                           *
*                                                                      *
*          Observe that most of the time Win and Wout will overlap.    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
      Real*8 Win(ijkla*kCar*lCar), Scrt(nScrt),
     &       Coeff3(kCar,kCar), Coeff4(lCar,lCar),
     &       Wout(mcd*ijkla)
      Logical Tr3, Pr3, Tr4, Pr4
*
      If (Tr3.and.Tr4) Then
*        Call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
*        Starting with IJKL,a,cd transforming to D,IJKL,a,c
*        Call xxDGeMul(Coeff4,lCar,'T',
*    &               Win,ijkla*kCar,'T',
*    &               Scrt,lSph,
*    &               lSph,lCar,ijkla*kCar)
         Call TTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)
*
*        Call RecPrt(' In CrSph: (a0|cD) ',' ',Scrt,lSph*ijkla,kCar)
*        Call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
*        Transform D,IJKL,a,c to CD,IJKL,a
*        Call xxDGeMul(Coeff3,kCar,'T',
*    &               Scrt,lSph*ijkla,'T',
*    &               Wout,kSph,
*    &               kSph,kCar,lSph*ijkla)
         Call TTMul(Coeff3,Scrt,Wout,kCar,kSph,lSph*ijkla)
      Else If (Tr4) Then
*        Call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
*        Starting with IJKL,a,cd transforming to D,IJKL,a,c
*        Call xxDGeMul(Coeff4,lCar,'T',
*    &               Win,ijkla*kCar,'T',
*    &               Scrt,lSph,
*    &               lSph,lCar,ijkla*kCar)
         Call TTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)
*        Transpose D,IJKL,a,c to cD,IJKL,a
         Call DGeTMO(Scrt,lSph*ijkla,lSph*ijkla,kCar,Wout,kCar)
      Else If (Tr3) Then
*        Transpose IJKL,a,c,d to d,IJKL,a,c
         Call DGeTMO(Win,ijkla*kCar,ijkla*kCar,lCar,Scrt,lCar)
*
*        Call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
*        Transform D,IJKL,a,c to CD,IJKL,a
*        Call xxDGeMul(Coeff3,kCar,'T',
*    &               Scrt,lSph*ijkla,'T',
*    &               Wout,kSph,
*    &               kSph,kCar,lSph*ijkla)
         Call TTMul(Coeff3,Scrt,Wout,kCar,kSph,lSph*ijkla)
      Else
*         Transpose IJKL,a,c,d to c,d,IJKL,a
          If (kCar*lCar.ne.1) Then
             call dcopy_(ijkla*kCar*lCar,Win,1,Scrt,1)
             Call DGeTMO(Scrt,ijkla,ijkla,kCar*lCar,Wout,kCar*lCar)
          Else
             call dcopy_(ijkla*kCar*lCar,Win,1,Scrt,1)
             call dcopy_(ijkla*kCar*lCar,Scrt,1,Wout,1)
          End If
      End If
*
*     Call RecPrt(' In CrSph1: (a0|CD)  ',' ',Wout,mcd,ijkla)
*     Call GetMem(' Exit CrSph1','CHECK','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Pr3)
         Call Unused_logical(Pr4)
      End If
      End
