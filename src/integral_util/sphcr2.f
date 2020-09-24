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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine SphCr2(Win,ijkl,ncd,
     &                  Scrt,nScrt,
     &                  Coeff1,iCar,iSph,Tr1,Pr1,
     &                  Coeff2,jCar,jSph,Tr2,Pr2,Wout,mab)
************************************************************************
*                                                                      *
* Object : to transform the two-electron integrals from cartesian      *
*          gaussians to real spherical harmonic gaussians.             *
*                                                                      *
*          Observe that most of the time Win and Wout will overlap.    *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DGEMM_   (ESSL)                                         *
*              RecPrt                                                  *
*              DGeTMO   (ESSL)                                         *
*              DCopy    (ESSL)                                         *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to back projection to cartesian gaussians,      *
*             January '92.                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
      Real*8 Win(ijkl*ncd*iSph*jSph), Scrt(nScrt),
     &       Coeff1(iCar,iCar), Coeff2(jCar,jCar),
     &       Wout(ijkl*ncd*mab)
      Logical Tr1, Pr1, Tr2, Pr2
*
      iRout = 60
      iPrint = nPrint(iRout)
*     Call RecPrt(' In SphCr2: P(AB|cd) ',' ',Win,ncd*ijkl,iSph*jSph)
      If (Tr1.and.Tr2) Then
*        Call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
*--------Starting with cd,IJKL,A,B transform to b,cd,IJKL,A
*        Call xxDGeMul(Coeff2,jCar,'N',
*    &               Win,ijkl*ncd*iSph,'T',
*    &               Scrt,jCar,
*    &               jCar,jSph,ijkl*ncd*iSph)
         Call NTMul(Coeff2,Win,Scrt,jCar,jSph,ijkl*ncd*iSph)
*
*        Call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
*--------Transform b,cd,IJKL,A to ab,cd,IJKL
*        Call xxDGeMul(Coeff1,iCar,'N',
*    &               Scrt,jCar*ncd*ijkl,'T',
*    &               Wout,iCar,
*    &               iCar,iSph,jCar*ncd*ijkl)
         Call NTMul(Coeff1,Scrt,Wout,iCar,iSph,jCar*ncd*ijkl)
*        Transpose ab,cd,IJKL to IJKL,ab,cd
         call dcopy_(mab*ncd*ijkl,Wout,1,Scrt,1)
         Call DGeTMO(Scrt,mab*ncd,mab*ncd,ijkl,Wout,ijkl)
      Else If (Tr2) Then
*        Call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
*--------Starting with cd,IJKL,a,B transform to b,cd,IJKL,a
*        Call xxDGeMul(Coeff2,jCar,'N',
*    &               Win,ncd*ijkl*iCar,'T',
*    &               Scrt,jCar,
*    &               jCar,jSph,ncd*ijkl*iCar)
         Call NTMul(Coeff2,Win,Scrt,jCar,jSph,ncd*ijkl*iCar)
*--------Transpose b,cd,IJKL,a to IJKL,ab,cd
         Call DGeTMO(Scrt,jCar*ncd,jCar*ncd,ijkl*iCar,Wout,ijkl*iCar)
      Else If (Tr1) Then
*--------Transpose cd,IJKL,A,b to b,cd,IJKL,A
         Call DGeTMO(Win,ncd*ijkl*iSph,ncd*ijkl*iSph,jCar,Scrt,jCar)
*
*        Call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
*--------Transform b,cd,IJKL,A to ab,cd,IJKL
*        Call xxDGeMul(Coeff1,iCar,'N',
*    &               Scrt,jCar*ncd*ijkl,'T',
*    &               Wout,iCar,
*    &               iCar,iSph,jCar*ncd*ijkl)
         Call NTMul(Coeff1,Scrt,Wout,iCar,iSph,jCar*ncd*ijkl)
*--------Transpose ab,cd,IJKL to IJKL,ab,cd
         call dcopy_(iCar*jCar*ncd*ijkl,Wout,1,Scrt,1)
         Call DGeTMO(Scrt,iCar*jCar*ncd,iCar*jCar*ncd,ijkl,Wout,ijkl)
      Else
*---------Transpose cd,IJKL,ab to IJKL,ab,cd
          If (ncd.ne.1) Then
             call dcopy_(ncd*ijkl*iCar*jCar,Win,1,Scrt,1)
             Call DGeTMO(Scrt,ncd,ncd,ijkl*iCar*jCar,Wout,
     &                                ijkl*iCar*jCar)
          Else
             call dcopy_(ncd*ijkl*iCar*jCar,Win,1,Scrt,1)
             call dcopy_(ncd*ijkl*iCar*jCar,Scrt,1,Wout,1)
          End If
      End If
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In SphCr2: P(ab|cd)',' ',Wout,ijkl,ncd*mab)
      End If
*     Call GetMem(' Exit SphCr2','CHECK','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Pr1)
         Call Unused_logical(Pr2)
      End If
      End
