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
      SubRoutine CrSph2(Win,ijkl,ncd,
     &                  Scrt,nScrt,
     &                  Coeff1,iCar,iSph,Tr1,Pr1,
     &                  Coeff2,jCar,jSph,Tr2,Pr2,Wout,mab)
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
      Real*8 Win(ijkl*ncd*iCar*jCar), Scrt(nScrt),
     &       Coeff1(iCar,iCar), Coeff2(jCar,jCar),
     &       Wout(ijkl*ncd*mab)
      Logical Tr1, Pr1, Tr2, Pr2
*
      iRout = 60
      iPrint = nPrint(iRout)
*     Call RecPrt(' In CrSph2: (ab|CD) ',' ',Win,ncd*ijkl,iCar*jCar)
      If (Tr1.and.Tr2) Then
*        Call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
*        Starting with CD,IJKL,a,b transform to B,CD,IJKL,a
*        Call xxDGeMul(Coeff2,jCar,'T',
*    &               Win,ijkl*ncd*iCar,'T',
*    &               Scrt,jSph,
*    &               jSph,jCar,ijkl*ncd*iCar)
         Call TTMul(Coeff2,Win,Scrt,jCar,jSph,ijkl*ncd*iCar)
*
*        Call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
*        Transform B,CD,IJKL,a to AB,CD,IJKL
*        Call xxDGeMul(Coeff1,iCar,'T',
*    &               Scrt,jSph*ncd*ijkl,'T',
*    &               Wout,iSph,
*    &               iSph,iCar,jSph*ncd*ijkl)
         Call TTMul(Coeff1,Scrt,Wout,iCar,iSph,jSph*ncd*ijkl)
*        Transpose AB,CD,IJKL to IJKL,AB,CD
         call dcopy_(mab*ncd*ijkl,Wout,1,Scrt,1)
         Call dGeTMO(Scrt,mab*ncd,mab*ncd,ijkl,Wout,ijkl)
      Else If (Tr2) Then
*        Call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
*        Starting with CD,IJKL,a,b transform to B,CD,IJKL,a
*        Call xxDGeMul(Coeff2,jCar,'T',
*    &               Win,ncd*ijkl*iCar,'T',
*    &               Scrt,jSph,
*    &               jSph,jCar,ncd*ijkl*iCar)
         Call TTMul(Coeff2,Win,Scrt,jCar,jSph,ijkl*ncd*iCar)
*        Transpose B,CD,IJKL,a to IJKL,aB,CD
         Call DGeTMO(Scrt,jSph*ncd,jSph*ncd,ijkl*iCar,Wout,ijkl*iCar)
      Else If (Tr1) Then
*        Transpose CD,IJKL,a,b to b,CD,IJKL,a
         Call DGeTMO(Win,ncd*ijkl*iCar,ncd*ijkl*iCar,jCar,Scrt,jCar)
*
*        Call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
*        Transform b,CD,IJKL,a to Ab,CD,IJKL
*        Call xxDGeMul(Coeff1,iCar,'T',
*    &               Scrt,jCar*ncd*ijkl,'T',
*    &               Wout,iSph,
*    &               iSph,iCar,jCar*ncd*ijkl)
         Call TTMul(Coeff1,Scrt,Wout,iCar,iSph,jCar*ncd*ijkl)
*        Transpose Ab,CD,IJKL to IJKL,Ab,CD
         call dcopy_(iSph*jCar*ncd*ijkl,Wout,1,Scrt,1)
         Call DGeTMO(Scrt,iSph*jCar*ncd,iSph*jCar*ncd,ijkl,Wout,ijkl)
      Else
*         Transpose CD,IJKL,a,b to IJKL,ab,CD
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
*     Call RecPrt(' In CrSph2: (AB|CD)',' ',Wout,ijkl,ncd*mab)
*     Call GetMem(' Exit CrSph2','CHECK','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Pr1)
         Call Unused_logical(Pr2)
      End If
      End
