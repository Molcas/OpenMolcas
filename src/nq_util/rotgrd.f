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
* Copyright (C) 2003, Roland Lindh                                     *
************************************************************************
      Subroutine RotGrd(RA,ZA,O,dOdx,d2Odx2,nAtoms,Do_Grad,Do_Hess)
************************************************************************
*                                                                      *
*     Object: Compute the principle axis system and to optionally      *
*             evaluate the gradient of the principle axis system.      *
*                                                                      *
*             See: B. G. Johnson et al., CPL, 220, 377 (1994).         *
*                                                                      *
*     Author: R. Lindh, Dept. of Chem. Phys., Univ. of Lund, Sweden.   *
*                                                                      *
*             Created on board M/S Polarlys on voyage from Tromso to   *
*             Trondheim, Sept. 2003.                                   *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 RA(3,nAtoms), ZA(nAtoms),
     &       O(3,3), dOdx(3,3,nAtoms,3), d2Odx2(3,3,nAtoms,3,nAtoms,3),
     &               dMdx(3,3), dMdy(3,3),
     &       EVal(3), T(3),
     &       Px(3,3), Py(3,3)
      Logical Do_Grad, Rot_Corr, Do_Hess
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('RotGrd: RA',' ',RA,3,nAtoms)
      Call RecPrt('RotGrd: ZA',' ',ZA,1,nAtoms)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the total charge
*
      Z_Tot=DDot_(nAtoms,[One],0,ZA,1)
*
*---- Form the center of nuclear charge
*
      Call Compute_T(Z_Tot,T,ZA,RA,nAtoms)
*
*---- Form O
*
      Call Compute_O(ZA,RA,nAtoms,Z_Tot,T,O,EVal)
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.Do_Grad) Return
*                                                                      *
************************************************************************
*                                                                      *
*     Turn off rotational correction to the gradient if the eigen
*     vectors are close to degeneracy!
*
      Rot_Corr=.True.
      If (Abs(Eval(1)-Eval(2))/(EVal(1)+EVal(2)).lt.0.001D0) Then
         Write (6,*) 'Rotational correction to the DFT gradient is '
     &             //'turned off due to close-to-degeneracy problems!'
         Rot_Corr=.False.
      End If
      If (Abs(Eval(1)-Eval(3))/(EVal(1)+EVal(3)).lt.0.001D0) Then
         Write (6,*) 'Rotational correction to the DFT gradient is '
     &             //'turned off due to close-to-degeneracy problems!'
         Rot_Corr=.False.
      End If
      If (Abs(Eval(2)-Eval(3))/(EVal(2)+EVal(3)).lt.0.001D0) Then
         Write (6,*) 'Rotational correction to the DFT gradient is '
     &             //'turned off due to close-to-degeneracy problems!'
         Rot_Corr=.False.
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the gradient
*
      Do iAtom = 1, nAtoms
         dTdRAi=ZA(iAtom)/Z_Tot
         Do iCar = 1, 3
*
*---------- Form dO/dx
*
            Call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,
     &                        iAtom,iCar,dTdRAi,dMdx,
     &                        dOdx(1,1,iAtom,iCar),Px)
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.Do_Hess) Return
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the Hessian
*
      Do iAtom = 1, nAtoms
         dTdRAi=ZA(iAtom)/Z_Tot
         Do iCar = 1, 3
            Call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,
     &                        iAtom,iCar,dTdRAi,dMdx,
     &                        dOdx(1,1,iAtom,iCar),Px)
*
            Do jAtom = 1, iAtom
               dTdRAj=ZA(jAtom)/Z_Tot
               jCar_Max = 3
               If (iAtom.eq.jAtom) jCar_Max = iCar
               Do jCar = 1, jCar_Max
                  Call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,
     &                              jAtom,jCar,dTdRAj,dMdy,
     &                              dOdx(1,1,jAtom,jCar),Py)
*
*---------------- Form d2O/dx2
*
                  Call Compute_d2Odx2(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,
     &                                iAtom,iCar,dTdRAi,dMdx,
     &                                dOdx(1,1,iAtom,iCar),Px,
     &                                jAtom,jCar,dTdRAj,dMdy,
     &                                dOdx(1,1,jAtom,jCar),Py,
     &                                d2Odx2(1,1,iAtom,iCar,jAtom,jCar))
*
                  If (iAtom.ne.jAtom .or.
     &                (iAtom.eq.jAtom.and.iCar.ne.jCar)) Then
                      call dcopy_(9,d2Odx2(1,1,iAtom,iCar,jAtom,jCar),1,
     &                             d2Odx2(1,1,jAtom,jCar,iAtom,iCar),1)
                  End If
*
               End Do  ! jCar
            End Do     ! iCar
*
         End Do        ! jAtom
      End Do           ! iAtom
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
