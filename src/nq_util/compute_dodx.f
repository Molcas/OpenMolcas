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
      Subroutine Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,
     &                        iAtom,iCar,dTdRAi,dMdx,dOdx,Px)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 ZA(nAtoms), RA(3,nAtoms), T(3), O(3,3), EVal(3),
     &       dOdx(3,3), dMdx(3,3), Px(3,3)
      Logical Rot_Corr
*---- Local arrays
      Real*8 OtMx(3,3), OtMxO(3,3)
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.Rot_Corr) Then
         Call FZero(dOdx,9)
         Return
      End If
*
*---- Differentiate the nuclear charge moment tensor, M.
*
      Call Compute_dMdx(ZA,RA,nAtoms,T,iAtom,iCar,dTdRAi,dMdx)
*
*
*     Form O(t) dMdRAi O
*
      Call DGEMM_('T','N',
     &            3,3,3,
     &            1.0d0,O,3,
     &            dMdx,3,
     &            0.0d0,OtMx,3)
      Call DGEMM_('N','N',
     &            3,3,3,
     &            1.0d0,OtMx,3,
     &            O,3,
     &            0.0d0,OtMxO,3)
*
*     Compute the off diagonal elements in Px
*
      If (Abs(OtMxO(2,3)+OtMxO(3,2)).lt.1.0D-14) Then
         If (abs(EVal(2)-EVal(3)).lt.1.0D-14) Then
             Alpha=One
         Else
             Alpha=Zero
         End If
      Else
         If (abs(EVal(2)-EVal(3)).lt.1.0D-14) Then
            Alpha=1.0D99
         Else
            Alpha =
     &        -(OtMxO(2,3)+OtMxO(3,2))/(Two*(EVal(2)-EVal(3)))
         End If
      End If
*
      If (Abs(OtMxO(1,3)+OtMxO(3,1)).lt.1.0D-14) Then
         If (abs(EVal(3)-EVal(1)).lt.1.0D-14) Then
            Beta=One
         Else
            Beta=Zero
         End If
      Else
         If (abs(EVal(3)-EVal(1)).lt.1.0D-14) Then
            Beta=1.0D99
         Else
            Beta  =
     &        -(OtMxO(1,3)+OtMxO(3,1))/(Two*(EVal(3)-EVal(1)))
         End If
      End If
      If (Abs(OtMxO(1,2)+OtMxO(2,1)).lt.1.0D-14) Then
         If (abs(EVal(1)-EVal(2)).lt.1.0D-14) Then
            Gamma = One
         Else
            Gamma = Zero
         End If
      Else
         If (abs(EVal(1)-EVal(2)).lt.1.0D-14) Then
            Gamma=1.0D99
         Else
            Gamma =
     &        -(OtMxO(1,2)+OtMxO(2,1))/(Two*(EVal(1)-EVal(2)))
         End If
      End If
*
      Call FZero(Px,9)
      Px(1,2)= Gamma
      Px(2,1)=-Gamma
      Px(1,3)=-Beta
      Px(3,1)= Beta
      Px(2,3)= Alpha
      Px(3,2)=-Alpha
*
*     Finally evaluate dO/dRAi
*
      Call DGEMM_('N','N',
     &            3,3,3,
     &            1.0d0,O,3,
     &            Px,3,
     &            0.0d0,dOdx(1,1),3)
#ifdef _DEBUGPRINT_
C     Call RecPrt('M','(3G20.10)',M,3,3)
C     Call RecPrt('RotGrd: O','(3G20.10)',O,3,3)
C     Call RecPrt('RotGrd: EVal',' ',EVal,3,1)
C     Call RecPrt('RotGrd: dMdx','(3G20.10)',dMdx,3,3)
C     Call RecPrt('RotGrd: OtMxO','(3G20.10)',OtMxO,3,3)
C     Write (*,*) 'A,B,G=',Alpha,Beta,Gamma
C     Call RecPrt('RotGrd: Px','(3G20.10)',Px,3,3)
C     Call RecPrt('dOdx','(3G20.10)',dOdx,3,3)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
