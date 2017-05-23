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
      Subroutine Compute_dMdx(ZA,RA,nAtoms,T,iAtom,iCar,dTdRai,dMdx)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 ZA(nAtoms), RA(3,nAtoms), T(3), dMdx(3,3)
#ifdef _DEBUG_
      Real*8 M(3,3)
#endif
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Z_Tot=DDot_(nAtoms,One,0,ZA,1)
      delta=1.0D-4
      temp = RA(iCar,iAtom)
*
      RA(iCar,iAtom) = Temp + Delta
      Call Compute_M(ZA,nAtoms,RA,Z_Tot,T,M)
*
      RA(iCar,iAtom) = Temp - Delta
      Call Compute_M(ZA,nAtoms,RA,Z_Tot,T,dMdx)
*
      RA(iCar,iAtom) = Temp
*
      Do i = 1, 3
         Do j = 1, 3
            dMdx(i,j) = (M(i,j)-dMdx(i,j))/(2.0D0*Delta)
         End Do
      End Do
      Call RecPrt('dMdx(Numerical)',' ',dMdx,3,3)
#endif

      Call FZero(dMdx,3*3)
      Do jAtom = 1, nAtoms
         ZB=ZA(jAtom)
         If (iAtom.eq.jAtom) Then
            tmp=(One-dTdRAi)*ZB
         Else
            tmp=(   -dTdRAi)*ZB
         End If
*
         RTx=RA(1,jAtom)-T(1)
         RTy=RA(2,jAtom)-T(2)
         RTz=RA(3,jAtom)-T(3)
         If (iCar.eq.1) Then
            dMdx(2,2) = dMdx(2,2) + Two*tmp*RTx
            dMdx(3,3) = dMdx(3,3) + Two*tmp*RTx
            dMdx(1,2) = dMdx(1,2) -     tmp*RTy
            dMdx(2,1) = dMdx(2,1) -     RTy*tmp
            dMdx(1,3) = dMdx(1,3) -     tmp*RTz
            dMdx(3,1) = dMdx(3,1) -     RTz*tmp
         End If
         If (iCar.eq.2) Then
            dMdx(1,1) = dMdx(1,1) + Two*tmp*RTy
            dMdx(3,3) = dMdx(3,3) + Two*tmp*RTy
            dMdx(1,2) = dMdx(1,2) -     RTx*tmp
            dMdx(2,1) = dMdx(2,1) -     tmp*RTx
            dMdx(2,3) = dMdx(2,3) -     tmp*RTz
            dMdx(3,2) = dMdx(3,2) -     RTz*tmp
         End If
         If (iCar.eq.3) Then
            dMdx(1,1) = dMdx(1,1) + Two*tmp*RTz
            dMdx(2,2) = dMdx(2,2) + Two*tmp*RTz
            dMdx(1,3) = dMdx(1,3) -     RTx*tmp
            dMdx(3,1) = dMdx(3,1) -     tmp*RTx
            dMdx(2,3) = dMdx(2,3) -     RTy*tmp
            dMdx(3,2) = dMdx(3,2) -     tmp*RTy
         End If
      End Do
*
*     Remove noise
*
      Do i = 1, 3
         Do j = 1, 3
            If (Abs(dMdx(i,j)).lt.1.0D-14) dMdx(i,j)=Zero
         End Do
      End Do
#ifdef _DEBUG_
      Call RecPrt('dMdx',' ',dMdx,3,3)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
