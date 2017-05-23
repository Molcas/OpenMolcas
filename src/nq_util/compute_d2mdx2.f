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
      Subroutine Compute_d2Mdx2(ZA,nAtoms,iAtom,iCar,dTdRAi,
     &                                    jAtom,jCar,dTdRaj,d2Mdx2)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 ZA(nAtoms),d2Mdx2(3,3)
*                                                                      *
************************************************************************
*                                                                      *
      Call FZero(d2Mdx2,9)
      Do kAtom = 1, nAtoms
         ZB=ZA(kAtom)
         If (kAtom.eq.iAtom) Then
            tmpi=One-dTdRAi
         Else
            tmpi=   -dtdRAi
         End If
         If (kAtom.eq.jAtom) Then
            tmpj=One-dTdRAi
         Else
            tmpj=   -dtdRAi
        End If
*
        If (iCar.eq.1.and.jCar.eq.1) Then
           d2Mdx2(2,2) = d2Mdx2(2,2) + Two*ZB*tmpi*tmpj
           d2Mdx2(3,3) = d2Mdx2(3,3) + Two*ZB*tmpi*tmpj
        End If
        If (iCar.eq.1.and.jCar.eq.2) Then
           d2Mdx2(1,2) = d2Mdx2(1,2) -     ZB*tmpi*tmpj
           d2Mdx2(2,1) = d2Mdx2(2,1) -     ZB*tmpj*tmpi
         End If
         If (iCar.eq.1.and.jCar.eq.3) Then
            d2Mdx2(1,3) = d2Mdx2(1,3) -     ZB*tmpi*tmpj
            d2Mdx2(3,1) = d2Mdx2(3,1) -     ZB*tmpj*tmpi
         End If
         If (iCar.eq.2.and.jCar.eq.1) Then
            d2Mdx2(1,2) = d2Mdx2(1,2) -     ZB*tmpj*tmpi
            d2Mdx2(2,1) = d2Mdx2(2,1) -     ZB*tmpi*tmpj
         End If
         If (iCar.eq.2.and.jCar.eq.2) Then
            d2Mdx2(1,1) = d2Mdx2(1,1) + Two*ZB*tmpi*tmpj
            d2Mdx2(3,3) = d2Mdx2(3,3) + Two*ZB*tmpi*tmpj
         End If
         If (iCar.eq.2.and.jCar.eq.3) Then
            d2Mdx2(2,3) = d2Mdx2(2,3) -     ZB*tmpi*tmpj
            d2Mdx2(3,2) = d2Mdx2(3,2) -     ZB*tmpj*tmpi
         End If
         If (iCar.eq.3.and.iCar.eq.1) Then
            d2Mdx2(1,3) = d2Mdx2(1,3) -     ZB*tmpj*tmpi
            d2Mdx2(3,1) = d2Mdx2(3,1) -     ZB*tmpi*tmpj
         End If
         If (iCar.eq.3.and.jCar.eq.2) Then
            d2Mdx2(2,3) = d2Mdx2(2,3) -     ZB*tmpj*tmpi
            d2Mdx2(3,2) = d2Mdx2(3,2) -     ZB*tmpi*tmpj
         End If
         If (iCar.eq.3.and.jCar.eq.3) Then
            d2Mdx2(1,1) = d2Mdx2(1,1) + Two*ZB*tmpi*tmpj
            d2Mdx2(2,2) = d2Mdx2(2,2) + Two*ZB*tmpi*tmpj
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(dTdRaj)
      End
