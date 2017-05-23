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
      Subroutine ZMatConv(LuWr,iAtom,iErr)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension u1(3),u2(3),u3(3),u4(3),vj(3),vp(3)
#include "g_zmatconv.fh"

      iErr = 0

      torad = 3.1415926536d0 / 180.0d0
      dCBiAtom = COS(torad*Zmat(iAtom,2))
      dSBiAtom = SIN(torad*Zmat(iAtom,2))
      dCTiAtom = COS(torad*Zmat(iAtom,3))
      dSTiAtom = SIN(torad*Zmat(iAtom,3))
      If (ABS(dCBiAtom).LT.1.0d-10) dCBiAtom = 0.0d0
      If (ABS(dSBiAtom).LT.1.0d-10) dSBiAtom = 0.0d0
      If (ABS(dCTiAtom).LT.1.0d-10) dCTiAtom = 0.0d0
      If (ABS(dSTiAtom).LT.1.0d-10) dSTiAtom = 0.0d0

      Call vec(1.0d-06,u1,iZmat(iAtom,2),iZmat(iAtom,3),iErr)
! Vet. u1(NB-NT)
      If (iErr.NE.0) GoTo 9990 ! Vettore u1 nullo
      Call vec(1.0d-06,u2,iZmat(iAtom,1),iZmat(iAtom,2),iErr)
! Vet. u2(NA,NB)
      If (iErr.NE.0) GoTo 9990 ! Vettore u2 nullo
      arg = 1.0d0 - ( u1(1)*u2(1) + u1(2)*u2(2) + u1(3)*u2(3) )**2
      If (arg.LT.0.0d0) GoTo 9990 ! u1 e u2 allineati
      r = sqrt(arg) ! (u1.u2)^0.5
      If (r .LT. 1.0d-6) GoTo 9990 ! r piccolo

      Call crprod(u1,u2,vp) ! Vettore piano perpendicolare u1 x u2
      Do i=1,3
        u3(i) = vp(i)/r
      EndDo
      Call crprod(u3,u2,u4)
      Do i=1,3
        vj(i) = Zmat(iAtom,1)*(-u2(i)*dCBiAtom+u4(i)*dSBiAtom*dCTiAtom
     &                                       +u3(i)*dSBiAtom*dSTiAtom)
        Coords(iAtom,i) = Coords(iZmat(iAtom,1),i) + vj(i)
      EndDo

      GoTo 9999

9990  iErr=1
      Write(LuWr,*) ' [Z-Mat_Conv] Incipient floating point error ',
     &'detected for atom ',iAtom

9999  Return
      End
