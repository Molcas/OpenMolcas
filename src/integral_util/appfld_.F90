!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine AppFld_(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)
      use Constants, only: One, Two
      Implicit None
      Integer lMax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
      Real*8 Radius, Eps, EpsInf
      Logical NonEq

      Integer l, ip
      Real*8 RInv, Fact, rPoti, DblFac, F
!
!     Statement function
!
      f(Eps,l)=(DBLE(1+l)*(Eps-One))/(DBLE(1+l)*Eps+DBLE(l))
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Multipole Moments',' ',Cavxyz,                       &
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif
!
!-----Backtransform from cartesian to spherical harmonics
!
      Call Tranca(Cavxyz,Cavsph,lmax,.True.)
#ifdef _DEBUGPRINT_
      Call RecPrt(' CavSph',' ',Cavsph,(lMax+1)**2,1)
#endif
!
!-----Evaluate the electric field components at the origin.
!     This is identical to the charge distribution on the
!     boundary of the cavity!
!
      ip = 1
      If (NonEq) Then
         Do l=0,lmax
            rinv=One/radius**(2*l+1)
            fact=(Two*f(EpsInf,l)-f(EpsInf,l)**2/f(Eps,l))
            rpoti=rinv*fact*DblFac(2*l-1)
            Call DScal_(2*l+1,rpoti,Cavsph(ip),1)
            ip = ip + 2*l+1
         End Do
      Else
         Do l=0,lmax
            rinv=One/radius**(2*l+1)
            fact=f(Eps,l)
            rpoti=rinv*fact*DblFac(2*l-1)
            Call DScal_(2*l+1,rpoti,Cavsph(ip),1)
            ip = ip + 2*l+1
         End Do
      End If
!
!-----Transform electric field components from spherical harmonics
!     to cartesians.
!
      Call Tranca(Cavxyz,Cavsph,lmax,.False.)
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Electric Field',' ',Cavxyz,                          &
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif
!
      Return
      End Subroutine AppFld_
