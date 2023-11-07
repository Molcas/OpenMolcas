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
      Subroutine AppFld(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer lMax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Real*8 Radius, Eps, EpsInf
      Logical NonEq
!
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
!
      Return
      End Subroutine AppFld

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
      Call RecPrt('Multipole Moments',' ',Cavxyz,
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
      Call RecPrt('Electric Field',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif
!
      Return
      End Subroutine AppFld_

      Subroutine AppFld_NonEq_1(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer lMax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Real*8 Radius, Eps, EpsInf
      Logical NonEq
!
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_1(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
!
      Return
      End Subroutine AppFld_NonEq_1

      Subroutine AppFld_1(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)
      use Constants, only: One
      Implicit None
      Integer lMax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
      Real*8 Radius, Eps, EpsInf
      Logical NonEq

      Integer l, ip
      Real*8 f, rInv, rPoti, Fact, DblFac
!
!     Statement function
!
      f(Eps,l)=(DBLE(1+l)*(Eps-One))/(DBLE(1+l)*Eps+DBLE(l))
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Multipole Moments',' ',Cavxyz,
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
      Do l=0,lmax
         rinv=One/radius**(2*l+1)
         fact=F(Eps,l)*(One-F(EpsInf,l)/F(Eps,l))**2
         rpoti=rinv*fact*DblFac(2*l-1)
         Call DScal_(2*l+1,rpoti,Cavsph(ip),1)
         ip = ip + 2*l+1
      End Do
!
!-----Transform electric field components from spherical harmonics
!     to cartesians.
!
      Call Tranca(Cavxyz,Cavsph,lmax,.False.)
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Electric Field',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif
!
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_logical(NonEq)
      End Subroutine AppFld_1

      Subroutine AppFld_NonEq_2(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer lmax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Real*8 radius, EPS, EPSInf
      Logical NonEq
!
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_2(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
!
      Return
      End Subroutine AppFld_NonEq_2

      Subroutine AppFld_2(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)
      use Constants, only: One
      Implicit None
      Integer lMax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
      Real*8 Radius, Eps, EpsInf
      Logical NonEq

      Integer ip, l
      Real*8 rInv, Fact, rPoti, DblFac, F
!
!     Statement function
!
      f(Eps,l)=(DBLE(1+l)*(Eps-One))/(DBLE(1+l)*Eps+DBLE(l))
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Multipole Moments',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif
!
!-----Backtransform from cartesian to spherical harmonics
!
      Call Tranca(Cavxyz,Cavsph,lmax,.True.)
#ifdef _DEBUGPRINT_
      Call RecPrt(' CavSph',' ',
     &                              Cavsph,(lMax+1)**2,1)
#endif
!
!-----Evaluate the electric field components at the origin.
!     This is identical to the charge distribution on the
!     boundary of the cavity!
!
      ip = 1
      Do l=0,lmax
         rinv=One/radius**(2*l+1)
         fact = F(Eps,l) - F(EpsInf,l)
     &        - (F(EpsInf,l)-F(EpsInf,l)**2/F(Eps,l))
         rpoti=rinv*fact*DblFac(2*l-1)
         Call DScal_(2*l+1,rpoti,Cavsph(ip),1)
         ip = ip + 2*l+1
      End Do
!
!-----Transform electric field components from spherical harmonics
!     to cartesians.
!
      Call Tranca(Cavxyz,Cavsph,lmax,.False.)
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Electric Field',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif
!
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_logical(NonEq)
      End Subroutine AppFld_2
