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
      Subroutine AppFld(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Logical NonEq
*
      nCavSph=(lMax+1)**2
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
*
      Return
      End
      Subroutine AppFld_(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
      Logical NonEq
*
*     Statement function
*
      f(Eps,l)=(DBLE(1+l)*(Eps-One))/(DBLE(1+l)*Eps+DBLE(l))
*
      iRout=2
      iPrint=nPrint(iRout)
*
      If (iPrint.ge.99) Call RecPrt('Multipole Moments',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
*
*-----Backtransform from cartesian to spherical harmonics
*
      Call Tranca(Cavxyz,Cavsph,lmax,.True.)
      If (iPrint.ge.99) Call RecPrt(' CavSph',' ',
     &                              Cavsph,(lMax+1)**2,1)
*
*-----Evaluate the electric field components at the origin.
*     This is identical to the charge distribution on the
*     boundary of the cavity!
*
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
*
*-----Transform electric field components from spherical harmonics
*     to cartesians.
*
      Call Tranca(Cavxyz,Cavsph,lmax,.False.)
*
      If (iPrint.ge.99) Call RecPrt('Electric Field',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
*
      Return
      End
      Subroutine AppFld_NonEq_1(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Logical NonEq
*
      nCavSph=(lMax+1)**2
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_1(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
*
      Return
      End
      Subroutine AppFld_1(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
      Logical NonEq
*
*     Statement function
*
      f(Eps,l)=(DBLE(1+l)*(Eps-One))/(DBLE(1+l)*Eps+DBLE(l))
*
      iRout=2
      iPrint=nPrint(iRout)
*
      If (iPrint.ge.99) Call RecPrt('Multipole Moments',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
*
*-----Backtransform from cartesian to spherical harmonics
*
      Call Tranca(Cavxyz,Cavsph,lmax,.True.)
      If (iPrint.ge.99) Call RecPrt(' CavSph',' ',
     &                              Cavsph,(lMax+1)**2,1)
*
*-----Evaluate the electric field components at the origin.
*     This is identical to the charge distribution on the
*     boundary of the cavity!
*
      ip = 1
      Do l=0,lmax
         rinv=One/radius**(2*l+1)
         fact=F(Eps,l)*(One-F(EpsInf,l)/F(Eps,l))**2
         rpoti=rinv*fact*DblFac(2*l-1)
         Call DScal_(2*l+1,rpoti,Cavsph(ip),1)
         ip = ip + 2*l+1
      End Do
*
*-----Transform electric field components from spherical harmonics
*     to cartesians.
*
      Call Tranca(Cavxyz,Cavsph,lmax,.False.)
*
      If (iPrint.ge.99) Call RecPrt('Electric Field',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(NonEq)
      End
      Subroutine AppFld_NonEq_2(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Logical NonEq
*
      nCavSph=(lMax+1)**2
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_2(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
*
      Return
      End
      Subroutine AppFld_2(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
      Logical NonEq
*
*     Statement function
*
      f(Eps,l)=(DBLE(1+l)*(Eps-One))/(DBLE(1+l)*Eps+DBLE(l))
*
      iRout=2
      iPrint=nPrint(iRout)
*
      If (iPrint.ge.99) Call RecPrt('Multipole Moments',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
*
*-----Backtransform from cartesian to spherical harmonics
*
      Call Tranca(Cavxyz,Cavsph,lmax,.True.)
      If (iPrint.ge.99) Call RecPrt(' CavSph',' ',
     &                              Cavsph,(lMax+1)**2,1)
*
*-----Evaluate the electric field components at the origin.
*     This is identical to the charge distribution on the
*     boundary of the cavity!
*
      ip = 1
      Do l=0,lmax
         rinv=One/radius**(2*l+1)
         fact = F(Eps,l) - F(EpsInf,l)
     &        - (F(EpsInf,l)-F(EpsInf,l)**2/F(Eps,l))
         rpoti=rinv*fact*DblFac(2*l-1)
         Call DScal_(2*l+1,rpoti,Cavsph(ip),1)
         ip = ip + 2*l+1
      End Do
*
*-----Transform electric field components from spherical harmonics
*     to cartesians.
*
      Call Tranca(Cavxyz,Cavsph,lmax,.False.)
*
      If (iPrint.ge.99) Call RecPrt('Electric Field',' ',Cavxyz,
     &                              (lMax+1)*(lMax+2)*(lMax+3)/6,1)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(NonEq)
      End
