!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************
      Subroutine ContEI(I,mDeg,ESIT,ix,iy,iz,temp)
      implicit None
      Integer mDeg, ix, iy, iz
      Integer I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg)
      Real*8 ESIT((mDeg+1)*(mDeg+2)/2), Temp

      Integer n, ip, a, b, c
!
! Purpose: Express the interaction tensor, defined by the
! quantities T(a,b,c) as functions of the vector R=(x,y,z),
! where a,b, and c are nonnegative integers and
! T(a,b,c)=((d/dx)**a)((d/dy)**b)((d/dz)**c) 1/R, in terms
! of a polynomial:
! T(a,b,c)=
!  (sum over p,q,r of I(a,b,c,p,q,r) x**p y**q z**r)/(R**(2*n+1)),
! where n=a+b+c.
! The polynomial coefficients are integers, and are 0 unless
! p+q+r=n.
!
!
!----- Statement function
!
!      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
!
! write out only elements with a>=b>=c. The others are obtained
! by index permutation.
! This restriction has been removed! (Roland Lindh)
!     n=mDeg
!     do 200 a=n,0,-1
!     do 200 b=n-a,0,-1
!     c=n-a-b
!     write(*,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,
!    &     Ind(n,a,c)
!     do 150 p=n,0,-1
!     do 150 q=n-p,0,-1
!     r=n-p-q
!     coef=I(a,b,c,p,q,r)
!     if(coef.eq.0) goto 150
!     write(*,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)')
!    &  coef,p,q,r,Ind(n,p,r)
!150  continue
!200  continue
!
!     Write (*,*) ' Temp,ix,iy,iz=',temp,ix,iy,iz
      n=mDeg
      ip = 0
      Do a=n,0,-1
         Do b=n-a,0,-1
            c=n-a-b
            ip = ip + 1
!           Write (*,*) ip, I(a,b,c,ix,iy,iz)
            If (I(a,b,c,ix,iy,iz).ne.0)                                 &
     &      ESIT(ip)=ESIT(ip)+DBLE(I(a,b,c,ix,iy,iz))*temp
         End Do
      End Do
!
      Return
      End Subroutine ContEI
