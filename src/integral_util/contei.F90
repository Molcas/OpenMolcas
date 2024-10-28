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

subroutine ContEI(I,mDeg,ESIT,ix,iy,iz,temp)

use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mDeg, I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg), ix, iy, iz
real(kind=wp), intent(inout) :: ESIT(nTri_Elem1(mDeg))
real(kind=wp), intent(in) :: Temp
integer(kind=iwp) :: a, b, c, ip, n

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

! write out only elements with a>=b>=c. The others are obtained
! by index permutation.
! This restriction has been removed! (Roland Lindh)
!n = mDeg
!do a=n,0,-1
!  do b=n-a,0,-1
!    c = n-a-b
!    write(6,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,Ind(n,a,c)
!    do p=n,0,-1
!      do q=n-p,0,-1
!        r = n-p-q
!        coef = I(a,b,c,p,q,r)
!        if (coef /= 0) write(6,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)') coef,p,q,r,Ind(n,p,r)
!      end do
!    end do
!  end do
!end do

!write(u6,*) ' Temp,ix,iy,iz=',temp,ix,iy,iz
n = mDeg
ip = 0
do a=n,0,-1
  do b=n-a,0,-1
    c = n-a-b
    ip = ip+1
    !write(u6,*) ip, I(a,b,c,ix,iy,iz)
    if (I(a,b,c,ix,iy,iz) /= 0) ESIT(ip) = ESIT(ip)+real(I(a,b,c,ix,iy,iz),kind=wp)*temp
  end do
end do

return

end subroutine ContEI
