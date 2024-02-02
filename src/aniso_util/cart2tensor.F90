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

subroutine cart2tensor(Jc,Jt)
! this subroutine computed tensorial parameters for Hexch= S1.Jcart.S2

implicit none
integer, parameter :: wp = kind(0.d0)
real(kind=8), intent(in) :: Jc(3,3)
complex(kind=8), intent(out) :: Jt(-1:1,-1:1)
! local variables
complex(kind=8) :: c2, cx2, i
complex(kind=8) :: xx, xy, xz, yx, yy, yz, zx, zy, zz

jt(-1:1,-1:1) = (0.0_wp,0.0_wp)

i = (0.0_wp,1.0_wp)
c2 = (0.5_wp,0.0_wp)
cx2 = cmplx(sqrt(0.5_wp),0.0_wp,wp)

xx = cmplx(jc(1,1),0.0_wp,wp)
xy = cmplx(jc(1,2),0.0_wp,wp)
xz = cmplx(jc(1,3),0.0_wp,wp)
yx = cmplx(jc(2,1),0.0_wp,wp)
yy = cmplx(jc(2,2),0.0_wp,wp)
yz = cmplx(jc(2,3),0.0_wp,wp)
zx = cmplx(jc(3,1),0.0_wp,wp)
zy = cmplx(jc(3,2),0.0_wp,wp)
zz = cmplx(jc(3,3),0.0_wp,wp)

jt(1,1) = c2*(xx-i*xy-i*yx-yy)
jt(-1,-1) = c2*(xx+i*xy+i*yx-yy)
jt(1,-1) = c2*(-xx-i*xy+i*yx-yy)
jt(-1,1) = c2*(-xx+i*xy-i*yx-yy)

jt(1,0) = cx2*(-xz+i*yz)
jt(-1,0) = cx2*(xz+i*yz)

jt(0,1) = cx2*(-zx+i*zy)
jt(0,-1) = cx2*(zx+i*zy)

jt(0,0) = zz

return

end subroutine cart2tensor
