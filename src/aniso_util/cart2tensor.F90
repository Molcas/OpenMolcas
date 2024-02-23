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

use Constants, only: Half, cOne, Onei
use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: Jc(3,3)
complex(kind=wp), intent(out) :: Jt(-1:1,-1:1)
complex(kind=wp) :: xx, xy, xz, yx, yy, yz, zx, zy, zz
complex(kind=wp), parameter :: c2 = Half*cOne, cx2 = sqrt(Half)*cOne

xx = jc(1,1)*cOne
xy = jc(1,2)*cOne
xz = jc(1,3)*cOne
yx = jc(2,1)*cOne
yy = jc(2,2)*cOne
yz = jc(2,3)*cOne
zx = jc(3,1)*cOne
zy = jc(3,2)*cOne
zz = jc(3,3)*cOne

jt(1,1) = c2*(xx-Onei*xy-Onei*yx-yy)
jt(-1,-1) = c2*(xx+Onei*xy+Onei*yx-yy)
jt(1,-1) = c2*(-xx-Onei*xy+Onei*yx-yy)
jt(-1,1) = c2*(-xx+Onei*xy-Onei*yx-yy)

jt(1,0) = cx2*(-xz+Onei*yz)
jt(-1,0) = cx2*(xz+Onei*yz)

jt(0,1) = cx2*(-zx+Onei*zy)
jt(0,-1) = cx2*(zx+Onei*zy)

jt(0,0) = zz

return

end subroutine cart2tensor
