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

subroutine tensor2cart_minus(Jt,Jc)
! Hexch = -Jsph.S1.S2

use Constants, only: Half, Onei
use Definitions, only: wp

implicit none
complex(kind=wp), intent(in) :: Jt(-1:1,-1:1)
real(kind=wp), intent(out) :: Jc(3,3)
complex(kind=wp) :: mm, mp, mz, pm, pp, pz, zm, zp, zz

pp = jt(1,1)
pz = jt(1,0)
pm = jt(1,-1)
zp = jt(0,1)
zz = jt(0,0)
zm = jt(0,-1)
mp = jt(-1,1)
mz = jt(-1,0)
mm = jt(-1,-1)

jc(1,1) = Half*real(-mm+mp+pm-pp)
jc(2,2) = Half*real(mm+mp+pm+pp)

jc(1,2) = Half*real(Onei*mm-Onei*mp+Onei*pm-Onei*pp)
jc(2,1) = Half*real(Onei*mm+Onei*mp-Onei*pm-Onei*pp)

jc(1,3) = sqrt(Half)*real(-zm+zp)
jc(3,1) = sqrt(Half)*real(-mz+pz)

jc(2,3) = sqrt(Half)*real(Onei*zm+Onei*zp)
jc(3,2) = sqrt(Half)*real(Onei*mz+Onei*pz)

jc(3,3) = real(-zz)

return

end subroutine tensor2cart_minus
