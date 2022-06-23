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

subroutine kinemat(ndim,evtkin,type1,type2,Energy)
!bs this routine generates the kinematic A-factors=sqrt((E+mc^2)/(2E))
!bs (type1) and   c*A/(E+mc^2) (type2)
!bs The c in the second kinematic factor comes from Jan Almloef and
!bs Odd Gropen in Rev in Comp.Chem. 8(1996)

use Constants, only: Zero, One, Half, speed => c_in_au
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndim
real(kind=wp), intent(in) :: evtkin(ndim)
real(kind=wp), intent(out) :: type1(ndim), type2(ndim), Energy(ndim)
integer(kind=iwp) :: irun
real(kind=wp), parameter :: fine = One/speed, speed2 = speed**2, speed4 = speed2**2

! E= sqrt(p**2 c**2 + m**2 c**4)
! p**2= 2*m*TKIN
! with m = 1
do irun=1,ndim
  if (evtkin(irun) < Zero) call SysAbendMsg('kinemat','strange kinetic energy ',' ')
  Energy(irun) = sqrt((evtkin(irun)+evtkin(irun))*speed2+speed4)
end do
! sqrt((E+mc^2)/(2E)):
type1(:) = sqrt(Half*(One+speed2/Energy(:)))
! c*A/(E+mc^2)
type2(:) = speed*type1(:)/(Energy(:)+speed2)

return

end subroutine kinemat
