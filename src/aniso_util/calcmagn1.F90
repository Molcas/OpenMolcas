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

subroutine calcmagn1(N,E,M,T,MT,Z)
!-----------------------------------------------------------------------
! Input data:
!   N =  number of the states to be considered, scalar integer, input
!   E =  Zeeman energy of these states, array(N), real, input
!   M =  momentum of these states ( diagonal), complex array (N,N), input
!        ( only diagonal elements are used )
!   T =  scalar real number denoting temperature in K, real, input
! Output data:
!   MT=  computed momentum at temperature T, scalar real, output
!   Z =  Boltzmann statistical sum, real scalar
!-----------------------------------------------------------------------

use Constants, only: Zero, cm_s, hPlanck, kBoltzmann
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: E(N), T
complex(kind=wp), intent(in) :: M(N,N)
real(kind=wp), intent(out) :: MT, Z
integer(kind=iwp) :: i, im, mp1
real(kind=wp), parameter :: kB = kBoltzmann/(cm_s*hPlanck) ! in cm-1*K-1
!----------------------------------------------------------------------

Z = Zero
MT = Zero
im = 0
mp1 = 0
im = mod(N,11)

if (im /= 0) then
  do i=1,im !N   ! im
    Z = Z+exp(-(E(i)-E(1))/kB/T)
    MT = MT+exp(-(E(i)-E(1))/kB/T)*real(M(I,I))
  end do
  if (N < 11) then
    MT = MT/Z
    return
  end if
end if

mp1 = im+1
do i=mp1,N,11
  Z = Z+exp(-(E(i)-E(1))/kB/T)+exp(-(E(i+1)-E(1))/kB/T)+exp(-(E(i+2)-E(1))/kB/T)+exp(-(E(i+3)-E(1))/kB/T)+ &
      exp(-(E(i+4)-E(1))/kB/T)+exp(-(E(i+5)-E(1))/kB/T)+exp(-(E(i+6)-E(1))/kB/T)+exp(-(E(i+7)-E(1))/kB/T)+ &
      exp(-(E(i+8)-E(1))/kB/T)+exp(-(E(i+9)-E(1))/kB/T)+exp(-(E(i+10)-E(1))/kB/T)

  MT = MT+exp(-(E(i)-E(1))/kB/T)*real(M(i,i))+exp(-(E(i+1)-E(1))/kB/T)*real(M(i+1,i+1))+exp(-(E(i+2)-E(1))/kB/T)*real(M(i+2,i+2))+ &
       exp(-(E(i+3)-E(1))/kB/T)*real(M(i+3,i+3))+exp(-(E(i+4)-E(1))/kB/T)*real(M(i+4,i+4))+ &
       exp(-(E(i+5)-E(1))/kB/T)*real(M(i+5,i+5))+exp(-(E(i+6)-E(1))/kB/T)*real(M(i+6,i+6))+ &
       exp(-(E(i+7)-E(1))/kB/T)*real(M(i+7,i+7))+exp(-(E(i+8)-E(1))/kB/T)*real(M(i+8,i+8))+ &
       exp(-(E(i+9)-E(1))/kB/T)*real(M(i+9,i+9))+exp(-(E(i+10)-E(1))/kB/T)*real(M(i+10,i+10))
end do

MT = MT/Z

return

end subroutine calcmagn1
