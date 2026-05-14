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

subroutine Dipol_Exchange(N1,N2,vec,dist,M1,M2,HDIP)
! this Subroutine computes the dipolar coupling between the two moments

use Constants, only: Zero, Three, cZero, cOne, cm_s, hPlanck, mBohr
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N1, N2
real(kind=wp), intent(in) :: vec(3), dist
complex(kind=wp), intent(in) :: M1(3,N1,N1), M2(3,N2,N2)
complex(kind=wp), intent(out) :: HDIP(N1,N1,N2,N2)
integer(kind=iwp) :: i1, i2, j1, j2, m
complex(kind=wp) :: d3, HL, mb2c, p1, p2a, p2b, vec1(3)
real(kind=wp), parameter :: MB2 = 1.0e23_wp*mBohr**2/(cm_s*hPlanck) ! -- the value of (mu_Bohr*mu_Bohr)/(angstrom^3)  in cm-1
complex(kind=wp), parameter :: threeC = Three*cOne

if ((N1 <= 0) .or. (N2 <= 0)) return
if (dist == Zero) then
  HDIP(:,:,:,:) = cZero
  write(u6,'(A)') 'DIPOL_EXCHANGE::  dist = 0'
  write(u6,'(A)') 'this is not normal. Stop.'
  return
end if
! switch to Complex arithmetic:
d3 = dist*dist*dist*cOne
mb2c = MB2*cOne
vec1(:) = vec(:)*cOne
! calculate the dipolar coupling:
do i1=1,N1
  do j1=1,N1
    do i2=1,N2
      do j2=1,N2
        p2a = cZero
        p2b = cZero
        p1 = cZero
        do m=1,3
          p2a = p2a+M1(m,i1,j1)*vec1(m)
          p2b = p2b+M2(m,i2,j2)*vec1(m)
          p1 = p1+M1(m,i1,j1)*M2(m,i2,j2)
        end do
        HL = p1-threeC*p2a*p2b
        HDIP(i1,j1,i2,j2) = mb2c*HL/d3
      end do
    end do
  end do
end do

return

end subroutine Dipol_Exchange
