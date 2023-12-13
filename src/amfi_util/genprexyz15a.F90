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

subroutine genprexyz15a(icheckxy,icheckz,interxyz)

use AMFI_global, only: Lmax
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: icheckxy(0:Lmax,0:Lmax,0:Lmax,0:Lmax), icheckz(0:Lmax,0:Lmax,0:Lmax,0:Lmax)
integer(kind=iwp), intent(inout) :: interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax)
integer(kind=iwp) :: irun, M1, M2, M3, M4

!bs the following M values are the ones from the cartesian
!bs linear combinations. interxyz gives the sign sequence
!bs for interacting spherical functions, starting with
!bs type 1 (++++) and ending with type 16 (-++-)
do M4=0,Lmax
  do M3=0,Lmax
    do M2=0,Lmax
      do M1=0,Lmax
        irun = 0
        if (icheckxy(m1,m2,m3,m4)+icheckz(m1,m2,m3,m4) > 0) then
          if (abs(m1+m2-m3-m4) <= 1) then
            irun = irun+1
            interxyz(irun,m1,m2,m3,m4) = 1    ! + + + +
            if ((m1 > 0) .and. (m2 > 0) .and. (m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 2  ! - - - -
            end if
          end if
          if (abs(m1+m2-m3+m4) <= 1) then
            if (m4 > 0) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 3  ! + + + -
            end if
            if ((m1 > 0) .and. (m2 > 0) .and. (m3 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 4  ! - - - +
            end if
          end if
          if (abs(m1+m2+m3-m4) <= 1) then
            if (m3 > 0) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 5  ! + + - +
            end if
            if ((m1 > 0) .and. (m2 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 6  ! - - + -
            end if
          end if
          if (abs(m1-m2-m3-m4) <= 1) then
            if (m2 > 0) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 7  ! + - + +
            end if
            if ((m1 > 0) .and. (m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 8  ! - + - -
            end if
          end if
          if (abs(-m1+m2-m3-m4) <= 1) then
            if (m1 > 0) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 9  ! - + + +
            end if
            if ((m2 > 0) .and. (m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 10 ! + - - -
            end if
          end if
          if (abs(m1+m2+m3+m4) <= 1) then
            if ((m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 11 ! + + - -
            end if
            if ((m1 > 0) .and. (m2 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 12 ! - - + +
            end if
          end if
          if (abs(m1-m2-m3+m4) <= 1) then
            if ((m2 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 13 ! + - + -
            end if
            if ((m1 > 0) .and. (m3 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 14 ! - + - +
            end if
          end if
          if (abs(m1-m2+m3-m4) <= 1) then
            if ((m2 > 0) .and. (m3 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 15 ! + - - +
            end if
            if ((m1 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,m1,m2,m3,m4) = 16 ! - + + -
            end if
          end if
        end if
      end do
    end do
  end do
end do

return

end subroutine genprexyz15a
