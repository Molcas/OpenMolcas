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

implicit real*8(a-h,o-z)
dimension icheckxy(*), icheckz(*), interxyz(16,*)
#include "para.fh"
#include "Molcas.fh"

!bs the following M values are the ones from the cartesian
!bs linear combinations. interxyz gives the sign sequence
!bs for interacting spherical functions, starting with
!bs type 1 (++++) and ending with type 16 (-++-)
ilauf = 1
do M4=0,Lmax
  do M3=0,Lmax
    do M2=0,Lmax
      do M1=0,Lmax
        irun = 0
        if (icheckxy(ilauf)+icheckz(ilauf) > 0) then
          if (abs(m1+m2-m3-m4) <= 1) then
            irun = irun+1
            interxyz(irun,ilauf) = 1    ! + + + +
            if ((m1 > 0) .and. (m2 > 0) .and. (m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 2  ! - - - -
            end if
          end if
          if (abs(m1+m2-m3+m4) <= 1) then
            if (m4 > 0) then
              irun = irun+1
              interxyz(irun,ilauf) = 3  ! + + + -
            end if
            if ((m1 > 0) .and. (m2 > 0) .and. (m3 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 4  ! - - - +
            end if
          end if
          if (abs(m1+m2+m3-m4) <= 1) then
            if (m3 > 0) then
              irun = irun+1
              interxyz(irun,ilauf) = 5  ! + + - +
            end if
            if ((m1 > 0) .and. (m2 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 6  ! - - + -
            end if
          end if
          if (abs(m1-m2-m3-m4) <= 1) then
            if (m2 > 0) then
              irun = irun+1
              interxyz(irun,ilauf) = 7  ! + - + +
            end if
            if ((m1 > 0) .and. (m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 8  ! - + - -
            end if
          end if
          if (abs(-m1+m2-m3-m4) <= 1) then
            if (m1 > 0) then
              irun = irun+1
              interxyz(irun,ilauf) = 9  ! - + + +
            end if
            if ((m2 > 0) .and. (m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 10 ! + - - -
            end if
          end if
          if (abs(m1+m2+m3+m4) <= 1) then
            if ((m3 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 11 ! + + - -
            end if
            if ((m1 > 0) .and. (m2 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 12 ! - - + +
            end if
          end if
          if (abs(m1-m2-m3+m4) <= 1) then
            if ((m2 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 13 ! + - + -
            end if
            if ((m1 > 0) .and. (m3 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 14 ! - + - +
            end if
          end if
          if (abs(m1-m2+m3-m4) <= 1) then
            if ((m2 > 0) .and. (m3 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 15 ! + - - +
            end if
            if ((m1 > 0) .and. (m4 > 0)) then
              irun = irun+1
              interxyz(irun,ilauf) = 16 ! - + + -
            end if
          end if
        end if
        ilauf = ilauf+1
      end do
    end do
  end do
end do

return

end subroutine genprexyz15a
