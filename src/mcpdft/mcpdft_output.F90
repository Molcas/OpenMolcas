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

module mcpdft_output
  use definitions, only: wp, iwp

  implicit none
  private

  ! amount of output written
  integer(kind=iwp) :: silent  = 0
  integer(kind=iwp) :: terse   = 1
  integer(kind=iwp) :: usual   = 2
  integer(kind=iwp) :: verbose = 3
  integer(kind=iwp) :: debug   = 4
  Integer(kind=iwp) :: insane  = 5

  integer(kind=iwp), dimension(7) :: iPrLoc
  integer(kind=iwp) :: lf=6, iPrGlb

  public :: silent, terse, usual, verbose, debug, insane, lf, iPrGlb, iPrLoc
  public :: set_print_level

  contains
  subroutine set_print_level(global_print, local_print)
    implicit none

    logical, external :: reduce_prt

    integer(kind=iwp), intent(in) :: global_print
    integer(kind=iwp), dimension(7), intent(in) :: local_print

    integer :: i ! dummy loop variable

    iPrGlb = global_print
    if (iPrGlb == silent) then
      do i=1, 7
        iPrLoc(i) = 0
      end do
    else
      do i=1, 7
        iPrLoc(i) = 0
        if (local_print(i) > 0) then
          iPrLoc(i) = max(iPrGlb, iPrLoc(i))
        end if
      end do
    end if

    ! If inside an optimization loop, set down the print level
    ! unless we *really* want a lot of output
    if (reduce_prt()) then
      iPrGlb = max(iPrGlb - usual, silent)
      do i=1, 7
        iPrLoc(i) = max(iPrLoc(i)-usual, silent)
      end do
    end if

    if (iPrLoc(1) >= debug) then
      write(lf, *) ' set_print_level: Print levels have been set to'
      write(lf, *) '  Global print level iPrGlb=', iPrGlb
      write(lf, *) '  Individual sections print levels, iPrLoc:'
      write(lf, '(1x,7I5)') (iPrLoc(i), i=1,7)
    end if

  end subroutine set_print_level
end module mcpdft_output