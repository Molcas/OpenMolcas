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


! TODO(matthew-hennefarth): Remove iPrLoc from MC-PDFT module
module mcpdft_output
  use definitions, only: iwp

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
  integer(kind=iwp) :: lf=6, iPrGlb = 0

  public :: silent, terse, usual, verbose, debug, insane, lf, iPrGlb, iPrLoc
  public :: set_print_level

  contains
  subroutine set_print_level()
    ! Determines the global print level and local print levels
    ! Note, that it is impossible for iPrLoc to ever differ from
    ! iPrGlb, and therefore iPrLoc should just be removed from
    ! this module all together
    implicit none

    logical, external :: reduce_prt
    integer(kind=iwp), external :: iPrintLevel

    integer :: i ! dummy loop variable

    iPrGlb = iPrintLevel(-1)
    if (reduce_prt()) then
      ! If inside an optimization loop, set down the print level
      ! unless we *really* want a lot of output
      iPrGlb = max(iPrGlb - usual, silent)
    end if

    iPrLoc = iPrGlb

    if (iPrGlb >= debug) then
      write(lf, *) ' set_print_level: Print levels have been set to'
      write(lf, *) '  Global print level iPrGlb=', iPrGlb
      write(lf, *) '  Individual sections print levels, iPrLoc:'
      write(lf, '(1x,7I5)') (iPrLoc(i), i=1,7)
    end if

  end subroutine set_print_level
end module mcpdft_output