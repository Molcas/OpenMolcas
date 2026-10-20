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

use Definitions, only: iwp, u6

implicit none
private

integer(kind=iwp) :: iPrGlb = 0, iPrLoc(7)

public :: iPrGlb, iPrLoc, set_print_level

contains

! Determines the global print level and local print levels
! Note, that it is impossible for iPrLoc to ever differ from
! iPrGlb, and therefore iPrLoc should just be removed from
! this module all together
subroutine set_print_level()

  use PrintLevel, only: DEBUG, SILENT, USUAL

  integer(kind=iwp) :: i
  integer(kind=iwp), external :: iPrintLevel
  logical(kind=iwp), external :: reduce_prt

  iPrGlb = iPrintLevel(-1)
  ! If inside an optimization loop, set down the print level
  ! unless we *really* want a lot of output
  if (reduce_prt()) iPrGlb = max(iPrGlb-USUAL,SILENT)

  iPrLoc(:) = iPrGlb

  if (iPrGlb >= DEBUG) then
    write(u6,*) ' set_print_level: Print levels have been set to'
    write(u6,*) '  Global print level iPrGlb=',iPrGlb
    write(u6,*) '  Individual sections print levels, iPrLoc:'
    write(u6,'(1x,7I5)') (iPrLoc(i),i=1,7)
  end if

end subroutine set_print_level

end module mcpdft_output
