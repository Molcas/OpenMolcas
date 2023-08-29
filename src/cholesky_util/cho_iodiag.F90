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

subroutine CHO_IODIAG(DIAG,IOPT)
!
! Purpose: write/read a copy of diagonal to disk (1st reduced set).
!          The file is opened and closed here.
!          IOPT=1: write
!          IOPT=2: read

use Cholesky, only: LuPri, nnBstRT
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: DIAG(*)
integer(kind=iwp), intent(in) :: IOPT
integer(kind=iwp) :: IADR, LENGTH, LUNIT
character(len=*), parameter :: FNAME = 'CHODIAG', SECNAM = 'CHO_IODIAG_1'

LUNIT = 7
call DANAME(LUNIT,FNAME)
if ((IOPT == 1) .or. (IOPT == 2)) then
  LENGTH = NNBSTRT(1)
  IADR = 0
  call DDAFILE(LUNIT,IOPT,DIAG,LENGTH,IADR)
else   ! error
  write(LUPRI,*) SECNAM,': IOPT out of bounds: ',IOPT
  call CHO_QUIT('Error in '//SECNAM,104)
end if
call DACLOS(LUNIT)

end subroutine CHO_IODIAG
