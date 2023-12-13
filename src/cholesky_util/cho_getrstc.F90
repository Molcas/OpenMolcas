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

subroutine CHO_GETRSTC()
!
! Purpose: read and check decomposition restart info.

use Cholesky, only: LuPri, ModRst
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IFAIL
character(len=*), parameter :: SECNAM = 'CHO_GETRSTC'

! Read restart file, populating restart common block.
! ---------------------------------------------------

IFAIL = 0
call CHO_RDRSTC(IFAIL)
if (IFAIL /= 0) then
  write(LUPRI,'(A,A)') SECNAM,': error reading decomposition restart file.'
  write(LUPRI,'(A,A,I10)') SECNAM,': return code from reading routine:',IFAIL
  call CHO_QUIT('Error reading decomposition restart file',104)
end if

! Check system info .
! -------------------

IFAIL = 0
call CHO_RSTMOL(IFAIL)
if (IFAIL /= 0) then
  write(LUPRI,'(A,A)') SECNAM,': decomposition restart failure.'
  call CHO_QUIT('Decomposition restart failure in '//SECNAM,105)
end if

! Check decomposition configuration.
! ----------------------------------

IFAIL = 0
call CHO_RSTCNF(IFAIL)
if (IFAIL /= 0) then
  write(LUPRI,'(A,A,I6,A)') SECNAM,':',IFAIL,' configuration discrepancies detected.'
  if (MODRST == -1) then
    write(LUPRI,'(A)') 'Recovery: using configuration from restart file.'
    call CHO_RESETCNF()
  else if (MODRST == 0) then
    write(LUPRI,'(A)') 'Recovery: none, program stops.'
    call CHO_QUIT('Restart configuration error in '//SECNAM,105)
  else if (MODRST == 1) then
    write(LUPRI,'(A)') 'Recovery: using input configuration.'
  else
    write(LUPRI,'(A,A,I6,A)') SECNAM,': restart model,',MODRST,', not recognized.'
    call CHO_QUIT('Error in '//SECNAM,103)
  end if
end if

end subroutine CHO_GETRSTC
