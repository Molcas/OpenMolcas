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

subroutine CHO_RSTMOL(NERR)
!
! Purpose: check restart molecular info.

use Cholesky, only: LUPRI, NBAS, NNSHL, NSHELL, NSYM, XNBAS, XNNSHL, XNSHELL, XNSYM
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: NERR
integer(kind=iwp) :: ISYM

NERR = 0

if (XNSYM /= NSYM) then
  write(LUPRI,'(A,I3,A,I3)') 'RESTART ERROR: #irreps from restart file:',XNSYM,' Expected:',NSYM
  NERR = NERR+1
else
  do ISYM=1,NSYM
    if (XNBAS(ISYM) /= NBAS(ISYM)) then
      write(LUPRI,'(A,I2,A,I9,A,I9)') 'RESTART ERROR: #basis functions (sym.',ISYM,') from restart file:',XNBAS(ISYM), &
                                      ' Expected:',NBAS(ISYM)
      NERR = NERR+1
    end if
  end do
end if

if (XNSHELL /= NSHELL) then
  write(LUPRI,'(A,I9,A,I9)') 'RESTART ERROR: #shells from restart file:',XNSHELL,' Expected:',NSHELL
  NERR = NERR+1
end if

if (XNNSHL /= NNSHL) then
  write(LUPRI,'(A,I9,A,I9)') 'RESTART ERROR: #shell pairs from restart file:',XNNSHL,' Expected:',NNSHL
  NERR = NERR+1
end if

end subroutine CHO_RSTMOL
