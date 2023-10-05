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

subroutine CHO_MCA_INT_1_DBG(DIAG,LEVEL)
!
! Purpose: debug seward interface routine CHO_MCA_INT_1.
!
! LEVEL =  1: test diagonal, reduced set 1 (i.e. initial).
!          2: test diagonal, reduced set 2 (i.e. current).
!          3: test symmetry of integral matrix (shell quadruple-based)

use Cholesky, only: LuPri
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: Level
integer(kind=iwp) :: IRED
logical(kind=iwp) :: LOCDIAG, LOCSYM

call CHO_HEAD('Debugging CHO_MCA_INT_1','=',80,LUPRI)
write(LUPRI,'(A,I2)') 'Debug level',LEVEL

if (LEVEL == 1) then
  LOCDIAG = .true.
  LOCSYM = .false.
  IRED = 1
else if (LEVEL == 2) then
  LOCDIAG = .true.
  LOCSYM = .false.
  IRED = 2
else if (LEVEL == 3) then
  LOCDIAG = .false.
  LOCSYM = .true.
else
  LOCDIAG = .false.
  LOCSYM = .false.
  write(LUPRI,'(A)') 'Debug level not recognized --- debug cancelled!'
end if

if (LOCDIAG) call CHO_MCA_INT_1_DBG1(DIAG,IRED)

if (LOCSYM) call CHO_MCA_INT_1_DBG2()

end subroutine CHO_MCA_INT_1_DBG
