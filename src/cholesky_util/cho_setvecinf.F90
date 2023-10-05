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

subroutine CHO_SETVECINF(IVEC,ISYM,IAB,IPASS,ILOC)
!
! Purpose: set info for vector IVEC of sym. ISYM.

use Cholesky, only: InfVec, LuPri, MaxVec, nnBstR
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, ISYM, IAB, IPASS, ILOC
character(len=*), parameter :: SECNAM = 'CHO_SETVECINF'

if (IVEC > MAXVEC) then
  write(LUPRI,*) SECNAM,': too many Cholesky vectors!'
  write(LUPRI,*) SECNAM,': symmetry: ',ISYM
  write(LUPRI,*) SECNAM,': max. allowed is ',MAXVEC
  write(LUPRI,*) SECNAM,': please increase max. allowed'
  call CHO_QUIT('Too many Cholesky vectors in '//SECNAM,104)
else if (IVEC == MAXVEC) then ! no set next addr.
  INFVEC(IVEC,1,ISYM) = IAB   ! diag. index red. set 1
  INFVEC(IVEC,2,ISYM) = IPASS ! global red. set
else
  INFVEC(IVEC,1,ISYM) = IAB   ! diag. index red. set 1
  INFVEC(IVEC,2,ISYM) = IPASS ! global red. set
  INFVEC(IVEC+1,4,ISYM) = INFVEC(IVEC,4,ISYM)+NNBSTR(ISYM,ILOC) ! next addr.
end if

end subroutine CHO_SETVECINF
