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

!#define _DEBUGPRINT_
subroutine ZBASE(NVEC,IVEC,NCLASS)
! Some class division exists with NVEC(ICLASS) members in
! class ICLASS.
!
! Construct array IVEC(ICLASS) giving first element of
! class ICLASS in full addressing

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NCLASS, NVEC(NCLASS)
integer(kind=iwp), intent(out) :: IVEC(NCLASS)
integer(kind=iwp) :: ICLASS

do ICLASS=1,NCLASS
  if (ICLASS == 1) then
    IVEC(1) = 1
  else
    IVEC(ICLASS) = IVEC(ICLASS-1)+NVEC(ICLASS-1)
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,'(A)') '  ZBASE : NVEC and IVEC'
write(u6,'(A)') '  ====================='
call IWRTMA(NVEC,1,NCLASS,1,NCLASS)
call IWRTMA(IVEC,1,NCLASS,1,NCLASS)
#endif

return

end subroutine ZBASE
