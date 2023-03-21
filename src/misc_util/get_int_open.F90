!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine Get_Int_Open(iSymp,iSymq,iSymr,iSyms)

use GetInt_mod, only: LuCVec, pq1
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iSymp, iSymq, iSymr, iSyms
character(len=6) :: Fname
character(len=*), parameter :: BaseNm = 'CHFV'

! Open files.
LuCVec(1) = 7
write(Fname,'(A4,I1,I1)') BaseNm,iSymp,iSymq
call DANAME_MF_WA(LuCVec(1),Fname)
if (iSymp /= iSymr) then
  LuCVec(2) = 7
  write(Fname,'(A4,I1,I1)') BaseNm,iSymr,iSyms
  call DANAME_MF_WA(LuCVec(2),Fname)
else
  LuCVec(2) = -1
end if

pq1=1

end subroutine Get_Int_Open
