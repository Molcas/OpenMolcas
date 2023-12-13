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

subroutine cp_SpcInt()

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, iq, LuTmp1, LuTmp2, nq, nQQ
logical(kind=iwp) :: Exists
character(len=16) :: filnam
character(len=14) :: qLbl
real(kind=wp), allocatable :: rK(:)

!write(u6,*) ' Copy files'

! Copy file SPCINX to SPCIN1

filnam = 'SPCINX'
call f_Inquire(filnam,Exists)
if (Exists) then
  LuTmp1 = 11
  LuTmp2 = 12
  call molcas_binaryopen_vanilla(luTmp1,filnam)
  !open(luTmp1,File=filnam,Form='unformatted',Status='unknown')
  filnam = 'SPCIN1'
  call molcas_binaryopen_vanilla(luTmp2,filnam)
  !open(luTmp2,File=filnam,Form='unformatted',Status='unknown')
  rewind(LuTmp1)
  rewind(LuTmp2)

  read(LuTmp1) nq,nQQ
  write(LuTmp2) nq,nQQ
  call mma_allocate(rK,nQQ,Label='rK')
  do iq=1,nq
    read(LuTmp1) qLbl,(rK(i),i=1,nQQ)
    write(LuTmp2) qLbl,(rK(i),i=1,nQQ)
  end do
  call mma_deallocate(rK)

  close(LuTmp1)
  close(LuTmp2)
end if

return

end subroutine cp_SpcInt
