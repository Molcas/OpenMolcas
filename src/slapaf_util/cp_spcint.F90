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

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
character*14 qLbl, filnam*16
logical Exist
real*8, allocatable :: rK(:)

!write(6,*) ' Copy files'

! Copy file SPCINX to SPCIN1

filnam = 'SPCINX'
call f_Inquire(filnam,Exist)
if (Exist) then
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
