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

subroutine zasun_zr(i1,length,valn,jn,kn,ln)
! this routine writes one block of 3-indices and appropriate
! values of integrals into an open TEMP-file
!
! i1 - number of pivot index (I)
! length - number of valid integrals in block (I)
! this routine has also stattemp and tmpnam as inputs,
! but they are imported from ccsort_global

use ccsort_global, only: iokey, lrectemp, lunpublic, mbas, nrectemp, nsize, stattemp, tmpnam
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: i1, length, jn(nsize,mbas), kn(nsize,mbas), ln(nsize,mbas)
real(kind=wp), intent(_IN_) :: valn(nsize,mbas)
integer(kind=iwp) :: f_iostat !, iRec
logical(kind=iwp) :: is_error
integer(kind=iwp), allocatable :: jkl(:)
integer(kind=iwp), parameter :: constj = 1024*1024, constk = 1024

! pack indices

call mma_allocate(jkl,length,label='jkl')
jkl(1:length) = ln(1:length,i1)+constj*jn(1:length,i1)+constk*kn(1:length,i1)

! open corresponding TEMP file in corresponding form

if (iokey == 1) then

  ! Fortran IO

  if (stattemp(i1) == 0) then
    ! file will be opened first time, it must be opened
    ! with the pointer at then first position
    call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
    !open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown')
    stattemp(i1) = 1

  else
    ! file was already used in expansion of this block, it must
    ! be opened with the pointer at the end of the file
!#   ifdef _DECAXP_
    call molcas_open_ext2(lunpublic,tmpnam(i1),'append','unformatted',f_iostat,.false.,1,'unknown',is_error)
    !vv open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown',access='append')
!#   else
!    call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
!    open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown')
!    do iRec=1,nrectemp(i1)
!      Read(lunpublic) m2
!    end do
!#   endif

  end if

  write(lunpublic) valn(1:length,i1),jkl(1:length)
  close(lunpublic)

else

  ! MOLCAS IO

  call daname(lunpublic,tmpnam(i1))
  call ddafile(lunpublic,1,valn(:,i1),length,stattemp(i1))
  call idafile(lunpublic,1,jkl,length,stattemp(i1))
  call daclos(lunpublic)

end if

call mma_deallocate(jkl)

nrectemp(i1) = nrectemp(i1)+1
lrectemp(i1) = length

return

end subroutine zasun_zr
