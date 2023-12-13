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

subroutine zasun_pck(i1,length,valn,jn,kn,ln)
! this routine writes one block of 3-indices and appropriate
! values of integrals into an open TEMP-file
!
! i1 - number of pivot index (I)
! length - number of valid integrals in block (I)
! this routine has also jn,kn,ln,valn
! and stattemp and tmpnam as inputs, but they are
! imported from ccsort_global

use ccsort_global, only: iokey, lrectemp, lunpublic, mbas, nrectemp, nsize, stattemp, tmpnam
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, ItoB, RtoB

implicit none
integer(kind=iwp), intent(in) :: i1, length, jn(nsize,mbas), kn(nsize,mbas), ln(nsize,mbas)
real(kind=wp), intent(in) :: valn(nsize,mbas)
integer(kind=iwp) :: ihelp, iRec, m2
real(kind=wp) :: rhelp
character(len=RtoB+ItoB) :: pphelp
character(len=RtoB+ItoB), allocatable :: pp(:)
integer(kind=iwp), parameter :: constj = 1024**2, constk = 1024

! pack indices and integral values

call mma_allocate(pp,length,label='pp')

do m2=1,length
  rhelp = valn(m2,i1)
  ihelp = ln(m2,i1)+constj*jn(m2,i1)+constk*kn(m2,i1)
  pphelp(1:RtoB) = transfer(rhelp,pphelp(1:RtoB))
  pphelp(RtoB+1:) = transfer(ihelp,pphelp(RtoB+1:))
  pp(m2) = pphelp
end do

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
#   ifdef _DECAXP_
    call molcas_open_ext2(lunpublic,tmpnam(i1),'append','unformatted',f_iostat,.false.,1,'unknown',is_error)

    !open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown',access='append')
#   else
    call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
    !open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown')
    do iRec=1,nrectemp(i1)
      read(lunpublic) m2
    end do
#   endif

  end if

  write(lunpublic) pp
  close(lunpublic)

else

  ! MOLCAS IO

  call daname(lunpublic,tmpnam(i1))
  call cdafile(lunpublic,1,pp,(RtoB+ItoB)*length,stattemp(i1))
  call daclos(lunpublic)

end if

call mma_deallocate(pp)

nrectemp(i1) = nrectemp(i1)+1
lrectemp(i1) = length

return

end subroutine zasun_pck
