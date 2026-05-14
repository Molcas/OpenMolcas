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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine intchk_cvb(iarr,nmax,nread,ifc,a,lflag)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: iarr(*), lflag
integer(kind=iwp), intent(in) :: nmax, ifc
integer(kind=iwp), intent(out) :: nread
character(len=*), intent(in) :: a
integer(kind=iwp) :: i, ifrom, istr, ito(1), lf, ncnt, nr
integer(kind=iwp), parameter :: ncmp = 4, nkeyw = 3
character(len=*), parameter :: keyword(nkeyw) = ['NONE    ','ALL     ','TO      ']

nread = 0
lf = lflag
do
  call fstring_cvb(keyword,nkeyw,istr,ncmp,1)
  if (istr > 0) lf = lflag
  if (istr == 1) then
    ! 'NONE'
    nread = 0
  else if (istr == 2) then
    ! 'ALL'
    if (lflag == -1) then
      nread = nmax
      do i=1,nmax
        iarr(i) = i
      end do
    else
      nread = 0
      lf = 1-lflag
    end if
  else if (istr == 3) then
    ! 'TO'
    if (nread == nmax) then
      write(u6,'(3a)') ' Too many numbers specified in ',a,' keyword!'
      call abend_cvb()
    else if (nread == 0) then
      write(u6,'(3a)') ' No number before ',a,' -- TO keyword!'
      call abend_cvb()
    end if
    call int_cvb(ito,1,nr,ifc)
    if (nr == -1) then
      write(u6,'(3a)') ' No number after ',a,' -- TO keyword!'
      call abend_cvb()
    end if
    ifrom = iarr(nread)
    if (ifrom > ito(1)) then
      write(u6,*) ' From greater than to:',ifrom,ito(1)
      call abend_cvb()
    else if (nread+ito(1)-ifrom > nmax) then
      write(u6,'(3a)') ' Too many numbers specified in ',a,' keyword!'
      call abend_cvb()
    end if
    do i=ifrom+1,ito(1)
      nread = nread+1
      iarr(nread) = i
    end do
  else
    call int_cvb(iarr(1+nread),nmax-nread,nr,ifc)
    if (nread > 0) lf = lflag
    if (nr == -1) then
      write(u6,'(3a)') ' Too many numbers specified in ',a,' keyword!'
      call abend_cvb()
    end if
    nread = nread+nr
  end if
  if ((istr <= 0) .and. (nr <= 0)) exit
end do

if (lflag /= -1) lflag = lf

do i=1,nread
  if ((iarr(i) < 1) .or. (iarr(i) > nmax)) then
    write(u6,'(3a,i5)') ' Illegal ',a,' number read!',iarr(i)
    write(u6,'(a,i3)') ' Must be in the range 1 --',nmax
    call abend_cvb()
  end if
end do

! Sort numbers and ensure uniqueness:
call sorti_cvb(nread,iarr)
ncnt = 1
do i=2,nread
  if (iarr(i) /= iarr(ncnt)) then
    ncnt = ncnt+1
    iarr(ncnt) = iarr(i)
  end if
end do
nread = min(ncnt,nread)

return

end subroutine intchk_cvb
