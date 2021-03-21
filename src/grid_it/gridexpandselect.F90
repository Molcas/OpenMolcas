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

subroutine gridExpandSelect(SelectStr)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use grid_it_globals, only: iReq, MAXGRID, nReq
use Definitions, only: iwp, u6

implicit none
character(len=120), intent(in) :: SelectStr
integer(kind=iwp) :: i, ibr, ibrm, iend, ifirst, iibeg, iiend, ilen, istart, istatus, isymm
character(len=120) :: tmp

ifirst = 0
ilen = 0
tmp = ' '
do i=1,len(SelectStr)
  if ((ifirst == 0) .and. (SelectStr(i:i) /= ' ')) then
    ilen = ilen+1
    tmp(ilen:ilen) = SelectStr(i:i)
    ifirst = 1
  else if ((ifirst == 1) .and. (SelectStr(i:i) /= ' ')) then
    ilen = ilen+1
    tmp(ilen:ilen) = SelectStr(i:i)
  else if ((ifirst == 1) .and. (SelectStr(i:i) == ' ')) then
    ilen = ilen+1
    tmp(ilen:ilen) = SelectStr(i:i)
    ifirst = 0
  end if
end do
if (ilen < 2) then
  write(u6,*) 'SELEct section is incomplete'
  call Quit_OnUserError()
end if
!write(u6,*) 'current',tmp
nReq = 0
do
  istart = 1
  iend = index(tmp(istart:),' ')
  ibr = index(tmp(istart:iend),':')
  if (ibr == 0) then
    write(u6,*) 'Wrong format in SELEct section'
    write(u6,*) 'Expecting : sign in >',tmp(istart:iend),'<'
    call Quit_OnUserError()
  end if
  !write(u6,*) 'v01 >',tmp(istart:istart+ibr-2),'<'
  read(tmp(istart:istart+ibr-2),*,iostat=istatus) isymm
  if (istatus /= 0) then
    write(u6,*) 'Error in analyzing SELECT section'
    call Quit_OnUserError()
  end if
  ibrm = index(tmp(istart+ibr+1:iend),'-')
  if (ibrm == 0) then
    ! the only number
    !write(u6,*) 'v02 >',tmp(istart+ibr:iend),'<'
    read(tmp(istart+ibr:iend),*,iostat=istatus) iibeg
    if (istatus /= 0) then
      write(u6,*) 'Error in analyzing SELECT section'
      call Quit_OnUserError()
    end if
    iReq(2*nReq+1) = isymm
    iReq(2*nReq+2) = iibeg
    nReq = nReq+1
    if (nReq > MAXGRID) then
      write(u6,*) 'Too many Grids requested'
      call Quit_OnUserError()
    end if
  else
    !write(u6,*) 'v03 >',tmp(istart+ibr:istart+ibr+ibrm-1),'<'
    read(tmp(istart+ibr:istart+ibr+ibrm-1),*,iostat=istatus) iibeg
    if (istatus /= 0) then
      write(u6,*) 'Error in analyzing SELECT section'
      call Quit_OnUserError()
    end if
    !write(u6,*) 'v04 >',tmp(istart+ibr+ibrm+1:iend),'<'
    read(tmp(istart+ibr+ibrm+1:iend),*,iostat=istatus) iiend
    if (istatus /= 0) then
      write(u6,*) 'Error in analyzing SELECT section'
      call Quit_OnUserError()
    end if
    if (iiend < iibeg) then
      write(u6,*) 'Wrong data in SELEct section'
      call Quit_OnUserError()
    end if
    do i=iibeg,iiend
      iReq(2*nReq+1) = isymm
      iReq(2*nReq+2) = i
      nReq = nReq+1
      if (nReq > MAXGRID) then
        write(u6,*) 'Too many Grids requested'
        call Quit_OnUserError()
      end if
    end do
  end if
  !write(u6,*) 'current',tmp(iend+1:)
  tmp = tmp(iend+1:)
  if (tmp == ' ') exit
end do

return

end subroutine GridExpandSelect
