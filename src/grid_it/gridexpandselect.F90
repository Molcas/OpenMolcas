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

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "grid.fh"
character SelectStr*120, tmp*120

ifirst = 0
ilen = 0
tmp = ' '
do i=1,120
  if (ifirst == 0 .and. SelectStr(i:i) /= ' ') then
    ilen = ilen+1
    tmp(ilen:ilen) = SelectStr(i:i)
    ifirst = 1
    goto 10
  end if
  if (ifirst == 1 .and. SelectStr(i:i) /= ' ') then
    ilen = ilen+1
    tmp(ilen:ilen) = SelectStr(i:i)
  end if
  if (ifirst == 1 .and. SelectStr(i:i) == ' ') then
    ilen = ilen+1
    tmp(ilen:ilen) = SelectStr(i:i)
    ifirst = 0
  end if
10 continue
end do
if (ilen < 2) then
  write(6,*) 'SELEct section is incomplete'
  call Quit_OnUserError()
end if
!write(6,*) 'current',tmp
nReq = 0
1 istart = 1
iend = index(tmp(istart:),' ')
ibr = index(tmp(istart:iend),':')
if (ibr == 0) then
  write(6,*) 'Wrong format in SELEct section'
  write(6,*) 'Expecting : sign in >',tmp(istart:iend),'<'
  call Quit_OnUserError()
end if
!write(6,*) 'v01 >',tmp(istart:istart+ibr-2),'<'
read(tmp(istart:istart+ibr-2),*,err=20,end=20) isymm
ibrm = index(tmp(istart+ibr+1:iend),'-')
if (ibrm == 0) then
  ! the only number
  !write(6,*) 'v02 >',tmp(istart+ibr:iend),'<'
  read(tmp(istart+ibr:iend),*,err=20,end=20) iibeg
  iReq(2*nReq+1) = isymm
  iReq(2*nReq+2) = iibeg
  nReq = nReq+1
  if (nReq > MAXGRID) goto 30
else
  !write(6,*) 'v03 >',tmp(istart+ibr:istart+ibr+ibrm-1),'<'
  read(tmp(istart+ibr:istart+ibr+ibrm-1),*,err=20,end=20) iibeg
  !write(6,*) 'v04 >',tmp(istart+ibr+ibrm+1:iend),'<'
  read(tmp(istart+ibr+ibrm+1:iend),*,err=20,end=20) iiend
  if (iiend < iibeg) then
    write(6,*) 'Wrong data in SELEct section'
    call Quit_OnUserError()
  end if
  do i=iibeg,iiend
    iReq(2*nReq+1) = isymm
    iReq(2*nReq+2) = i
    nReq = nReq+1
    if (nReq > MAXGRID) goto 30
  end do
end if
!write(6,*) 'current',tmp(iend+1:)
tmp = tmp(iend+1:)
if (tmp /= ' ') goto 1

return

20 write(6,*) 'Error in analyzing SELECT section'
call Quit_OnUserError()
30 write(6,*) 'Too many Grids requested'
call Quit_OnUserError()

end subroutine GridExpandSelect
