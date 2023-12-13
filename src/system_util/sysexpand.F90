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
! Copyright (C) 2001, Valera Veryazov                                  *
!***********************************************************************

subroutine SysExpand(strin,strout,iRet)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     general purpose routine for expanding frequently used messages   *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     V.Veryazov University of Lund, 2001                              *
!                                                                      *
!***********************************************************************
!   character s*256
!   call sysexpand('MSG: open ',s,i)  - correct call
!   call sysexpand('MSG: opens ',s,i) - uncorrect, but could be fixed
!   call sysexpand('MSG: uNit ',s,i)  - mixed case
!   call sysexpand('MSG: blah ',s,i)  - will print BLAH  i=real len of s
!   call sysexpand('blah-blah ',s,i)  - will print string as is and return i=0
!***********************************************************************

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: strin
character(len=*), intent(out) :: strout
integer(kind=iwp), intent(out) :: iRet
integer(kind=iwp), save :: ifset = 0, itab(0:255)
integer(kind=iwp) :: i, ii, ik, ip, j
character(len=512) :: sstrin
character :: c
character(len=*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ ', &
                               lw = 'abcdefghijklmnopqrstuvwxyz ', &
                               printable = '1234567890-=~!@#$%^&*()_+<>,.?/\[]":;'
! FORTRAN hash :-)
integer(kind=iwp), parameter :: MAXlabel = 8
integer(kind=iwp), save :: MSGlen(MAXlabel)
character(len=*), parameter :: MSGlabel(MAXlabel) = [character(len=13) :: &
                                                     'OPEN', &
                                                     'CLOSE', &
                                                     'UNIT', &
                                                     'DELETE', &
                                                     'SEEK', &
                                                     'INVALIDOPTION', &
                                                     'USED', &
                                                     'NOTOPENED' &
                                                    ], &
  MSGtext(MAXlabel) = [character(len=128) :: &
                       'Premature abort while opening file',                            & !OPEN
                       'Premature abort while closing the file',                        & !CLOSE
                       'Invalid unit number (Lu<=0 or Lu>99)',                          & !UNIT
                       'Premature abort while removing the file',                       & !DELETE
                       'Premature abort while seeking the file',                        & !SEEK
                       'An invalid option or combination of options has been supplied', & !INVALIDOPTION
                       'Invalid unit number. The file is already opened',               & !USED
                       'File is not Opened'                                             & !NOTOPENED
                      ]

! preset of saved data
! this code uses a part of upcase routine
if (ifset == 0) then
  ifset = 1
  do i=0,255
    itab(i) = -1
  end do
  do ii=1,26
    i = ichar(up(ii:ii))
    j = ichar(lw(ii:ii))
    itab(j) = i
    itab(i) = i
  end do
  do i=1,MAXlabel
    do j=128,1,-1
      if (MSGtext(i)(j:j) /= ' ') exit
    end do
    MSGlen(i) = j
  end do
end if
! fixation of bug in G77
sstrin = strin
! no action if it's an ordinary message
if (sstrin(1:4) /= 'MSG:') then
  do ii=1,len(sstrin)
    c = sstrin(ii:ii)
      ip = index(up,c)+index(lw,c)+index(printable,c)
      !write(u6,*) 'ii',ii
      if (ip == 0) sstrin(ii:ii) = ' '
  end do
  iRet = 0
  return
end if
! uppercase with removing noncharacters.
ik = 0
do ii=5,len(sstrin)
  i = ichar(sstrin(ii:ii))
  j = itab(i)
  if (j >= 0) then
    ik = ik+1
    sstrin(ik:ik) = char(j)
  end if
end do
strout = sstrin(1:ik)
iRet = ik
do ii=1,MAXlabel
  if (sstrin(1:ik) == MSGlabel(ii)) then
    strout = MSGtext(ii)(1:MSGlen(ii))
    iRet = MSGlen(ii)
    return
  end if
end do
! try again...
do ii=1,MAXlabel
  if (sstrin(1:4) == MSGlabel(ii)(1:4)) then
    strout = MSGtext(ii)(1:MSGlen(ii))
    iRet = MSGlen(ii)
    return
  end if
end do

return

end subroutine SysExpand
