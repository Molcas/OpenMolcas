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

subroutine GenBMp(irc,X,n,m,Lunit,nStp,StpSiz,Color)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: irc, n, m, Lunit, nStp
real(kind=wp) :: X(n,m), StpSiz
character :: Color
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i0, i255, iBin, ipBin, ipBMp, iUpLim, j, nBin, nBin_Def, nCh
real(kind=wp) :: absX, Step, Step_Def
character :: Color_Def, myColor
character(len=*), parameter :: SecNam = 'GenBMp'
integer(kind=iwp), external :: iRnge

! Defaults and limits.
! --------------------

nBin_Def = 5
Step_Def = 1.0e-2_wp
Color_Def = 'R'
iUpLim = 999999
i0 = 0
i255 = 255

! Check input.
! ------------

irc = 0
if ((n < 1) .or. (m < 1)) then
  return
end if
if ((n > iUpLim) .or. (m > iUpLim)) then
  irc = 1
  return
end if
if (Lunit < 1) then
  irc = 2
  return
end if

if ((nStp < 2) .or. (nStp > 256)) then
  nBin = nBin_Def
else
  nBin = nStp
end if
if (StpSiz <= Zero) then
  Step = Step_Def
else
  Step = StpSiz
end if
myColor = Color
call UpCase(myColor)
if ((myColor /= 'R') .and. (myColor /= 'G') .and. (myColor /= 'B')) then
  myColor = Color_Def
end if

! Set up bins.
! ------------

nCh = 255/(nBin-1)

call GetMem('Bins','Allo','Real',ipBin,nBin)
call GetMem('iBMp','Allo','Inte',ipBMp,nBin)

Work(ipBin) = One
do iBin=2,nBin-1
  Work(ipBin-1+iBin) = Work(ipBin+iBin-2)*Step
end do
Work(ipBin-1+nBin) = -One

iWork(ipBMp-1+nBin) = 255
do iBin=nBin-1,1,-1
  iWork(ipBMp-1+iBin) = iWork(ipBMp+iBin)-nCh
end do

! Generate bitmap file.
! Note the special loop structure.
! --------------------------------

write(Lunit,'(2(1X,I6))') m,n ! #col,#row
if (myColor == 'R') then ! red
  do i=n,1,-1
    do j=1,m
      absX = abs(X(i,j))
      iBin = iRnge(absX,Work(ipBin),nBin)
      if (iWork(ipBmp+iBin-1) == 255) then
        write(Lunit,'(4(1X,I3))') i255,i255,i255,i0
      else
        write(Lunit,'(4(1X,I3))') iWork(ipBmp+iBin-1),i0,i0,i0
      end if
    end do
  end do
else if (myColor == 'G') then ! green
  do i=n,1,-1
    do j=1,m
      absX = abs(X(i,j))
      iBin = iRnge(absX,Work(ipBin),nBin)
      if (iWork(ipBmp+iBin-1) == 0) then
        write(Lunit,'(4(1X,I3))') i255,i255,i255,i0
      else
        write(Lunit,'(4(1X,I3))') i0,iWork(ipBmp+iBin-1),i0,i0
      end if
    end do
  end do
else if (myColor == 'B') then ! blue
  do i=n,1,-1
    do j=1,m
      absX = abs(X(i,j))
      iBin = iRnge(absX,Work(ipBin),nBin)
      if (iWork(ipBmp+iBin-1) == 0) then
        write(Lunit,'(4(1X,I3))') i255,i255,i255,i0
      else
        write(Lunit,'(4(1X,I3))') i0,i0,iWork(ipBmp+iBin-1),i0
      end if
    end do
  end do
else
  call SysAbendMsg(SecNam,'Logical error!','(Should never happen)')
end if

! De-allocations.
! ---------------

call GetMem('iBMp','Free','Inte',ipBMp,nBin)
call GetMem('Bins','Free','Real',ipBin,nBin)

end subroutine GenBMp
