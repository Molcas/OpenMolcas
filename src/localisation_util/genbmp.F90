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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: n, m, Lunit, nStp
real(kind=wp), intent(in) :: X(n,m), StpSiz
character, intent(in) :: Color
integer(kind=iwp) :: i, iBin, iUpLim, j, nBin, nBin_Def, nCh
real(kind=wp) :: absX, Step, Step_Def
character :: Color_Def, myColor
integer(kind=iwp), allocatable :: BMp(:)
real(kind=wp), allocatable :: Bin(:)
character(len=*), parameter :: SecNam = 'GenBMp'
integer(kind=iwp), external :: iRnge

! Defaults and limits.
! --------------------

nBin_Def = 5
Step_Def = 1.0e-2_wp
Color_Def = 'R'
iUpLim = 999999

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

call mma_allocate(Bin,nBin,label='Bins')
call mma_allocate(BMp,nBin,label='iBMp')

Bin(1) = One
do iBin=2,nBin-1
  Bin(iBin) = Bin(iBin-1)*Step
end do
Bin(nBin) = -One

BMp(nBin) = 255
do iBin=nBin-1,1,-1
  BMp(iBin) = BMp(iBin+1)-nCh
end do

! Generate bitmap file.
! Note the special loop structure.
! --------------------------------

write(Lunit,'(2(1X,I6))') m,n ! #col,#row
if (myColor == 'R') then ! red
  do i=n,1,-1
    do j=1,m
      absX = abs(X(i,j))
      iBin = iRnge(absX,Bin,nBin)
      if (BMp(iBin) == 255) then
        write(Lunit,'(4(1X,I3))') 255,255,255,0
      else
        write(Lunit,'(4(1X,I3))') BMp(iBin),0,0,0
      end if
    end do
  end do
else if (myColor == 'G') then ! green
  do i=n,1,-1
    do j=1,m
      absX = abs(X(i,j))
      iBin = iRnge(absX,Bin,nBin)
      if (BMp(iBin) == 0) then
        write(Lunit,'(4(1X,I3))') 255,255,255,0
      else
        write(Lunit,'(4(1X,I3))') 0,BMp(iBin),0,0
      end if
    end do
  end do
else if (myColor == 'B') then ! blue
  do i=n,1,-1
    do j=1,m
      absX = abs(X(i,j))
      iBin = iRnge(absX,Bin,nBin)
      if (BMp(iBin) == 0) then
        write(Lunit,'(4(1X,I3))') 255,255,255,0
      else
        write(Lunit,'(4(1X,I3))') 0,0,BMp(iBin),0
      end if
    end do
  end do
else
  call SysAbendMsg(SecNam,'Logical error!','(Should never happen)')
end if

! De-allocations.
! ---------------

call mma_deallocate(Bin)
call mma_deallocate(BMp)

end subroutine GenBMp
