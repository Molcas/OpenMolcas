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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Banner(Lines,nLines,nWidth)
!***********************************************************************
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!***********************************************************************

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nLines, nWidth
character(len=*), intent(in) :: Lines(nLines)
integer(kind=iwp), parameter :: MxWdth = 132
character(len=MxWdth-2) :: Line
character(len=72) :: frmt
integer(kind=iwp) :: i, iEnd, iFrst, j, jEnd, jFrst, Length, mWidth, nChar, nSplit

mWidth = nWidth
nChar = len(Lines(1))
if (nChar+2 > mWidth) mWidth = nChar+2
mWidth = min(MxWdth-2,mWidth)
write(frmt,'(A,i3,A)') '(1X,A',mWidth,')'
Line(1:mWidth) = repeat('*',mWidth)
write(u6,frmt) Line
Line(2:mWidth-1) = repeat(' ',mWidth-2)
write(u6,frmt) Line
do i=1,nLines
  do j=1,nChar
    if (Lines(i)(j:j) /= ' ') exit
  end do
  iFrst = j
  do j=nChar,iFrst,-1
    if (Lines(i)(j:j) /= ' ') exit
  end do
  iEnd = j
  Line(2:mWidth-1) = repeat(' ',mWidth-2)
  Length = iEnd-iFrst+1
  nSplit = (mWidth-2-Length)/2
  jFrst = 1+nSplit+1
  jEnd = jFrst+Length-1
  Line(jFrst:jEnd) = Lines(i)(iFrst:iEnd)
  write(u6,frmt) Line
end do

Line(2:mWidth-1) = repeat(' ',mWidth-2)
write(u6,frmt) Line
Line(2:mWidth-1) = repeat('*',mWidth-2)
write(u6,frmt) Line

return

end subroutine Banner
