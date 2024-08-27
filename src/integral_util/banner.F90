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
implicit none
integer nLines, nWidth
integer, parameter :: MxWdth = 132
character(len=*) Lines(nLines)
character(len=MxWdth-2) Line
character(len=72) format
integer mWidth, nChar, i, j, k, iFrst, iEnd, Length, nSplit, jFrst, jEnd

mWidth = nWidth
nChar = len(Lines(1))
if (nChar+2 > mWidth) mWidth = nChar+2
mWidth = min(MxWdth-2,mWidth)
write(format,'(A,i3,A)') '(1X,A',mWidth,')'
do i=1,mWidth
  Line(i:i) = '*'
end do
write(6,format) Line
do i=2,mWidth-1
  Line(i:i) = ' '
end do
write(6,format) Line
do i=1,nLines
  do j=1,nChar
    if (Lines(i)(j:j) /= ' ') Go To 21
  end do
21 continue
  iFrst = j
  do j=nChar,iFrst,-1
    if (Lines(i)(j:j) /= ' ') Go To 31
  end do
31 continue
  iEnd = j
  do k=2,mWidth-1
    Line(k:k) = ' '
  end do
  Length = iEnd-iFrst+1
  nSplit = (mWidth-2-Length)/2
  jFrst = 1+nSplit+1
  jEnd = jFrst+Length-1
  Line(jFrst:jEnd) = Lines(i)(iFrst:iEnd)
  write(6,format) Line
end do

do k=2,mWidth-1
  Line(k:k) = ' '
end do
write(6,format) Line
do k=2,mWidth-1
  Line(k:k) = '*'
end do
write(6,format) Line

return

end subroutine Banner
