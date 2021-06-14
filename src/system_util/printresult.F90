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

subroutine PrintResult(iUnit,FMT,STR,iCount,STR2,value,iRank)
! routine to print result in the form:
!      write(iUnit,FMT) STR,iCount,STR2,(Value(i),i=1,iRank)
!  or, if iCount=0
!      write(iUnit,FMT) STR,(Value(i),i=1,iRank)

character*(*) STR
character*(*) STR2
character*(*) FMT
character*120 TMP
character*2 Marker
real*8 value(iRank)
#include "icolorize.fh"

if (icolorize == 1) then
  Marker = '::'
  if (iCount == 0) then
    write(TMP,FMT) STR,(value(i),i=1,iRank)
  else
    write(TMP,FMT) STR,iCount,STR2,(value(i),i=1,iRank)
  end if
  init = 1
  if (TMP(1:3) == '   ') init = 3
  write(iUnit,'(a,a)') Marker,TMP(init:mylen(TMP))
else
  if (iCount == 0) then
    write(iUnit,FMT) STR,(value(i),i=1,iRank)
  else
    write(iUnit,FMT) STR,iCount,STR2,(value(i),i=1,iRank)
  end if
end if

return

end subroutine PrintResult
