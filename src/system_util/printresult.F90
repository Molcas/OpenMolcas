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

subroutine PrintResult(iUnit,FRMT,STR,iCount,STR2,Val,iRank)
! routine to print result in the form:
!      write(iUnit,FRMT) STR,iCount,STR2,(Val(i),i=1,iRank)
!  or, if iCount=0
!      write(iUnit,FRMT) STR,(Val(i),i=1,iRank)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iUnit, iCount, iRank
character(len=*), intent(in) :: FRMT, STR, STR2
real(kind=wp), intent(in) :: Val(iRank)
character(len=120) :: TMP
character(len=2), parameter :: Marker = '::'
#include "print.fh"

if (icolorize == 1) then
  if (iCount == 0) then
    write(TMP,FRMT) STR,Val(1:iRank)
  else
    write(TMP,FRMT) STR,iCount,STR2,Val(1:iRank)
  end if
  if (TMP(1:3) == '   ') TMP = TMP(3:)
  write(iUnit,'(a,a)') Marker,trim(TMP)
else
  if (iCount == 0) then
    write(iUnit,FRMT) STR,Val(1:iRank)
  else
    write(iUnit,FRMT) STR,iCount,STR2,Val(1:iRank)
  end if
end if

return

end subroutine PrintResult
