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

subroutine GenBMp_Loc(X,nRow,nCol,FilNam,Color)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nRow, nCol
real(kind=wp), intent(in) :: X(nRow,nCol)
character(len=*), intent(in) :: FilNam
character, intent(in) :: Color
integer(kind=iwp) :: irc, Lunit, nStp
real(kind=wp) :: StpSiz
character(len=80) :: Txt
character(len=*), parameter :: SecNam = 'GenBMp_Loc'
integer(kind=iwp), external :: isFreeUnit

Lunit = isFreeUnit(11)
call Molcas_Open(Lunit,FilNam)
irc = 0
nStp = -1
StpSiz = -One
call GenBMp(irc,X,nRow,nCol,Lunit,nStp,StpSiz,Color)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'GenBMp returned',irc
  call SysAbendMsg(SecNam,'GenBMp failed!',Txt)
end if
close(Lunit,Status='Keep')

end subroutine GenBMp_Loc
