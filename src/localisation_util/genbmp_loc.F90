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

implicit none
integer nRow, nCol
real*8 X(nRow,nCol)
character*(*) FilNam
character*1 Color
character*10 SecNam
parameter(SecNam='GenBMp_Loc')
integer isFreeUnit
external isFreeUnit
character*80 Txt
integer Lunit, irc, nStp
real*8 StpSiz

Lunit = isFreeUnit(11)
call Molcas_Open(Lunit,FilNam)
irc = 0
nStp = -1
StpSiz = -1.0d0
call GenBMp(irc,X,nRow,nCol,Lunit,nStp,StpSiz,Color)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'GenBMp returned',irc
  call SysAbendMsg(SecNam,'GenBMp failed!',Txt)
end if
close(Lunit,Status='Keep')

end subroutine GenBMp_Loc
