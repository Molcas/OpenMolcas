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

subroutine BasisReader(LuWr,nBase,iglobal,nxbas,xb_label,xb_bas,iErr)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuWr, iglobal
integer(kind=iwp), intent(out) :: nBase, iErr
integer(kind=iwp), intent(inout) :: nxbas
character(len=*), intent(in) :: xb_label(*), xb_bas(*)
#include "g_zmatconv.fh"
integer(kind=iwp) :: i, icount
logical(kind=iwp) :: Found
character(len=48) :: Line, TMP
character(len=2) :: SimbA

iErr = 0
nBase = 0
icount = 1

10 continue
!read(iZMUnit,'(A)',err=9904,end=9999) Line
if (iglobal == 0) then
  Line = trim(xb_label(icount))//'.'//trim(xb_bas(icount))
else
  Line = 'FF.'//trim(xb_bas(icount))
  nxbas = 1
end if
!write(u6,*) '>>',Line
Found = .false.
SimbA = Line(1:2)
if (SimbA(2:2) == '.') SimbA(2:2) = ' '
do i=1,Num_Elem
  if ((SimbA(1:1) <= 'z') .and. (SimbA(1:1) >= 'a')) SimbA(1:1) = char(ichar(SimbA(1:1))-32)
  if ((SimbA(2:2) <= 'Z') .and. (SimbA(2:2) >= 'A')) SimbA(2:2) = char(ichar(SimbA(2:2))+32)

  if (SimbA == adjustl(PTab(i)) .and. iglobal == 0) then
    Base(i) = Line
    BasAva(i) = .true.
    Found = .true.
    nBase = nBase+1
  end if
  if (iglobal == 1) then
    Base(i) = Line
    Base(i)(1:2) = adjustl(PTab(i))
    if (Base(i)(2:2) == ' ') then
      TMP = Base(i)(3:)
      Base(i)(2:) = TMP(1:47)
    end if
    BasAva(i) = .true.
    Found = .true.
    nBase = nBase+1
  end if

end do
if (.not. Found) then
  iErr = 1
  write(LuWr,*) ' [BasisReader]: Wrong symbol in line'
  write(LuWr,*) '                ',Line
  goto 9999
end if
icount = icount+1
if (icount <= nxbas) goto 10

!9904 iErr = 1
!write(LuWr,*) ' [BasisReader]: Unable to read z-matrix !'
9999 continue

return

end subroutine BasisReader
