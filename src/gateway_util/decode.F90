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
! Copyright (C) Bjorn O. Roos                                          *
!               1991, Roland Lindh                                     *
!***********************************************************************

subroutine Decode(LBL,string,N,Hit)
!***********************************************************************
! Object: to find the character string 'string' between                *
!         dots N-1 and N in the label LBL of length lLBL               *
!         blanks in string are removed                                 *
! Called from: Rdbsl                                                   *
! Subroutine calls: none                                               *
!                                                                      *
! Author: Bjoern Roos, University of Lund, Sweden                      *
!***********************************************************************

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: LBL
character(len=*), intent(out) :: string
integer(kind=iwp), intent(in) :: N
logical(kind=iwp), intent(inout) :: Hit
integer(kind=iwp) :: i, i1, idot, lLBL, lstring
character(len=80) :: xstring
character, parameter :: dot = '.'

!write(u6,'(1x,a)') LBL
i1 = 1
idot = 0
lstring = 0
lLBL = len(LBL)
do i=1,lLBL
  if (LBL(i:i) /= dot) cycle
  idot = idot+1
  if (idot == N-1) i1 = i+1
  if (idot == N) then
    xstring = ' '
    !write(u6,'(1x,A,/,1X,A)') ' xstring=',xstring
    if (i > i1) xstring = LBL(i1:i-1)
    !write(u6,'(1x,A,/,1X,A)') ' xstring=',xstring
    lstring = i-i1
    exit
  end if
end do
if (idot == N) then
  ! Pack the string

  Hit = .true.
  i1 = 0
  string = ' '
  do i=1,lstring
    if (xstring(i:i) == ' ') cycle
    i1 = i1+1
    string(i1:i1) = xstring(i:i)
  end do
else if (Hit) then
  call WarningMessage(2,'Decode: error in basis set label')
  write(u6,'(A,A)') 'LBL=',LBL
  call Abend()
end if

return

end subroutine Decode
