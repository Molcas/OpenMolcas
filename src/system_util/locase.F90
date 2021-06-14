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
! Copyright (C) Per Ake Malmqvist                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine lowercases a text string.                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-AAke Malmqvist                                          *
!          Lund University                                             *
!                                                                      *
!***********************************************************************

subroutine LoCase(string)

character*(*) string
character*26 up, lw
dimension itab(0:255)
save up, lw, ifset, itab
data up/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
data lw/'abcdefghijklmnopqrstuvwxyz'/
data ifset/0/

if (ifset == 0) then
  ifset = 1
  do i=0,255
    itab(i) = i
  end do
  do ii=1,26
    i = ichar(up(ii:ii))
    j = ichar(lw(ii:ii))
    itab(i) = j
  end do
end if

do ii=1,len(string)
  i = ichar(string(ii:ii))
  j = itab(i)
  string(ii:ii) = char(j)
end do

return

end subroutine LoCase
