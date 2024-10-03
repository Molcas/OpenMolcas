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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003-2005, Valera Veryazov                             *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine StlLst(LLink)

use LnkLst, only: nLList
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LLink
integer(kind=iwp) :: iRoot

!return
write(u6,*)
write(u6,*) '*********** Status of Linked List *************'
write(u6,*)
write(u6,*) ' LLink:',LLink
write(u6,*)
write(u6,*) ' CNOD data'
write(u6,*) 'Error code:                       ',nLList(LLink,0)
write(u6,*) 'Pointer to first NODE in the list:',nLList(LLink,1)
write(u6,*) 'Actual length of list:            ',nLList(LLink,2)
write(u6,*) '# of vectors in core:             ',nLList(LLink,3)
write(u6,*)
iRoot = nLList(LLink,1)
do while (iRoot /= 0)
  write(u6,*) ' NODE data'
  write(u6,*) 'NODE @:                         ',iRoot
  write(u6,*) 'Pointer to next NODE:           ',nLList(iRoot,0)
  write(u6,*) 'Pointer to stored vector:       ',nLList(iRoot,1)
  if (nLList(iRoot,5) >= 1) then
    write(u6,*) 'Vector status:                  in Core'
  else
    write(u6,*) 'Vector status:                  on Disk'
  end if
  write(u6,*) 'Next free position:             ',nLList(iRoot,2)
  write(u6,*) 'Length of vector:               ',nLList(iRoot,3)
  write(u6,*) 'Iteration number:               ',nLList(iRoot,4)
  write(u6,*)
  iRoot = nLList(iRoot,0)
end do
write(u6,*) '************ End of Status Report *************'
write(u6,*)

return

end subroutine StlLst
