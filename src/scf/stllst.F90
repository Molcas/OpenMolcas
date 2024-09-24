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

implicit none
integer LLink, iRoot

!return
write(6,*)
write(6,*) '*********** Status of Linked List *************'
write(6,*)
write(6,*) ' LLink:',LLink
write(6,*)
write(6,*) ' CNOD data'
write(6,*) 'Error code:                       ',nLList(LLink,0)
write(6,*) 'Pointer to first NODE in the list:',nLList(LLink,1)
write(6,*) 'Actual length of list:            ',nLList(LLink,2)
write(6,*) '# of vectors in core:             ',nLList(LLink,3)
write(6,*)
iRoot = nLList(LLink,1)
do while (iRoot /= 0)
  write(6,*) ' NODE data'
  write(6,*) 'NODE @:                         ',iRoot
  write(6,*) 'Pointer to next NODE:           ',nLList(iRoot,0)
  write(6,*) 'Pointer to stored vector:       ',nLList(iRoot,1)
  if (nLList(iRoot,5) >= 1) then
    write(6,*) 'Vector status:                  in Core'
  else
    write(6,*) 'Vector status:                  on Disk'
  end if
  write(6,*) 'Next free position:             ',nLList(iRoot,2)
  write(6,*) 'Length of vector:               ',nLList(iRoot,3)
  write(6,*) 'Iteration number:               ',nLList(iRoot,4)
  write(6,*)
  iRoot = nLList(iRoot,0)
end do
write(6,*) '************ End of Status Report *************'
write(6,*)

return

end subroutine StlLst
