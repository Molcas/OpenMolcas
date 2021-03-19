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

subroutine PrintLine(unit,line,len,isBinLuscus)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"
character line*128
integer unit
integer len, ll, li

if (ISLUSCUS == 1) then
  ll = len
  li = isBinLuscus
  !write(6,*) 'before pl ',line,' ',ll,li
  call prt_lusc(unit,line,ll,li)
else
  if (isBinary == 1) then
    write(unit) line(1:len)
  else
    write(unit,'(A)') line(1:len)
  end if
end if

return

end subroutine PrintLine
