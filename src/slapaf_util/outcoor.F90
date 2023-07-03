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

subroutine OutCoor(TEXT,Char,nDim,FI,N1,N2,lAngstroms)
!***********************************************************************
!                                                                      *
!     Object: To generate a cartesian output with atomic labels        *
!             N1 and N2 are the real limits of dummy FI                *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
#include "angstr.fh"
character*(*) TEXT
character*(*) char(nDim)
logical lAngstroms
real*8 FI(N1,N2)

Lu = 6
write(Lu,*)
write(Lu,*) '*********************************************************'
write(Lu,*) Text
write(Lu,*) '*********************************************************'
write(Lu,*) ' ATOM              X               Y               Z     '
do I=1,NDIM
  if (lAngstroms) then
    write(Lu,300) char(I),(FI(J,I)*angstr,J=1,3)
  else
    write(Lu,300) char(I),(FI(J,I),J=1,3)
  end if
300 format(2X,A,3X,3F16.6)
end do

write(Lu,*)

return

end subroutine OutCoor
