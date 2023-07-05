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

subroutine PrList(TEXT,Char,nDim,FI,N1,N2)
!***********************************************************************
!                                                                      *
!     Object: To generate a cartesian output with atomic labels        *
!             N1 and N2 are the real limits of dummy FI                *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
character*(*) TEXT
character*(*) char(nDim)
real*8 FI(N1,N2)

Lu = 6
write(Lu,100) Text
write(Lu,200)
do I=1,NDIM
  if (N1 == 3) then
    write(Lu,300) char(I),(FI(J,I),J=1,3)
  else
    write(Lu,300) char(I),(FI(I,J),J=1,3)
  end if
end do

return

100 format(//,1X,A,/)
200 format(5X,'ATOM',21X,'X',19X,'Y',19X,'Z',/)
300 format(5X,A,3X,3F20.10)

end subroutine PrList
