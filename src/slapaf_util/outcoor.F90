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

subroutine OutCoor(TEXT,Chr,nDim,FI,N1,N2,lAngstroms)
!***********************************************************************
!                                                                      *
!     Object: To generate a cartesian output with atomic labels        *
!             N1 and N2 are the real limits of dummy FI                *
!                                                                      *
!***********************************************************************

use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDim, N1, N2
character(len=*), intent(in) :: TEXT, Chr(nDim)
real(kind=wp), intent(in) :: FI(N1,N2)
logical(kind=iwp), intent(in) :: lAngstroms
integer(kind=iwp) :: I, Lu

Lu = u6
write(Lu,*)
write(Lu,*) '*********************************************************'
write(Lu,*) Text
write(Lu,*) '*********************************************************'
write(Lu,*) ' ATOM              X               Y               Z     '
do I=1,NDIM
  if (lAngstroms) then
    write(Lu,300) Chr(I),FI(1:3,I)*Angstrom
  else
    write(Lu,300) Chr(I),FI(1:3,I)
  end if
end do

write(Lu,*)

return

300 format(2X,A,3X,3F16.6)

end subroutine OutCoor
