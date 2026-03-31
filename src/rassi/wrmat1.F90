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

subroutine WRMAT1(ND1,ND2,XMAT)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ND1, ND2
real(kind=wp), intent(in) :: XMAT(ND1,ND2)
integer(kind=iwp) :: I, IBL, J, JEND, JSTA, NBL
integer(kind=iwp), parameter :: NCOL = 5

! NCOL=NR OF PRINTING COLUMNS.
NBL = (ND2+NCOL-1)/NCOL
do IBL=1,NBL
  JSTA = 1+NCOL*(IBL-1)
  JEND = min(NCOL*IBL,ND2)
  write(u6,'(//,5(8X,I8),/)') (J,J=JSTA,JEND)
  do I=1,ND1
    write(u6,'(1X,I3,5(1X,G16.9))') I,(XMAT(I,J),J=JSTA,JEND)
  end do
end do

end subroutine WRMAT1
