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

subroutine TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
! (NOT A SIMPLIFIED VERSION OF TETRAPAK)
!
! REFORMATING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
!
! IWAY = 1 : FULL TO PACKED
! IWAY = 2 : PACKED TO FULL FORM

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IWAY, MATDIM, NDIM
real(kind=wp), intent(inout) :: AUTPAK(MATDIM,MATDIM), APAK(*)
integer(kind=iwp) :: I, IJ, NTEST

if (IWAY == 1) then
  IJ = 0
  do I=1,NDIM
    APAK(IJ+1:IJ+I) = AUTPAK(1:I,I)
    IJ = IJ+I
  end do
else if (IWAY == 2) then
  IJ = 0
  do I=1,NDIM
    AUTPAK(I,1:I-1) = APAK(IJ+1:IJ+I-1)
    AUTPAK(1:I,I) = APAK(IJ+1:IJ+I)
    IJ = IJ+I
  end do
end if

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' AUTPAK AND APAK FROM TRIPAK'
  call WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
  call PRSYM(APAK,NDIM)
end if

end subroutine TRIPAK
