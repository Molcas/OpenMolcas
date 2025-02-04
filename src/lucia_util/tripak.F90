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

use Definitions, only: u6

implicit real*8(A-H,O-Z)
dimension AUTPAK(MATDIM,MATDIM), APAK(*)

if (IWAY == 1) then
  IJ = 0
  do I=1,NDIM
    do J=1,I
      APAK(IJ+J) = AUTPAK(J,I)
    end do
    IJ = IJ+I
  end do
else if (IWAY == 2) then
  IJ = 0
  do I=1,NDIM
    do J=1,I
      AUTPAK(I,J) = APAK(IJ+J)
      AUTPAK(J,I) = APAK(IJ+J)
    end do
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
