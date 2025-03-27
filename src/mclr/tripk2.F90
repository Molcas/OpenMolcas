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

!#define _DEBUGPRINT_
subroutine TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
! REFORMATTING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
!
! IWAY = 1 : FULL TO PACKED
!            LOWER HALF OF AUTPAK IS STORED IN APAK
! IWAY = 2 : PACKED TO FULL FORM
!            APAK STORED IN LOWER HALF
!             SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
! NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS

#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit real*8(A-H,O-Z)
dimension AUTPAK(MATDIM,MATDIM), APAK(*)

if (IWAY == 1) then
  IJ = 0
  do J=1,NDIM
    do I=J,NDIM
      APAK(IJ+I) = AUTPAK(I,J)
    end do
    IJ = IJ+NDIM-J
  end do
else if (IWAY == 2) then
  IJ = 0
  do J=1,NDIM
    do I=J,NDIM
      AUTPAK(J,I) = SIGN*APAK(IJ+I)
      AUTPAK(I,J) = APAK(IJ+I)
    end do
    IJ = IJ+NDIM-J
  end do
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' AUTPAK AND APAK FROM TRIPK2'
call WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
call PRSM2(APAK,NDIM)
#endif

return

end subroutine TRIPK2
