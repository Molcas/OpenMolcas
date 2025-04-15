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
subroutine TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SGN)
! REFORMATTING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
!
! IWAY = 1 : FULL TO PACKED
!            LOWER HALF OF AUTPAK IS STORED IN APAK
! IWAY = 2 : PACKED TO FULL FORM
!            APAK STORED IN LOWER HALF
!             SGN * APAK TRANSPOSED IS STORED IN UPPPER PART
! NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IWAY, MATDIM, NDIM
real(kind=wp) :: AUTPAK(MATDIM,MATDIM), APAK(*), SGN
integer(kind=iwp) :: IJ, J

if (IWAY == 1) then
  IJ = 0
  do J=1,NDIM
    APAK(IJ+J:IJ+NDIM) = AUTPAK(J:NDIM,J)
    IJ = IJ+NDIM-J
  end do
else if (IWAY == 2) then
  IJ = 0
  do J=1,NDIM
    AUTPAK(J,J:NDIM) = SGN*APAK(IJ+J:IJ+NDIM)
    AUTPAK(J:NDIM,J) = APAK(IJ+J:IJ+NDIM)
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
