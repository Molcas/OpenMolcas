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
subroutine TRIPK2_1(AUTPAK,APAK,MATDIM,NDIM)
! REFORMATTING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
!
! FULL TO PACKED
! LOWER HALF OF AUTPAK IS STORED IN APAK
! NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: MATDIM, NDIM
real(kind=wp), intent(in) :: AUTPAK(MATDIM,MATDIM)
real(kind=wp), intent(_OUT_) :: APAK(*)
integer(kind=iwp) :: IJ, J

IJ = 0
do J=1,NDIM
  APAK(IJ+J:IJ+NDIM) = AUTPAK(J:NDIM,J)
  IJ = IJ+NDIM-J
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' AUTPAK AND APAK FROM TRIPK2_1'
call WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
call PRSM2(APAK,NDIM)
#endif

end subroutine TRIPK2_1
