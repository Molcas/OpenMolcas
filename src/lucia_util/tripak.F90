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
subroutine TRIPAK(AUTPAK,APAK,MATDIM,NDIM)
! (NOT A SIMPLIFIED VERSION OF TETRAPAK)
!
! REFORMATING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
!
! FULL TO PACKED

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: MATDIM, NDIM
real(kind=wp), intent(in) :: AUTPAK(MATDIM,MATDIM)
real(kind=wp), intent(_OUT_) :: APAK(*)
integer(kind=iwp) :: I, IJ

IJ = 0
do I=1,NDIM
  APAK(IJ+1:IJ+I) = AUTPAK(1:I,I)
  IJ = IJ+I
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' AUTPAK AND APAK FROM TRIPAK'
call WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
call PRSYM(APAK,NDIM)
#endif

end subroutine TRIPAK
