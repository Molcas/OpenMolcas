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

subroutine ZSYM1(NIRREP)

use symmetry_info, only: Mul
use csm_data, only: ADASX, ADSXA, ASXAD, ITSDX, ITSSX, NSMCI, NSMDX, NSMST, NSMSX, SXDXSX, SXSXDX
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: NIRREP

NSMSX = NIRREP
NSMDX = NIRREP
NSMST = NIRREP
NSMCI = NIRREP
ITSSX = 1
ITSDX = 1

ADASX(:,:) = 0
ADASX(1:8,1:8) = Mul(:,:)
ADSXA(:,:) = 0
ADSXA(1:8,1:8) = Mul(:,:)
ASXAD(:,:) = 0
ASXAD(1:8,1:8) = Mul(:,:)
SXSXDX(:,:) = 0
SXSXDX(1:8,1:8) = Mul(:,:)
SXDXSX(:,:) = 0
SXDXSX(1:8,1:8) = Mul(:,:)

end subroutine ZSYM1
