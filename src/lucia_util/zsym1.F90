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
use lucia_data, only: MXPOBS
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

call ICPMT2(Mul,ADASX,8,8,MXPOBS,MXPOBS,1)
call ICPMT2(Mul,ADSXA,8,8,MXPOBS,2*MXPOBS,1)
call ICPMT2(Mul,ASXAD,8,8,MXPOBS,2*MXPOBS,1)
call ICPMT2(Mul,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
call ICPMT2(Mul,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)

end subroutine ZSYM1
