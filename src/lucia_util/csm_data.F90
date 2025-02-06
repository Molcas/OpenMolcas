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

module CSM_data

! NSMSX  : Number of symmetries single ex    (nIrrep)
! NSMDX  : Number of symmetries double ex    (nIrrep)
! NSMST  : Number of symmetries string       (nIrrep)
! NSMCI  : Number of symmetries CI Space     (nIrrep)
! ITSSX  : Total symmetrix single excitation (1)
! ITSDX  : Total symmetrix double excitation (1)
!
! ADASX  : symmetry operator initialized in syminf
! ASXAD  :                   initialized in syminf
! ADSXA  :                   initialized in syminf
! SXSXDX :                   initialized in syminf
! SXDXSX :                   initialized in syminf

use lucia_data, only: MXPOBS
use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: ADASX(MXPOBS,MXPOBS), ADSXA(MXPOBS,2*MXPOBS), ASXAD(MXPOBS,2*MXPOBS), ITSDX, ITSSX, NSMCI, NSMDX, NSMST, &
                     NSMSX, SXDXSX(2*MXPOBS,4*MXPOBS), SXSXDX(2*MXPOBS,2*MXPOBS)

public :: ADASX, ADSXA, ASXAD, ITSDX, ITSSX, NSMCI, NSMDX, NSMST, NSMSX, SXDXSX, SXSXDX

end module CSM_data
