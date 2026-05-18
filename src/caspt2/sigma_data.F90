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

module Sigma_data

use definitions, only: iwp, wp

implicit none

real(kind=wp) :: VAL1(2), VAL2(2)
integer(kind=iwp) :: INCX1, INCX2, INCX3, INCF1, INCF2, INCY1, INCY2, INCY3, LEN1, LEN2, NLST1, NLST2
integer(kind=iwp) :: NFSCA, NFDXP, NFMV, NFR1, IFTEST

end module Sigma_data
