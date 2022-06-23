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

subroutine Assemble_mGauss(As,Ad,nAs)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAs
real(kind=wp) :: As(nAs), Ad(nAs,6)

call DaXpY_(nAs,One,Ad(1,1),1,As,1)
call DaXpY_(nAs,One,Ad(1,4),1,As,1)
call DaXpY_(nAs,One,Ad(1,6),1,As,1)

return

end subroutine Assemble_mGauss
