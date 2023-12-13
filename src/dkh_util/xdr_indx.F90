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

subroutine xdr_indx(N,indx)
! read atomic/block information for local transformation

use Definitions, only: iwp

integer(kind=iwp), intent(in) :: N
integer(kind=iwp), intent(out) :: indx(n)

call get_iarray('Ctr Index Prim',indx,N)

end subroutine xdr_indx
