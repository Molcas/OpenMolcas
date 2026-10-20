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

subroutine par_range(n,i,j)
! Distribute 'n' evenly over processes and return
! the range (i,j) of this particular process.
! If there is no valid range, then j<i will be returned,
! so it is possible to use it consistently for looping.

use Para_Info, only: MyRank, nProcs
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(out) :: i, j
integer(kind=iwp) :: nqot, nrem

nqot = n/nprocs
nrem = n-nqot*nprocs
if (myrank < nrem) then
  i = myrank*(nqot+1)+1
  j = i+nqot
else
  i = nrem*(nqot+1)+(myrank-nrem)*nqot+1
  j = i+nqot-1
end if

end subroutine par_range
