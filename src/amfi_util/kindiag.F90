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

subroutine kindiag(TKIN,ndim,evec,eval,breit)
!bs determines eigenvectors and -values of TKIN

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndim
real(kind=wp), intent(in) :: tkin(ndim,ndim)
real(kind=wp), intent(out) :: evec(ndim,ndim), eval(ndim)
logical(kind=iwp), intent(in) :: breit
integer(kind=iwp) :: irun, irun1, irun2, itria, JRUN
real(kind=wp) :: fact
real(kind=wp), allocatable :: TKINTRIA(:)

call mma_allocate(TKINTRIA,ndim*(ndim+1)/2,label='TKINTRIA')

!bs move symmetric matrix to triangular matrix
itria = 1
do irun2=1,ndim
  do irun1=1,irun2
    TKINTRIA(itria) = TKIN(irun1,irun2)
    itria = itria+1
  end do
end do
call unitmat(evec,ndim)
!bs now diagonalize
call Jacob(TKINTRIA,evec,ndim,ndim)
!bs get the eigenvalues
if (breit) then
  eval(:) = Zero
else
  do irun=1,ndim
    eval(irun) = TKINTRIA(irun*(irun+1)/2)
  end do
end if
call mma_deallocate(TKINTRIA)
!bs ensure normalization of the vectors.
do IRUN=1,ndim
  fact = Zero
  do JRUN=1,ndim
    fact = fact+evec(JRUN,IRUN)*evec(JRUN,IRUN)
  end do
  fact = One/sqrt(fact)
  evec(:,IRUN) = fact*evec(:,IRUN)
end do

return

end subroutine kindiag
