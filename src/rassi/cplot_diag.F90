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

subroutine CPLOT_DIAG(MATR,MATI,nDIM,EIGVECR,EIGVECI)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDIM
real(kind=wp), intent(inout) :: MATR(nDIM*(nDIM+1)/2), MATI(nDIM*(nDIM+1)/2)
real(kind=wp), intent(out) :: EIGVECR(nDIM,nDIM), EIGVECI(nDIM,nDIM)

integer(kind=iwp) :: I, INFO, J
real(kind=wp), allocatable :: CEIGVAL(:), RWORK(:)
complex(kind=wp), allocatable :: CEIGVEC(:,:), MATFULL(:), ZWORK(:)

call mma_allocate(MATFULL,nDIM*(nDIM+1)/2,Label='MATFULL')
call mma_allocate(CEIGVAL,nDIM,Label='CEIGVAL')
call mma_allocate(CEIGVEC,nDIM,nDIM,Label='CEIGVEC')
call mma_allocate(ZWORK,2*nDIM-1,Label='ZWORK')
call mma_allocate(RWORK,3*nDIM-2,Label='RWORK')

do J=1,(nDIM*(nDIM+1)/2)
  MATFULL(J) = cmplx(MATR(J),MATI(J),kind=wp)
end do

call zhpev_('V','U',nDIM,MATFULL,CEIGVAL,CEIGVEC,nDIM,ZWORK,RWORK,INFO)
call mma_deallocate(MATFULL)
call mma_deallocate(ZWORK)
call mma_deallocate(RWORK)

if (INFO /= 0) then
  write(u6,*) 'Error in diagonalization'
  write(u6,*) 'INFO: ',INFO
  call ABEND()
end if

do I=1,nDIM
  do J=1,nDIM
    EIGVECR(I,J) = real(CEIGVEC(I,J))
    EIGVECI(I,J) = aimag(CEIGVEC(I,J))
  end do
end do

call DCOPY_(nDIM*(nDIM+1)/2,[Zero],0,MATR,1)
call DCOPY_(nDIM*(nDIM+1)/2,[Zero],0,MATI,1)

do J=1,nDIM
  MATR((J*(J-1)/2)+J) = CEIGVAL(J)
end do

call mma_deallocate(CEIGVAL)
call mma_deallocate(CEIGVEC)

end subroutine CPLOT_DIAG
