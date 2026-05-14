!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Roland Lindh                                     *
!               2014, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine QNR(nInter,nIter,dq,H,g)
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nInter, nIter
real(kind=wp), intent(inout) :: dq(nInter,nIter)
real(kind=wp), intent(in) :: H(nInter,nInter), g(nInter,nIter+1)
integer(kind=iwp) :: Info, k, maxk, recomp
real(kind=wp) :: alpha, beta, rr, RelThr
real(kind=wp), allocatable :: Ap(:), Lo(:,:), p(:), r(:), y(:), z(:)
real(kind=wp), parameter :: Thr = 1.0e-20_wp
real(kind=wp), external :: ddot_

! Compute a new independent geometry by relaxation of
! the gradient vector.

call mma_allocate(r,nInter)
call mma_allocate(p,nInter)
call mma_allocate(Ap,nInter)
call mma_allocate(z,nInter)
call mma_allocate(y,nInter)

dq(:,nIter) = Zero

maxk = max(10,nInter**2)
recomp = max(50,int(nInter/Ten))
r(:) = g(:,nIter)

call mma_allocate(Lo,nInter,nInter,label='Lo')

! With a dense matrix, the preconditioner could be replaced with
! something else, otherwise this is just solving the system with
! a direct method
Lo(:,:) = H(:,:)
call dpotrf_('L',nInter,Lo,nInter,info)

call DSyMV('L',nInter,-One,H,nInter,dq(:,nIter),1,One,r,1)
z(:) = r(:)
call DPoTrS('L',nInter,1,Lo,nInter,z,nInter,info)
p(:) = z(:)
rr = DDot_(nInter,z,1,r,1)
RelThr = Thr*max(rr,One)
k = 1
do while ((abs(rr) >= RelThr) .and. (k <= maxk))
  call dGeMV_('nInter',nInter,nInter,One,H,nInter,p,1,Zero,Ap,1)
  alpha = rr/DDot_(nInter,p,1,Ap,1)
  dq(:,nIter) = dq(:,nIter)+alpha*p(:)
  beta = rr
  if (mod(k,recomp) == 0) then
    r(:) = g(:,nIter)
    call dGeMV_('N',nInter,nInter,-One,H,nInter,dq(:,nIter),1,One,r,1)
  else
    r(:) = r(:)-alpha*Ap(:)
  end if
  z(:) = r(:)
  call DPoTrS('L',nInter,1,Lo,nInter,z,nInter,info)
  rr = DDot_(nInter,z,1,r,1)
  p(:) = p(:)*rr/beta+z(:)
  k = k+1
end do
call mma_deallocate(Lo)

! Set the return value
if (k <= maxk) then
  info = k
else
  info = -1
end if

call mma_deallocate(r)
call mma_deallocate(p)
call mma_deallocate(Ap)
call mma_deallocate(z)
call mma_deallocate(y)

if (Info < 0) call SysAbendMsg('QNR','Conjugate gradients not converged',' ')

#ifdef _DEBUGPRINT_
write(u6,*) 'CG converged in ',Info,' iterations'
call RecPrt(' H ',' ',H,nInter,nInter)
call RecPrt(' dq',' ',dq(1,nIter),1,nInter)
#endif

return

end subroutine QNR
