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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine Get_Nat_Lorb(Occ,FOcc,nO,nX,jOrb,Umat)

#include "intent.fh"

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Occ(*)
real(kind=wp), intent(_OUT_) :: FOcc(*)
real(kind=wp), intent(inout) :: Umat(*)
integer(kind=iwp), intent(in) :: nO, nX, jOrb(nO)
integer(kind=iwp) :: i, ii, j, nOx, nXx
real(kind=wp), allocatable :: eta(:), Z(:), ZZ(:)

if (nO < 1) return

call mma_allocate(eta,nX**2,label='eta_ik')
call mma_allocate(Z,nX**2,label='Z')
call mma_allocate(ZZ,nX,label='ZZ')
eta(:) = Zero
do i=1,nX
  ii = nX*(i-1)+i
  eta(ii) = Occ(i)
end do
nXx = max(1,nX)
nOx = max(1,nO)
call DGEMM_('N','N',nX,nO,nX,One,eta,nXx,Umat,nXx,Zero,Z,nXx)
call DGEMM_('T','N',nO,nO,nX,One,Umat,nXx,Z,nXx,Zero,eta,nOx)

call Eigen_Molcas(nO,eta,ZZ,Z)

call dcopy_(nO**2,eta,1,Umat,1)
do i=1,nO
  j = jOrb(i)
  FOcc(j) = ZZ(i)
end do
call mma_deallocate(eta)
call mma_deallocate(Z)
call mma_deallocate(ZZ)

return

end subroutine Get_Nat_Lorb
