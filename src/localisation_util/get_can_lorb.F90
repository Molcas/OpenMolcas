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

subroutine Get_Can_Lorb(Ene,Fock,nO,nX,jOrb,Umat)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Ene(*)
real(kind=wp), intent(inout) :: Fock(*), Umat(*)
integer(kind=iwp), intent(in) :: nO, nX, jOrb(nO)
integer(kind=iwp) :: i, nOx, nXx
real(kind=wp), allocatable, target :: eta_t(:), Zt(:)
real(kind=wp), pointer :: eta(:,:), Z(:,:)

if (nO < 1) return

call mma_allocate(eta_t,nX**2,label='eta_ik')
eta(1:nX,1:nX) => eta_t(1:nX**2)
call mma_allocate(Zt,max(nX,2)*nO,label='Zt')
Z(1:nX,1:nO) => Zt(1:nX*nO)
eta(:,:) = Zero
do i=1,nX
  eta(i,i) = Ene(i)
end do

nXx = max(1,nX)
nOx = max(1,nO)
call DGEMM_('N','N',nX,nO,nX,One,eta,nXx,Umat,nXx,Zero,Z,nXx)
eta(1:nO,1:nO) => eta_t(1:nO**2)
call DGEMM_('T','N',nO,nO,nX,One,Umat,nXx,Z,nXx,Zero,eta,nOx)

Z(1:nO,1:2) => Zt(1:nO*2)
call Eigen_Molcas(nO,eta,Z(:,1),Z(:,2))

call dcopy_(nO**2,eta,1,Umat,1)
do i=1,nO
  Fock(jOrb(i)) = Z(i,1)
end do

nullify(eta)
nullify(Z)
call mma_deallocate(eta_t)
call mma_deallocate(Zt)

return

end subroutine Get_Can_Lorb
