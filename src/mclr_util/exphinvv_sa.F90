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

subroutine ExpHinvv_SA(n,DIAG,Vec1,Vec3,E0)
! The CI preconditioner of the SA response for the state of energy E0:
!
!                    -1
! |Vec3> = |H  - E0|   |Vec1>
!            0
!
! H0 is the explicit Hamiltonian in the subspace iSel (ExplH_SA, through
! its eigenvectors ExplV and eigenvalues ExplE) and the diagonal DIAG
! elsewhere, as in the Davidson procedure of RASSCF (David5). The absolute
! values keep the preconditioner positive definite for excited states
! (the CI block of the SA Hessian, 2 w_i (H-E_i), is positive definite on
! the space orthogonal to the SA states). The same operator must be used
! for S (SA_PREC2) and in DMinvCI_SA.

use MCLR_Data, only: ExplE, ExplV, iSel, nSel
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: DIAG(n), Vec1(n), E0
real(kind=wp), intent(out) :: Vec3(n)
integer(kind=iwp) :: i, iConf
real(kind=wp) :: Z
real(kind=wp), allocatable :: Scr1(:), Scr3(:)
real(kind=wp), parameter :: ThrDiag = 1.0e-16_wp, ThrExpl = 1.0e-3_wp

! the diagonal part
do i=1,n
  Z = abs(E0-DIAG(i))
  if (Z < ThrDiag) Z = ThrDiag
  Vec3(i) = Vec1(i)/Z
end do

! the elements in the subspace of the explicit Hamiltonian
if (nSel > 1) then
  call mma_allocate(Scr1,nSel,Label='Scr1')
  call mma_allocate(Scr3,nSel,Label='Scr3')
  do i=1,nSel
    iConf = iSel(i)
    Scr1(i) = Vec1(iConf)
  end do
  call DGEMM_('T','N',nSel,1,nSel,One,ExplV,nSel,Scr1,nSel,Zero,Scr3,nSel)
  do i=1,nSel
    Z = abs(E0-ExplE(i))
    if (Z < ThrExpl) Z = ThrExpl
    Scr3(i) = Scr3(i)/Z
  end do
  call DGEMM_('N','N',nSel,1,nSel,One,ExplV,nSel,Scr3,nSel,Zero,Scr1,nSel)
  do i=1,nSel
    iConf = iSel(i)
    Vec3(iConf) = Scr1(i)
  end do
  call mma_deallocate(Scr1)
  call mma_deallocate(Scr3)
end if

end subroutine ExpHinvv_SA
