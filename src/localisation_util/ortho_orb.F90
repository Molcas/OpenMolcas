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

subroutine Ortho_Orb(Xmo,Smat,nBas,nOrb2Loc,nPass,Test)
! Purpose: orthonormalization of orbitals according to
!
!          V = X^T*S*X
!          X <- X*V^(-1/2)
!
!          where S is the AO overlap matrix.
!          The orthonormalization is carried out nPass times.
!          After this routine, X will satisfy X^T*S*X=1.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: Xmo(*)
real(kind=wp), intent(in) :: Smat(*)
integer(kind=iwp), intent(in) :: nBas, nOrb2Loc, nPass
logical(kind=iwp), intent(in) :: Test
integer(kind=iwp) :: i, iPass, iTask, l_Scr, nB, nErr, nO2L
real(kind=wp) :: xNrm
real(kind=wp), allocatable :: Scr(:), V(:,:), VISqrt(:,:), VSqrt(:,:)
real(kind=wp), parameter :: Tol = 1.0e-10_wp
character(len=*), parameter :: SecNam = 'Ortho_Orb'
real(kind=wp), external :: ddot_

! Check for quick return.
! -----------------------

if (nPass < 1) return

! Allocations.
! ------------

l_Scr = 2*(nBas**2)+nBas*(nBas+1)/2 ! needed in SqrtMt
call mma_allocate(V,nOrb2Loc,nOrb2Loc,label='V')
call mma_allocate(VSqrt,nOrb2Loc,nOrb2Loc,label='VSqrt')
call mma_allocate(VISqrt,nOrb2Loc,nOrb2Loc,label='VISqrt')
call mma_allocate(Scr,l_Scr,label='Scr')

! Orthonormalization passes.
! --------------------------

do iPass=1,nPass

  ! Compute V = X^T*S*X.
  ! --------------------

  call GetUmat_Localisation(V,Xmo,Smat,Xmo,Scr,nBas,nOrb2Loc)

  ! Compute V^(-1/2).
  ! -----------------

  iTask = 2 ! compute sqrt as well as inverse sqrt
  call SqrtMt(V,nOrb2Loc,iTask,VSqrt,VISqrt,Scr)

  ! Compute orthonormal X <- X*V^(-1/2).
  ! ------------------------------------

  nB = max(nBas,1)
  nO2L = max(nOrb2Loc,1)
  l_Scr = nBas*nOrb2Loc
  Scr(1:l_Scr) = Xmo(1:l_Scr)
  call DGEMM_('N','N',nBas,nOrb2Loc,nOrb2Loc,One,Scr,nB,VISqrt,nO2L,Zero,Xmo,nB)

end do

! Test orthonormalization (i.e. V=1?).
! ------------------------------------

if (Test) then
  nErr = 0
  call GetUmat_Localisation(V,Xmo,Smat,Xmo,Scr,nBas,nOrb2Loc)
  do i=1,nOrb2Loc
    V(i,i) = V(i,i)-One
  end do
  xNrm = sqrt(dDot_(nOrb2Loc**2,V,1,V,1))
  if (xNrm > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: ||X^TSX - 1|| = ',xNrm
    nErr = nErr+1
  end if
  if (nErr /= 0) then
    write(u6,*) SecNam,': failure after ',nPass,' passes'
    call SysAbendMsg(SecNam,'Orthonormalization failure!',' ')
  end if
end if

! De-allocations.
! ---------------

call mma_deallocate(V)
call mma_deallocate(VSqrt)
call mma_deallocate(VISqrt)
call mma_deallocate(Scr)

return

end subroutine Ortho_Orb
