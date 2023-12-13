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

subroutine ElRed2(nq,nx,Gmtrx,EVal,EVec,nK,uMtrx,g12K,Thr,BM,iBM,nB_Tot,nqB)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nq, nx, nB_Tot, iBM(nB_Tot), nqB(nq)
real(kind=wp), intent(out) :: Gmtrx(nq,nq), EVal(nTri_Elem(nq)), EVec(nq,nq)
integer(kind=iwp), intent(out) :: nK
real(kind=wp), intent(in) :: uMtrx(nX), Thr, BM(nB_Tot)
logical(kind=iwp), intent(in) :: g12K
integer(kind=iwp) :: i, iB, ii, ijTri, ik, Info, iQ, j, jB, jl, k, l, LDZ, mB, nB
real(kind=wp) :: B_ik, B_jl, rSum
logical(kind=iwp) :: Diagonal
real(kind=wp), allocatable :: W(:), Work(:)
real(kind=wp), parameter :: Zero_Approx = 1.0e-10_wp

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *

#ifdef _DEBUGPRINT_
call RecPrt('ElRed2: The u matrix','(5es21.12)',umtrx,nx,1)
#endif
if (nq == 0) then
  nK = 0
  return
end if

!                         T
! Form the G matrix, G=BuB

GMtrx(:,:) = Zero
iB = 0
do i=1,nq
  nB = nqB(i)
  do k=1,nB
    iB = iB+1
    ik = iBM(iB)
    B_ik = BM(iB)

    jB = 0
    do j=1,nq
      mB = nqB(j)
      do l=1,mB
        jB = jB+1
        jl = iBM(jB)
        if (ik == jl) then
          B_jl = BM(jB)
          GMtrx(i,j) = GMtrx(i,j)+B_ik*umtrx(ik)*B_jl
        end if
      end do
    end do

  end do
end do

Diagonal = .true.
do i=1,nq
  rSum = Zero
  do j=1,nq
    if (abs(Gmtrx(i,j)) < 1.0e-10_wp) Gmtrx(i,j) = Zero
    if (j /= i) rSum = rSum+GMtrx(i,j)
  end do
  Diagonal = Diagonal .and. (rSum == Zero)
end do

#ifdef _DEBUGPRINT_
call RecPrt('ElRed2: The G Matrix (nq x nq)','(5es21.12)',Gmtrx,nq,nq)
write(u6,*) 'Diagonal=',Diagonal
#endif

! Set up a unit matrix

call unitmat(EVec,nq)

! Set up the Hessian in lower triangular form, the elements are symmetrized.

do i=1,nQ
  do j=1,i
    ijTri = nTri_Elem(i-1)+j
    EVal(ijTri) = Half*(Gmtrx(i,j)+Gmtrx(j,i))
  end do
end do
#ifdef _DEBUGPRINT_
call TriPrt('Eval prediagonalization',' ',EVal,nQ)
#endif

!                                                    |g 0|
! Compute eigenvalues and eigenvectors G(K L) = (K L)|0 0|
! K: nonredundant vectors with eigenvalues g
! L: redundant vectors

if (.not. Diagonal) then
  LDZ = max(1,nq)
  call mma_allocate(Work,3*nq,Label='Work')
  Work(:) = Zero
  call mma_allocate(W,nq,Label='W')
  W(:) = Zero
  Info = 0
  call dspev_('V','U',nq,Eval,W,EVec,LDZ,Work,Info)
  if (Info /= 0) then
    write(u6,*) 'Info /= 0'
    write(u6,*) 'Info=',Info
    call Abend()
  end if
  EVal(:) = Zero
  do i=1,nq
    EVal(nTri_Elem(i)) = W(i)
  end do
  call mma_deallocate(W)
  call mma_deallocate(Work)
end if
EVal(:) = -Eval(:)
call JacOrd(EVal,EVec,nQ,nQ)
!Fix standard direction.
do iQ=1,nQ
  call VecPhase(EVec(:,iQ),nQ)
end do
EVal(:) = -Eval(:)
#ifdef _DEBUGPRINT_
call RecPrt('ElRed2: Eigenvectors',' ',EVec,nQ,nQ)
call TriPrt('ElRed2: Eigenvalues',' ',EVal,nQ)
#endif

!                                    -1/2
! Remove redundant vectors and form g     K

nK = 0
do i=1,nQ
  ii = nTri_Elem(i)
  if (EVal(ii) > Thr) nK = nK+1
  EVal(i) = EVal(ii)
  !if (g12K .and. (abs(EVal(i)) > Zero))
  if (g12K .and. (abs(EVal(i)) > Zero_Approx)) EVec(:,i) = EVec(:,i)/sqrt(EVal(i))
end do
#ifdef _DEBUGPRINT_
call RecPrt('ElRed2: The NonRedundant eigenvectors','(5es21.12)',EVec,nQ,nK)
call RecPrt('ElRed2: eigenvalues ','(8ES12.4)',EVal,1,nK)
#endif

return

end subroutine ElRed2
