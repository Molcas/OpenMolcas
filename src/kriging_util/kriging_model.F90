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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine kriging_model()

use kriging_mod, only: blAI, blaAI, blavAI, blvAI, detR, dy, full_R, Index_PGEK, Kv, lh, m_t, mblAI, nPoints, nD, nInter_Eff, &
                       ordinary, Rones, sb, sbmev, sbO, variance, y
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#define _DPOSV_
#ifdef _DPOSV_
use Constants, only: Two
#endif

implicit none
integer(kind=iwp) :: i_eff, is, ie, ise, iee, i, INFO ! ipiv the pivot indices that define the permutation matrix
integer(kind=iwp), allocatable :: IPIV(:)
real(kind=wp), allocatable :: B(:), A(:,:)
real(kind=wp), external :: dDot_

! Prediagonalize the part of the matrix corresponing to the value-value block

#define _PREDIAG_
#ifdef _PREDIAG_
integer(kind=iwp) :: j, ij
real(kind=wp), allocatable :: U(:,:), HTri(:), UBIG(:,:), C(:,:), D(:)
real(kind=wp) :: temp
#endif

call mma_Allocate(B,m_t,label='B')
call mma_Allocate(A,m_t,m_t,label='A')
call mma_Allocate(IPIV,m_t,label='IPIV')

! Initiate B according to Eq. (6) of ref.
B(1:nPoints) = One
B(nPoints+1:) = Zero

! Initiate A according to Eq. (2) of ref.

#ifdef _DEBUGPRINT_
call RecPrt('f',' ',B,1,m_t)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      Code with prediagonalization of the value-value block

#ifdef _PREDIAG_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First diagonalize the value-value block
!
! U will contain the eigenvectors

call mma_allocate(U,nPoints,nPoints,label='U')
U(:,:) = Zero
do i=1,nPoints
  U(i,i) = One
end do
call mma_allocate(HTri,nPoints*(nPoints+1)/2,label='HTri')
do i=1,nPoints
  do j=1,i
    ij = i*(i-1)/2+j
    HTri(ij) = Full_R(i,j)
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('U',' ',U,nPoints,nPoints)
#endif
call nidiag_new(HTri,U,nPoints,nPoints,0)
call Jacord(HTri,U,nPoints,nPoints)

! Introduce canonical phase factor

do i=1,nPoints
  Temp = DDot_(nPoints,[One],0,U(1,i),1)
  U(1:nPoints,i) = U(1:nPoints,i)*sign(One,Temp)
end do
#ifdef _DEBUGPRINT_
call RecPrt('U',' ',U,nPoints,nPoints)
call TriPrt('HTri',' ',HTri,nPoints)
#endif

! Now set up an eigenvector matrix for the whole space.

call mma_Allocate(UBIG,m_t,m_t,label='UBig')
UBIG(:,:) = Zero
UBIG(1:nPoints,1:nPoints) = U(:,:)
do i=nPoints+1,m_t
  UBIG(i,i) = One
end do
!call RecPrt('UBIG',' ',UBig,m_t,m_t)

! Transform the covariance matrix to this basis

call mma_Allocate(C,m_t,m_t,label='C')
C = Zero
call dgemm_('N','N',m_t,m_t,m_t,One,Full_R,m_t,UBIG,m_t,Zero,C,m_t)
call dgemm_('T','N',m_t,m_t,m_t,One,UBIG,m_t,C,m_t,Zero,A,m_t)

! For safety measures set the value-value to zero and then reintroduce the
! eigenvalues from the previous diagonalization.

A(1:nPoints,1:nPoints) = Zero
do i=1,nPoints
  A(i,i) = HTri(i*(i+1)/2)
end do
#ifdef _DEBUGPRINT_
call RecPrt('U^TAU',' ',A,m_t,m_t)
#endif

! Set up the F-vector with ones and zeros

call mma_allocate(D,m_t,label='D')
D(1:nPoints) = One
D(nPoints+1:) = Zero

! Transform the vector to the new basis

B(:) = Zero
call dgemm_('T','N',m_t,1,m_t,One,UBIG,m_t,D,m_t,Zero,B,m_t)
#ifdef _DEBUGPRINT_
call RecPrt('U^TB',' ',B,1,m_t)
#endif

call mma_deallocate(HTri)
call mma_deAllocate(C)
call mma_deallocate(U)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! here we do the same stuff but without prediagonalization

B(1:nPoints) = One
B(nPoints+1:) = Zero
A(:,:) = full_r(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Now form A_1 B (to be used for the computation of the dispersion)

#ifdef _DPOSV_
call DPOSV_('U',m_t,1,A,m_t,B,m_t,INFO)
#else
call DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO)
#endif
if (INFO /= 0) then
  write(u6,*) 'kriging_model: INFO.ne.0'
  write(u6,*) 'kriging_model: INFO=',INFO
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'Info=',Info
call RecPrt('PSI^{-1}',' ',A,m_t,m_t)
call RecPrt('X=PSI^{-1}f',' ',B,1,m_t)
call RecPrt('rones',' ',rones,1,m_t)
write(u6,*) 'nPoints=',nPoints
#endif

! In case of prediagonalization backtransform to the original basis

#ifdef _PREDIAG_
! Call RecPrt('(U^TAU)^{-1}U^TB',' ',B,1,m_t)
D(:) = B(:)
B(:) = Zero
call DGEMM_('N','N',m_t,1,m_t,One,UBIG,m_t,D,m_t,Zero,B,m_t)
!call RecPrt('B',' ',B,1,m_t)
call mma_deallocate(D)
call mma_deAllocate(UBIG)
#endif

rones(:) = B(:)     ! Move result over to storage for later use, R^-1F

! Compute contribution to the expression for the liklihood function.
!
! Now A contains the factors L and U from the factorization A = P*L*U as computed by DGESV
! Where L in the lower triangular matrix with 1 in the diagonal and U is the upper
! triangular matrix of A thus the determinant of A is giving by multipling its diagonal

detR = Zero
do i=1,m_t
#ifdef _DPOSV_
  detR = detR+Two*log(A(i,i))
#else
  detR = detR+log(abs(A(i,i)))
#endif
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Now work on the vector with the values and gradients, the generalized y vector
!
! Trend Function (baseline)

if (blaAI) then

  ! Make sure the base line is above any data point

  sb = -huge(sb)
  do i=1,nPoints
    sb = max(sb,y(i)+blavAI)
  end do
else if (mblAI) then
  sb = sbmev
else if (blAI) then
  sb = blvAI
else
  ordinary = .true.

  ! Note dy only containes gradients for the set of points which we are considering
  ! B(:) = [y(:),dy(:)]  original code before PGEK implementation
  B(1:nPoints) = y(1:nPoints)
  do i_eff=1,nInter_eff
    i = Index_PGEK(i_eff)
    is = (i-1)*(nPoints-nD)+1
    ise = nPoints+(i_eff-1)*(nPoints-nD)+1
    ie = is+(nPoints-nD)-1
    iee = ise+(nPoints-nD)-1
    B(ise:iee) = dy(is:ie)
  end do
#ifdef _DEBUGPRINT_
  write(u6,*) DDot_(m_t,rones,1,B,1),DDot_(nPoints,rones,1,[One],0)
#endif
  ! sbO:  FR^-1y/(FR^-1F)
  sbO = DDot_(m_t,rones,1,B,1)/DDot_(nPoints,rones,1,[One],0)
  sb = sbO
end if

!**********************************************************************
!
! form the actual value vector (y-b_0F)

!B(1:m_t) = [y(1:nPoints)-sb,dy(1:nInter*(nPoints-nD))]
B(1:nPoints) = y(1:nPoints)-sb
do i_eff=1,nInter_eff
  i = Index_PGEK(i_eff)
  is = +(i-1)*(nPoints-nD)+1
  ise = nPoints+(i_eff-1)*(nPoints-nD)+1
  ie = is+(nPoints-nD)-1
  iee = ise+(nPoints-nD)-1
  B(ise:iee) = dy(is:ie)
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'sb,ln(det|PSI|)=',sb,detR
call RecPrt('[y-sb,dy]','(12(2x,E9.3))',B,1,m_t)
#endif

!#undef _PREDIAG_
#ifdef _PREDIAG_

! Diagonalize the energy block of Psi

call mma_allocate(U,nPoints,nPoints,label='U')
U(:,:) = Zero
do i=1,nPoints
  U(i,i) = One
end do
call mma_allocate(HTri,nPoints*(nPoints+2)/2,label='HTri')
do i=1,nPoints
  do j=1,i
    ij = i*(i-1)/2+j
    HTri(ij) = Full_R(i,j)
  end do
end do
call nidiag_new(HTri,U,nPoints,nPoints,0)
call Jacord(HTri,U,nPoints,nPoints)
! Standardized phase factor
do i=1,nPoints
  Temp = DDot_(nPoints,[One],0,U(1,i),1)
  U(1:nPoints,i) = U(1:nPoints,i)*sign(One,Temp)
end do
#ifdef _DEBUGPRINT_
call RecPrt('U',' ',U,nPoints,nPoints)
call TriPrt('HTri',' ',HTri,nPoints)
#endif

! Construct a transformation which will transform the
! covariance  matrix to this new basis.

call mma_Allocate(UBIG,m_t,m_t,label='UBig')
UBIG(:,:) = Zero
UBIG(1:nPoints,1:nPoints) = U(:,:)
do i=nPoints+1,m_t
  UBIG(i,i) = One
end do
! Call RecPrt('UBIG',' ',UBig,m_t,m_t)
!
! Transform the covariance matrix to the new basis.

call mma_Allocate(C,m_t,m_t,label='C')
C(:,:) = Zero
call dgemm_('N','N',m_t,m_t,m_t,One,Full_R,m_t,UBIG,m_t,Zero,C,m_t)
call dgemm_('T','N',m_t,m_t,m_t,One,UBIG,m_t,C,m_t,Zero,A,m_t)

! Cleanup the first block - should be prefectly diagonal.

A(1:nPoints,1:nPoints) = Zero
do i=1,nPoints
  A(i,i) = HTri(i*(i+1)/2)
end do
#ifdef _DEBUGPRINT_
call RecPrt('U^TAU',' ',A,m_t,m_t)
#endif

! transform Kv to the new basis (y-F)

call mma_allocate(D,m_t,label='D')
D(:) = B(:)
Kv = Zero
call dgemm_('T','N',m_t,1,m_t,One,UBIG,m_t,D,m_t,Zero,Kv,m_t)
#ifdef _DEBUGPRINT_
call RecPrt('U^TKv',' ',Kv,1,m_t)
#endif

call mma_deallocate(HTri)
call mma_deAllocate(C)
call mma_deallocate(U)
#else
Kv(:) = B(:)
A(:,:) = full_r(:,:)
#endif

#ifdef _DEBUGPRINT_
call RecPrt('A',' ',A,m_t,m_t)
call RecPrt('Kv',' ',Kv,1,m_t)
#endif
#ifdef _DPOSV_
call DPOSV_('U',m_t,1,A,m_t,Kv,m_t,INFO)
#else
call DGESV_(m_t,1,A,m_t,IPIV,Kv,m_t,INFO)
#endif
#ifdef _PREDIAG_
#ifdef _DEBUGPRINT_
call RecPrt('(U^TAU)^{-1}U^TKv',' ',Kv,1,m_t)
#endif
D(:) = Kv(:)
Kv(:) = Zero
call DGEMM_('N','N',m_t,1,m_t,One,UBIG,m_t,D,m_t,Zero,Kv,m_t)
!call RecPrt('Kv',' ',Kv,1,m_t)
call mma_deallocate(D)
call mma_deAllocate(UBIG)
#endif

!Likelihood function
variance = dot_product(B,Kv)/real(m_t,kind=wp)
lh = variance*exp(detR/real(m_t,kind=wp))
#ifdef _DEBUGPRINT_
write(u6,*) 'Variance=',Variance
write(u6,*) 'Info=',Info
call RecPrt('X=A^{-1}Kv','(5(E15.7,2X))',Kv,1,m_t)
#endif

call mma_Deallocate(B)
call mma_Deallocate(A)
call mma_Deallocate(IPIV)

end subroutine kriging_model
