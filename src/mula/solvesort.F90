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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine SolveSort(A,C,S,D,W1,W2,W0,C1,C2,C0,r01,r02,r00,icre,iann,nMat,nd1,nd2,OccNumMat,nOsc,nDimTot)
!  Purpose:
!    Solve the secular equation AC = SCD.
!
!  Input:
!    A         : Real two dimensional array
!    S         : Real two dimensional array
!    W0,W1,W2  : Real two dimensional arrays - eigenvectors scaled by square root of harmonic frequencies.
!    C0,C1,C2  : Real two dimensional arrays - inverses of W0,W1 and W2.
!    r00,
!    r01,
!    r02       : Real arrays - geometry in internal coordinates.
!    nMat      : Two dimensional integer array
!
!  Output:
!    C         : Real two dimensional array
!    D         : Real array
!    OccNumMat : Real two dimensional array - occupation numbers for the different modes.
!
!  Calls:
!    Jacob
!    Jacord
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Eight, Half, auTocm
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nd1, nd2, icre(0:nd1,nd2), iann(0:nd1,nd2), nMat(0:nd1,nd2), nOsc, nDimTot
real(kind=wp), intent(in) :: A(nDimTot,nDimTot), W1(nOsc,nOsc), W2(nOsc,nOsc), W0(nOsc,nOsc), C1(nOsc,nOsc), C2(nOsc,nOsc), &
                             C0(nOsc,nOsc), r01(nOsc), r02(nOsc), r00(nOsc)
real(kind=wp), intent(out) :: C(nDimTot,nDimTot), D(nDimTot), OccNumMat(0:nDimTot-1,nOsc)
real(kind=wp), intent(inout) :: S(nDimTot,nDimTot)
integer(kind=iwp) :: i, ii, iOsc, j, jOsc, k, l, la, lc, m, n
real(kind=wp) :: det
real(kind=wp), allocatable :: A1(:,:), A2(:,:), alpha(:,:), alpha1(:,:), alpha2(:,:), Asymm(:,:), B1(:,:), B2(:,:), beta(:,:), &
                              C_col(:), d1(:), d2(:), r_temp1(:), r_temp2(:), SC_col(:), Scr(:), T(:,:), temp1(:,:), temp2(:,:), &
                              Tmp1(:,:)
real(kind=wp), external :: Ddot_

! Initialize.
n = nDimTot
call mma_allocate(Scr,n*(n+1)/2,label='Scr')
call mma_allocate(T,n,n,label='T')
call mma_allocate(Tmp1,n,n,label='T')
call mma_allocate(Asymm,n,n,label='Asymm')

! Cholesky decomposition of S.
do i=1,n
  S(i,i) = S(i,i)+1.0e-6_wp
end do
call Cholesky(S,Tmp1,n)
call unitmat(T,n)
call Dool_MULA(Tmp1,n,n,T,n,n,det)

! Make A symmetric and transform it to lower packed storage in Scratch.
call DGEMM_('N','N',n,n,n,One,A,n,T,n,Zero,Tmp1,n)
call DGEMM_('T','N',n,n,n,One,T,n,Tmp1,n,Zero,Asymm,n)
k = 1
do i=1,n
  do j=1,i
    Scr(k) = Asymm(i,j)
    k = k+1
  end do
end do

! Diagonalize Scratch.
call unitmat(Tmp1,n)
call Jacob(Scr,Tmp1,n,n)
call Jacord(Scr,Tmp1,n,n)

! Store the eigenvalues in array D and calculate C.
do i=1,n
  ii = i*(i+1)/2
  D(i) = Scr(ii)
end do

call DGEMM_('N','N',n,n,n,One,T,n,Tmp1,n,Zero,C,n)
call mma_deallocate(Scr)
call mma_deallocate(T)
call mma_deallocate(Tmp1)
call mma_deallocate(Asymm)

! Calculate alpha1, alpha2 and alpha.
call mma_allocate(temp1,nOsc,nOsc,label='temp1')
call mma_allocate(temp2,nOsc,nOsc,label='temp2')

call mma_allocate(alpha,nOsc,nOsc,label='alpha')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')

call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C1,nOsc,C1,nOsc,Zero,alpha1,nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C2,nOsc,C2,nOsc,Zero,alpha2,nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C0,nOsc,C0,nOsc,Zero,alpha,nOsc)

! Calculate beta.
call mma_allocate(beta,nOsc,nOsc,label='beta')
call mma_allocate(r_temp1,nOsc,label='r_temp1')
call mma_allocate(r_temp2,nOsc,label='r_temp2')

temp1(:,:) = alpha1
temp2(:,:) = Two*alpha

call Dool_MULA(temp2,nOsc,nOsc,temp1,nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,alpha2,nOsc,temp1,nOsc,Zero,beta,nOsc)

! Calculate A, B and d matrices.
call mma_allocate(A1,nOsc,nOsc,label='A1')
call mma_allocate(B1,nOsc,nOsc,label='B1')
call mma_allocate(A2,nOsc,nOsc,label='A2')
call mma_allocate(B2,nOsc,nOsc,label='B2')
call mma_allocate(d1,nOsc,label='d1')
call mma_allocate(d2,nOsc,label='d2')
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C1,nOsc,W0,nOsc,Zero,A1,nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W1,nOsc,C0,nOsc,Zero,temp1,nOsc)
B1(:,:) = A1-temp1
r_temp1(:) = r01-r00
call DGEMM_('N','N',nOsc,1,nOsc,One,beta,nOsc,r_temp1,nOsc,Zero,r_temp2,nOsc)
call DGEMM_('T','N',nOsc,1,nOsc,sqrt(Eight),W1,nOsc,r_temp2,nOsc,Zero,d1,nOsc)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C2,nOsc,W0,nOsc,Zero,A2,nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W2,nOsc,C0,nOsc,Zero,temp2,nOsc)
B2(:,:) = A2-temp2
r_temp1(:) = r02-r00
call DGEMM_('N','N',nOsc,1,nOsc,One,beta,nOsc,r_temp1,nOsc,Zero,r_temp2,nOsc)
call DGEMM_('T','N',nOsc,1,nOsc,sqrt(Eight),W2,nOsc,r_temp2,nOsc,Zero,d2,nOsc)

call mma_deallocate(temp1)
call mma_deallocate(temp2)

call mma_deallocate(alpha)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)

call mma_deallocate(beta)
call mma_deallocate(r_temp1)
call mma_deallocate(r_temp2)

m = (n/2)-1

call mma_allocate(C_col,[0,n-1],label='C_col')
call mma_allocate(SC_col,[0,n-1],label='SC_col')
OccNumMat(:,:) = Zero
do k=1,n
  do jOsc=1,nOsc
    C_col(:) = Zero
    do l=0,m
      do iOsc=1,nOsc
        lc = iCre(l,iosc)
        la = iAnn(l,iosc)
        if (lc >= 0) C_col(lc) = C_col(lc)+sqrt(real(nmat(lc,iosc),kind=wp))*c(l+1,k)*B1(iosc,josc)
        if (la >= 0) C_col(la) = C_col(la)+sqrt(real(nmat(l,iosc),kind=wp))*c(l+1,k)*A1(iosc,josc)
      end do
      C_col(l) = C_col(l)-c(l+1,k)*d1(josc)
    end do
    do l=0,m
      do iOsc=1,nOsc
        lc = icre(l,iosc)
        la = iann(l,iosc)
        if (lc >= 0) C_col(m+lc+1) = C_col(m+lc+1)+sqrt(real(nmat(lc,iosc),kind=wp))*c(m+l+2,k)*B2(iosc,josc)
        if (la >= 0) C_col(m+la+1) = C_col(m+la+1)+sqrt(real(nmat(l,iosc),kind=wp))*c(m+l+2,k)*A2(iosc,josc)
      end do
      C_col(m+l+1) = C_col(m+l+1)-c(m+l+2,k)*d2(josc)
    end do
    call DGEMM_('N','N',n,1,n,One,S,n,C_col,n,Zero,SC_col,n)
    OccNumMat(k-1,jOsc) = Ddot_(n,C_col,1,SC_col,1)
  end do
end do
write(u6,*)
write(u6,*) 'Frequencies'
write(u6,'(20i6)') (int((D(i+1)-D(1))*auTocm),i=1,m)
write(u6,*)

call mma_deallocate(A1)
call mma_deallocate(B1)
call mma_deallocate(A2)
call mma_deallocate(B2)
call mma_deallocate(d1)
call mma_deallocate(d2)
call mma_deallocate(C_col)
call mma_deallocate(SC_col)

end subroutine SolveSort
