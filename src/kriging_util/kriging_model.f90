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
SUBROUTINE kriging_model(nPoints)
  use kriging_mod
  Implicit None
#include "stdalloc.fh"
  External DDOt_
  Real*8 DDOt_
  real*8, Allocatable:: B(:), A(:,:)
#define _PREDIAG_
#ifdef _PREDIAG_
  real*8, Allocatable:: U(:,:), HTri(:), UBIG(:,:), C(:,:), D(:)
  real*8 temp
  integer j,ij
#endif
  Integer, Allocatable:: IPIV(:)
  Integer i,INFO,nPoints ! ipiv the pivot indices that define the permutation matrix
!
  Call mma_Allocate(B,m_t,Label="B")
  Call mma_Allocate(A,m_t,m_t,Label="A")
  Call mma_Allocate(IPIV,m_t,Label="IPIV")
!
! Initiate B according to Eq. (6) of ref.
  B=0.0D0
  B(1:nPoints)=1.0D0
!
! Initiate A according to Eq. (2) of ref.
!
!#define _DEBUG_
#ifdef _DEBUG_
  Call RecPrt('f',' ',B,1,m_t)
#endif
!
#ifdef _PREDIAG_
  Call mma_allocate(U,nPoints,nPoints,Label='U')
  U = 0.0D0
  forall (i=1:nPoints) U(i,i)=1.0D0
  Call mma_allocate(HTri,nPoints*(nPoints+2)/2,Label='HTri')
  Do i = 1, nPoints
    Do j = 1, i
      ij=i*(i-1)/2 + j
      HTri(ij)=Full_R(i,j)
    End Do
  End Do
#ifdef _DEBUG_
  Call RecPrt('U',' ',U,nPoints,nPoints)
#endif
  Call nidiag_new(HTri,U,nPoints,nPoints,0)
  Call Jacord    (HTri,U,nPoints,nPoints)
  Do i = 1, nPoints
    Temp=DDot_(nPoints,[1.0D0],0,U(1,i),1)
    U(1:nPoints,i)= U(1:nPoints,i) * Sign(1.0D0,Temp)
  End Do
#ifdef _DEBUG_
  Call RecPrt('U',' ',U,nPoints,nPoints)
  Call TriPrt('HTri',' ',HTri,nPoints)
#endif
!
  Call mma_Allocate(UBIG,m_t,m_t,Label='UBig')
  UBIG=0.0D0
  forall (i=1:m_t) UBIG(i,i)=1.0D0
  UBIG(1:nPoints,1:nPoints)=U
!           Call RecPrt('UBIG',' ',UBig,m_t,m_t)
!
  Call mma_Allocate(C,m_t,m_t,Label='C')
  C=0.0D0
  Call dgemm_('N','N',m_t,m_t,m_t,    &
              1.0D0,Full_R,m_t,       &
                    UBIG,m_t,         &
              0.0D0,C,m_t)
  Call dgemm_('T','N',m_t,m_t,m_t,    &
              1.0D0,UBIG,m_t,         &
                    C,m_t,            &
              0.0D0,A,m_t)
  A(1:nPoints,1:nPoints)=0.0D0
  forall (i=1:nPoints) A(i,i)=HTri(i*(i+1)/2)
#ifdef _DEBUG_
  Call RecPrt('U^TAU',' ',A,m_t,m_t)
#endif
  Call mma_allocate(D,m_t,Label='D')
  D=0.0D0
  D(1:nPoints)=1.0D0
  B=0.0D0
  Call dgemm_('T','N',m_t,1,m_t,      &
              1.0D0,UBIG,m_t,         &
                    D,m_t,            &
              0.0D0,B,m_t)
#ifdef _DEBUG_
  Call RecPrt('U^TB',' ',B,1,m_t)
#endif
!
  Call mma_deallocate(HTri)
  Call mma_deAllocate(C)
  Call mma_deallocate(U)
#else
  B=0.0D0
  B(1:nPoints)=1.0D0
  A(:,:) = full_r(:,:)
#endif
!
#define _DPOSV_
#ifdef _DPOSV_
  CALL DPOSV_('U',m_t,1,A,m_t,B,m_t,INFO )
#else
  CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
#endif
!
#ifdef _PREDIAG_
! Call RecPrt('(U^TAU)^{-1}U^TB',' ',B,1,m_t)
  D(:)=B
  B=0.0D0
  Call DGEMM_('N','N',m_t,1,m_t,      &
              1.0D0,UBIG,m_t,         &
                    D,m_t,            &
              0.0D0,B,m_t)
! Call RecPrt('B',' ',B,1,m_t)
  Call mma_deallocate(D)
  Call mma_deAllocate(UBIG)
#endif
  rones(:)=B(:)
  ! Write (6,*) 'A old ',A
  !----------------------------
  If (INFO.ne.0) Then
    Write (6,*) 'k: INFO.ne.0'
    Write (6,*) 'k: INFO=',INFO
    Call Abend()
  End If
#ifdef _DEBUG_
  Write (6,*) 'Info=',Info
  Call RecPrt('PSI^{-1}',' ',A,m_t,m_t)
  Call RecPrt('X=PSI^{-1}f',' ',B,1,m_t)
  Call RecPrt('rones',' ',rones,1,m_t)
  Write (6,*) 'nPoints=',nPoints
#endif
!
! Now A contains the factors L and U from the factorization A = P*L*U as computed by DGESV
! Where L in the lower triangular matrix with 1 in the diagonal and U is the upper
! triangular matrix of A thus the determinant of A is giving by multipling its diagonal
!
  detR = 0.0d0
  do i=1,m_t
#ifdef _DPOSV_
    detR = detR + 2.0D0*log(A(i,i))
#else
    detR = detR + log(abs(A(i,i)))
#endif
  enddo
!
!Trend Function (baseline)
  ordinary = .False.
  if (blaAI) then
!
!   Make sure the base line is above any data point
!
    sb = -1.0D99
    Do i = 1, nPoints
      sb = Max( sb, y(i) + blavAI )
    End Do
  else if (mblAI) then
    sb = sbmev
  else if (blAI) Then
    sb = blvAI
  else
    ordinary = .True.
    sbO = dot_product(y,B(1:nPoints))/sum(B(1:nPoints))
    B(:) = [y,dy]
#ifdef _DEBUG_
    Write (6,*) DDot_(m_t,rones,1,B,1) , DDot_(nPoints,rones,1,[1.0D0],0)
#endif
    sbO = DDot_(m_t,rones,1,B,1) / DDot_(nPoints,rones,1,[1.0D0],0)
    sb = sbO
  endif
!
!**********************************************************************
!
  B(:) = [y-sb,dy]
#ifdef _DEBUG_
  Write (6,*) 'sb,ln(det|PSI|)=',sb,detR
  Call RecPrt('[y-sb,dy]',' ',B,1,m_t)
#endif
!#undef _PREDIAG_
#ifdef _PREDIAG_
!
! Diagonalize the energy block of Psi
!
  Call mma_allocate(U,nPoints,nPoints,Label='U')
  U = 0.0D0
  forall (i=1:nPoints) U(i,i)=1.0D0
  Call mma_allocate(HTri,nPoints*(nPoints+2)/2,Label='HTri')
  Do i = 1, nPoints
    Do j = 1, i
      ij=i*(i-1)/2 + j
      HTri(ij)=Full_R(i,j)
    End Do
  End Do
  Call nidiag_new(HTri,U,nPoints,nPoints,0)
  Call Jacord    (HTri,U,nPoints,nPoints)
  Do i = 1, nPoints
    Temp=DDot_(nPoints,[1.0D0],0,U(1,i),1)
    U(1:nPoints,i)= U(1:nPoints,i) * Sign(1.0D0,Temp)
  End Do
#ifdef _DEBUG_
  Call RecPrt('U',' ',U,nPoints,nPoints)
  Call TriPrt('HTri',' ',HTri,nPoints)
#endif
!
! Construct a transformation which will transform the
! correlation matrix to this new basis.
!
  Call mma_Allocate(UBIG,m_t,m_t,Label='UBig')
  UBIG=0.0D0
  forall (i=1:m_t) UBIG(i,i)=1.0D0
  UBIG(1:nPoints,1:nPoints)=U
! Call RecPrt('UBIG',' ',UBig,m_t,m_t)
!
! Transform the correlation matrix to the new basis.
!
  Call mma_Allocate(C,m_t,m_t,Label='C')
  C=0.0D0
  Call dgemm_('N','N',m_t,m_t,m_t,    &
              1.0D0,Full_R,m_t,       &
                    UBIG,m_t,         &
              0.0D0,C,m_t)
  Call dgemm_('T','N',m_t,m_t,m_t,    &
              1.0D0,UBIG,m_t,         &
                    C,m_t,            &
              0.0D0,A,m_t)
!
! Cleanup the first block - should be prefectly diagonal.
!
  A(1:nPoints,1:nPoints)=0.0D0
  forall (i=1:nPoints) A(i,i)=HTri(i*(i+1)/2)
#ifdef _DEBUG_
  Call RecPrt('U^TAU',' ',A,m_t,m_t)
#endif
!
! transform Kv to the new basis
!
  Call mma_allocate(D,m_t,Label='D')
  D(:)=B
  Kv=0.0D0
  Call dgemm_('T','N',m_t,1,m_t,      &
              1.0D0,UBIG,m_t,         &
                    D,m_t,            &
              0.0D0,Kv,m_t)
#ifdef _DEBUG_
  Call RecPrt('U^TKv',' ',Kv,1,m_t)
#endif
!
  Call mma_deallocate(HTri)
  Call mma_deAllocate(C)
  Call mma_deallocate(U)
#else
  Kv(:)=B
  A(:,:)=full_r
#endif
!
#ifdef _DEBUG_
  Call RecPrt('A',' ',A,m_t,m_t)
  Call RecPrt('Kv',' ',Kv,1,m_t)
#endif
#ifdef _DPOSV_
  CALL DPOSV_('U',m_t,1,A,m_t,Kv,m_t,INFO)
#else
  CALL DGESV_(m_t,1,A,m_t,IPIV,Kv,m_t,INFO)
#endif
#ifdef _PREDIAG_
#ifdef _DEBUG_
  Call RecPrt('(U^TAU)^{-1}U^TKv',' ',Kv,1,m_t)
#endif
  D(:)=Kv
  Kv=0.0D0
  Call DGEMM_('N','N',m_t,1,m_t,      &
              1.0D0,UBIG,m_t,         &
                    D,m_t,            &
              0.0D0,Kv,m_t)
!           Call RecPrt('Kv',' ',Kv,1,m_t)
  Call mma_deallocate(D)
  Call mma_deAllocate(UBIG)
#endif
!
!Likelihood function
  variance = dot_product(B,Kv)/dble(m_t)
  lh = variance*exp(detR/dble(m_t))
#ifdef _DEBUG_
  Write (6,*) 'Variance=',Variance
  Write (6,*) 'Info=',Info
  Call RecPrt('X=A^{-1}Kv','(5(E15.7,2X))',Kv,1,m_t)
#endif
!
  Call mma_Deallocate(B)
  Call mma_Deallocate(A)
  Call mma_Deallocate(IPIV)
END SUBROUTINE kriging_model
