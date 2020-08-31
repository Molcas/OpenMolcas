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
SUBROUTINE predict(gh,iter,nInter)
  use kriging_mod
#include "stdalloc.fh"
  real*8 tsum
  integer INFO
  integer i,j,iter,nInter,gh ! ipiv the pivot indices that define the permutation matrix
  real*8, Allocatable :: B(:), A(:,:)
  integer, Allocatable :: IPIV(:)
!
  Call mma_allocate(B,m_t,label="B")
!
  do j=1,npx
    if (gh.eq.0) then
      !A contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
      Call mma_allocate(A,m_t,m_t,label="A")
      Call mma_allocate(IPIV,m_t,label="IPIV")
      ! calculations of Energy and dispersion
      A(:,:) = full_R
      B(:) = cv(:,j,1,1)
      pred(j) = sb + dot_product(B,Kv)
      CALL DGESV_(m_t, 1,A,m_t,IPIV,B,m_t,INFO )
      var(j) = 1d0 - dot_product(B,CV(:,j,1,1))
      if (ordinary) Then
        tsum = sum(rones(1:iter))
        B(:) = cv(:,j,1,1)
        var(j)=max(var(j)+(1d0-dot_product(B,rones))**2/tsum,0d0)
      end if
      sigma(j)=sqrt(var(j)*variance)
      Call mma_deallocate(A)
      Call mma_deallocate(IPIV)
    else if (gh.eq.1) then
      do k=1,nInter
        B(:) = cv(:,j,k,1)
        gpred(j,k) = dot_product(B,Kv)
      enddo
    else if (gh.eq.2) then
      do k=1,nInter
        do i=k,nInter
          B(:) = cv(:,j,i,k)
          hpred(j,k,i) = dot_product(B, Kv)
          if (i.ne.k) hpred(j,i,k) = hpred(j,k,i)
        enddo
      enddo
    endif
  enddo
!
  Call mma_deallocate(B)
!
END SUBROUTINE predict
