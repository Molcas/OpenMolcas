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
SUBROUTINE predict(gh)
  use kriging_mod
#include "stdalloc.fh"
  real*8 tsum
  integer INFO
  integer i,j,gh ! ipiv the pivot indices that define the permutation matrix
  real*8, Allocatable :: B(:), A(:,:)
  integer, Allocatable :: IPIV(:)
!
  Call mma_allocate(B,m_t,label="B")
!
nInter=nInter_save
    if (gh.eq.0) then

      !A contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
      Call mma_allocate(A,m_t,m_t,label="A")
      Call mma_allocate(IPIV,m_t,label="IPIV")
      ! calculations of Energy and dispersion
      A(:,:) = full_R
      B(:) = cv(:,1,1,1)
      pred(1) = sb + dot_product(B,Kv)
      CALL DGESV_(m_t, 1,A,m_t,IPIV,B,m_t,INFO )
      var(1) = 1d0 - dot_product(B,CV(:,1,1,1))

      if (ordinary) Then
        tsum = sum(rones(1:m_t))
        B(:) = cv(:,1,1,1)
        var(1)=max(var(1)+(1d0-dot_product(B,rones))**2/tsum,0d0)
      end if

      sigma(1)=sqrt(var(1)*variance)
      Call mma_deallocate(A)
      Call mma_deallocate(IPIV)

    else if (gh.eq.1) then

      do k=1,nInter
        B(:) = cv(:,1,k,1)
        gpred(1,k) = dot_product(B,Kv)
      enddo

    else if (gh.eq.2) then

      do k=1,nInter
        do i=k,nInter
          B(:) = cv(:,1,i,k)
          hpred(1,k,i) = dot_product(B, Kv)
          if (i.ne.k) hpred(1,i,k) = hpred(1,k,i)
        enddo
      enddo

    endif
!
  Call mma_deallocate(B)
!
END SUBROUTINE predict
