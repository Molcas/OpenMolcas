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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module ROMP2_Semicanonical
use Definitions, only: wp, iwp
use MBPT2_Global, only: MBPT2_All_Finite
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
implicit none
private
public :: ROMP2_Semi
contains

subroutine ROMP2_Semi(C,F,n,no,nFrozen,eps,E1,ierr)
  integer(kind=iwp), intent(in) :: n,no,nFrozen
  real(kind=wp), intent(inout) :: C(n,n)
  real(kind=wp), intent(in) :: F(n,n)
  real(kind=wp), intent(out) :: eps(n),E1
  integer(kind=iwp), intent(out) :: ierr
  real(kind=wp), allocatable :: M(:,:),U(:,:),W(:),Cnew(:,:),T(:,:)
  real(kind=wp) :: den
  integer(kind=iwp) :: k,first,last,nBlock,info,i,a
  ierr=0
  E1=Zero
  eps(:)=Zero
  ! Frozen-core convention: diagonalize the FULL occupied block first,
  ! then exclude its nFrozen lowest eigenvectors from correlation.
  ! Alpha and beta frozen subspaces may differ after semicanonicalization.
  ! This is not a fixed common-spatial-core projection convention.
  if ((n<1).or.(no<0).or.(no>n).or.(nFrozen<0).or.(nFrozen>no)) then
    ierr=4
    return
  end if
  if ((.not.MBPT2_All_Finite(C)).or.(.not.MBPT2_All_Finite(F))) then
    ierr=3
    return
  end if
  call mma_allocate(M,n,n,label='ROMP2 MO Fock')
  call mma_allocate(Cnew,n,n,label='ROMP2 rotated orbitals')
  call mma_allocate(T,n,n,label='ROMP2 Fock work')
  call DGEMM_('N','N',n,n,n,One,F,n,C,n,Zero,T,n)
  call DGEMM_('T','N',n,n,n,One,C,n,T,n,Zero,M,n)
  ! Diagonalize each spin's occupied and virtual blocks separately.
  ! No occupied/virtual mixing is permitted: the determinant is preserved.
  do k=1,2
    if (k==1) then
      first=1
      last=no
    else
      first=no+1
      last=n
    end if
    nBlock=last-first+1
    if (nBlock==0) cycle
    call mma_allocate(U,nBlock,nBlock,label='ROMP2 eigenvectors')
    call mma_allocate(W,max(1,3*nBlock),label='ROMP2 eigensolver work')
    U(:,:)=M(first:last,first:last)
    call DSYEV_('V','U',nBlock,U,nBlock,eps(first),W,size(W),info)
    if (info/=0) then
      ierr=1
      call mma_deallocate(U)
      call mma_deallocate(W)
      exit
    end if
    call DGEMM_('N','N',n,nBlock,nBlock,One,C(1,first),n,U,nBlock,Zero,Cnew(1,first),n)
    call mma_deallocate(U)
    call mma_deallocate(W)
  end do
  if (ierr==0) then
    C(:,:)=Cnew
    call DGEMM_('N','N',n,n,n,One,F,n,C,n,Zero,T,n)
    call DGEMM_('T','N',n,n,n,One,C,n,T,n,Zero,M,n)
    ! Core electrons remain in the ROHF density/Fock/reference energy.
    ! Only their correlation excitations are excluded, including singles.
    do i=nFrozen+1,no
      do a=no+1,n
        den=eps(i)-eps(a)
        if ((.not. MBPT2_All_Finite(den)) .or. (abs(den)<=1.0e-12_wp)) then
          ierr=2
          exit
        end if
        E1=E1+M(i,a)**2/den
      end do
      if (ierr/=0) exit
    end do
    if (.not. MBPT2_All_Finite(eps)) ierr=3
    if (.not. MBPT2_All_Finite(E1)) ierr=3
  end if
  call mma_deallocate(M)
  call mma_deallocate(Cnew)
  call mma_deallocate(T)
end subroutine ROMP2_Semi
end module ROMP2_Semicanonical
