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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Virt_Space(C_Occ,C_Virt,Ovrlp,nBas,nOcc,nVirt)
!***********************************************************************
!     The generation of starting orbitals suffers from the fact that   *
!     the virtual orbitals are not well defined. This routine is       *
!     supposed to generate a set of well-defined virtual orbitals from *
!     a set of well-defined occupied orbitals.                         *
!***********************************************************************

use Index_Functions, only: nTri_Elem
#ifdef _DEBUGPRINT_
use Index_Functions, only: iTri
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas, nOcc, nVirt
real(kind=wp), intent(in) :: C_Occ(nBas,nOcc), Ovrlp(nTri_Elem(nBas))
real(kind=wp), intent(inout) :: C_Virt(nBas,nVirt)
integer(kind=iwp) :: iBas, iOcc, iVirt, jBas, kBas, lBas, mVirt
real(kind=wp) :: tmp
real(kind=wp), allocatable :: C_tmp(:,:), EVa(:), EVe(:,:), Ovrlp_Sq(:,:), P(:,:), PNew(:)
real(kind=wp), parameter :: thr = 1.0e-14_wp
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ij

write(u6,*) 'nBas,nOcc,nVirt=',nBas,nOcc,nVirt
call RecPrt('C_Occ',' ',C_Occ,nBas,nOcc)
call RecPrt('C_Virt',' ',C_Virt,nBas,nVirt)
call TriPrt('Ovrlp',' ',Ovrlp,nBas)
do iOcc=1,nOcc
  tmp = Zero
  do iBas=1,nBas
    do jBas=1,nBas
      ij = iTri(iBas,jBas)
      tmp = tmp+C_Occ(iBas,iOcc)*Ovrlp(ij)*C_Occ(jBas,iOcc)
    end do
  end do
  write(u6,*) 'iOcc,tmp=',iOcc,tmp
end do
do iVirt=1,nVirt
  tmp = Zero
  do iBas=1,nBas
    do jBas=1,nBas
      ij = iTri(iBas,jBas)
      tmp = tmp+C_Virt(iBas,iVirt)*Ovrlp(ij)*C_Virt(jBas,iVirt)
    end do
  end do
  write(u6,*) 'iVirt,tmp=',iOcc,tmp
end do
#endif

! Compute S^{1/2}

call mma_allocate(Ovrlp_Sq,nBas,nBas,Label='Ovrlp_Sq')
call mma_allocate(EVa,nTri_Elem(nBas),Label='EVa')
call mma_allocate(EVe,nBas,nBas,Label='EVe')
call unitmat(EVe,nBas)
EVa(:) = Ovrlp(1:nTri_Elem(nBas))
call NIdiag(EVa,EVe,nBas,nBas)

do iBas=2,nBas
  EVa(iBas) = EVa(nTri_Elem(iBas))
end do

do iBas=1,nBas
  do jBas=1,nBas
    Ovrlp_Sq(iBas,jBas) = sum(EVe(iBas,:)*sqrt(EVa(1:nBas))*EVe(jBas,:))
  end do
end do

! Now transform the CMOs of the occupied orbitals to an orthonormal basis.

call mma_allocate(C_tmp,nBas,nOcc,Label='C_tmp')
call DGEMM_('N','N',nBas,nOcc,nBas,One,Ovrlp_Sq,nBas,C_Occ,nBas,Zero,C_tmp,nBas)

call mma_allocate(P,nBas,nBas,Label='P')
call mma_allocate(PNew,nBas,Label='PNew')

! Form the P matrix = 1 - |C(occ)><C(occ)|

call unitmat(P,nBas)
call DGEMM_('N','T',nBas,nBas,nOcc,-One,C_tmp,nBas,C_tmp,nBas,One,P,nBas)
#ifdef _DEBUGPRINT_
call RecPrt('C_tmp',' ',C_tmp,nBas,nOcc)
call RecPrt('P-mat',' ',P,nBas,nBas)
do iBas=1,nBas
  write(u6,*) 'iBas,P(iBas,iBas)=',iBas,P(iBas,iBas)
end do
#endif

! Now, use Cholesky Decomposition to reduce
! the size of the P-matrix from a nBas x nBas size
! down to something that is nBas x nVirt

mVirt = 0
do kBas=1,nBas

  tmp = Zero
  lBas = 0
  do iBas=1,nBas
    if (P(iBas,iBas) > tmp) then
      tmp = P(iBas,iBas)
      lBas = iBas
    end if
  end do

  ! Pick up a vector and project it against the previous

  PNew(:) = P(:,lBas)

  ! Normalize PNew

  tmp = sum(PNew(:)**2)

  ! Skip if this is a null vector.

  if (tmp < thr) cycle

  PNew(:) = PNew(:)/sqrt(tmp)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'New Trial vector'
  write(u6,*) '================'
  write(u6,*) 'kBas,lBas,tmp=',kBas,lBas,tmp
  call RecPrt('Normalized PNew',' ',PNew,nBas,1)
# endif

  do iOcc=1,nOcc

    ! From the trial vector eliminate the occupied space

    tmp = sum(PNew(:)*C_tmp(:,iOcc))
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iOcc,tmp=',iOcc,tmp
#   endif
    ! Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)

    PNew(:) = PNew(:)-tmp*C_Occ(:,iOcc)
  end do
  do iVirt=1,mVirt

    ! From the trial vector eliminate parts which are already expressed.

    tmp = sum(PNew(:)*C_Virt(:,iVirt))
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iVirt,tmp=',iVirt,tmp
#   endif

    ! Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)

    PNew(:) = PNew(:)-tmp*C_Virt(:,iVirt)
  end do

  ! Test that it is not a null vector!

  tmp = sum(PNew(:)**2)
# ifdef _DEBUGPRINT_
  write(u6,*) 'Norm after projection:',tmp
  write(u6,*)
# endif
  if (tmp > thr) then
    PNew(:) = PNew(:)/sqrt(tmp)

    mVirt = mVirt+1
    C_Virt(:,mVirt) = PNew(:)

    ! Update the P-matrix.

    do jBas=1,nBas
      P(:,jBas) = P(:,jBas)-PNew(:)*PNew(jBas)
    end do
  end if

  if (mVirt == nVirt) exit

end do

call mma_deallocate(C_tmp)
call mma_deallocate(P)
call mma_deallocate(PNew)

if (mVirt /= nVirt) then
  write(u6,*) 'mVirt /= nVirt'
  write(u6,*) 'mVirt,nVirt=',mVirt,nVirt
  call Abend()
end if

! Form S^(-1/2)

do iBas=1,nBas
  do jBas=1,nBas
    Ovrlp_Sq(iBas,jBas) = sum(EVe(iBas,:)*EVe(jBas,:)/sqrt(EVa(1:nBas)))
  end do
end do
call mma_deallocate(EVe)
call mma_deallocate(EVa)

call mma_allocate(C_tmp,nBas,nVirt,Label='C_tmp')
C_tmp(:,:) = C_Virt(:,:)
call DGEMM_('N','N',nBas,nVirt,nBas,One,Ovrlp_Sq,nBas,C_tmp,nBas,Zero,C_Virt,nBas)
call mma_deallocate(C_tmp)
call mma_deallocate(Ovrlp_Sq)

#ifdef _DEBUGPRINT_
call RecPrt('C_Virt(New)',' ',C_Virt,nBas,nVirt)
#endif

end subroutine Virt_Space
