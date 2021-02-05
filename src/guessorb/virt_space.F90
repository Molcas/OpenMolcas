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

!***********************************************************************
!     The generation of starting orbitals suffers from the fact that   *
!     the virtual orbitals are not well defined. This routine is       *
!     supposed to generate a set of well-defined virtual orbitals from *
!     a set of well-defined occupied orbitals.                         *
!***********************************************************************
subroutine Virt_Space(C_Occ,C_Virt,Ovrlp,nBas,nOcc,nVirt)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas, nOcc, nVirt
real(kind=wp), intent(in) :: C_Occ(nBas,nOcc), Ovrlp(nBas*(nBas+1)/2)
real(kind=wp), intent(inout) :: C_Virt(nBas,nVirt)
integer(kind=iwp) :: iBas, jBas, iOcc, kBas, iVirt, lBas, mVirt
logical(kind=iwp) :: Polished
real(kind=wp) :: tmp
real(kind=wp), allocatable :: P(:,:), Ovrlp_Sq(:,:), EVe(:,:), C_tmp(:,:), PNew(:), EVa(:)
real(kind=wp), parameter :: thr = 1.0e-14_wp

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, j, ij

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

if (nVirt == 0) call Abend()

! Compute S^{1/2}

call mma_allocate(Ovrlp_Sq,nBas,nBas,Label='Ovrlp_Sq')
call mma_allocate(EVa,nBas*(nBas+1)/2,Label='EVa')
call mma_allocate(EVe,nBas,nBas,Label='EVe')
call FZero(EVe,nBas**2)
call DCopy_(nBas,[One],0,EVe,nBas+1)
call DCopy_(nBas*(nBas+1)/2,Ovrlp,1,EVa,1)
call NIdiag(EVa,EVe,nBas,nBas,0)

do iBas=2,nBas
  EVa(iBas) = EVa(iBas*(iBas+1)/2)
end do

do iBas=1,nBas
  do jBas=1,nBas
    tmp = Zero
    do kBas=1,nBas
      tmp = tmp+EVe(iBas,kBas)*sqrt(EVa(kBas))*EVe(jBas,kBas)
    end do
    Ovrlp_Sq(iBas,jBas) = tmp
  end do
end do

! Now transform the the CMOs of the occupied orbitals to an
! orthonormal basis.

call mma_allocate(C_tmp,nBas,nOcc,Label='C_tmp')
call FZero(C_tmp,nBas*nOcc)
call DGEMM_('N','N',nBas,nOcc,nBas,One,Ovrlp_Sq,nBas,C_Occ,nBas,Zero,C_tmp,nBas)

call mma_allocate(P,nBas,nBas,Label='P')
call mma_allocate(PNew,nBas,Label='PNew')

! Form the P matrix = 1 - |C(occ)><C(occ)|

do iBas=1,nBas
  do jBas=1,nBas
    tmp = Zero
    if (iBas == jBas) tmp = One
    do iOcc=1,nOcc
      tmp = tmp-C_tmp(iBas,iOcc)*C_tmp(jBas,iOcc)
    end do
    P(iBas,jBas) = tmp
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('P-mat',' ',P,nBas,nBas)
do iBas=1,nBas
  write(u6,*) 'iBas,P(iBas,iBas)=',iBas,P(iBas,iBas)
end do
#endif

! Now, use Cholesky Decomposition to reduce
! the size of the P-matrix from a nBas x nBas size
! down to something that is nBas x nVirt

mVirt = 0
kBas_1_nBas: do kBas=1,nBas

  tmp = Zero
  lBas = 0
  do iBas=1,nBas
    if (P(iBas,iBas) > tmp) then
      tmp = P(iBas,iBas)
      lBas = iBas
    end if
  end do

  ! Pick up a vector and project it against the previous

  call DCopy_(nBas,P(1,lBas),1,PNew,1)

  ! Normalize PNew

  tmp = Zero
  do iBas=1,nBas
    tmp = tmp+PNew(iBas)**2
  end do

  ! Skip if this is a null vector.

  if (tmp < thr) cycle

  tmp = One/sqrt(tmp)
  call DScal_(nBas,tmp,PNew,1)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'New Trial vector'
  write(u6,*) '================'
  write(u6,*) 'kBas,lBas,tmp=',kBas,lBas,tmp
  call RecPrt('Normalized PNew',' ',PNew,nBas,1)
# endif

  !Polished = .false.
  Polished = .true.
  do_Polished: do

    do iOcc=1,nOcc

      ! From the trial vector eliminate the occupied space

      tmp = Zero
      do iBas=1,nBas
        tmp = tmp+PNew(iBas)*C_tmp(iBas,iOcc)
      end do
#     ifdef _DEBUGPRINT_
      write(u6,*) 'iOcc,tmp=',iOcc,tmp
#     endif
      ! Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)

      call DaXpY_(nBas,-tmp,C_Occ(1,iOcc),1,PNew,1)
    end do
    do iVirt=1,mVirt

      ! From the trial vector eliminate parts which already expressed.

      tmp = Zero
      do iBas=1,nBas
        tmp = tmp+PNew(iBas)*C_Virt(iBas,iVirt)
      end do
#     ifdef _DEBUGPRINT_
      write(u6,*) 'iVirt,tmp=',iVirt,tmp
#     endif

      ! Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)

      call DaXpY_(nBas,-tmp,C_Virt(1,iVirt),1,PNew,1)
    end do

    ! Test that it is not a null vector!

    tmp = Zero
    do iBas=1,nBas
      tmp = tmp+PNew(iBas)**2
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Norm after projection:',tmp
    write(u6,*)
#   endif
    if (tmp > thr) then
      tmp = One/sqrt(tmp)
      call DScal_(nBas,One/sqrt(tmp),PNew,1)
      if (.not. Polished) then
        Polished = .true.
        cycle do_Polished
      end if

      if (tmp > thr) then
        if (mVirt+1 > nVirt) then
          write(u6,*) 'mVirt.gt.nVirt'
          write(u6,*) 'mVirt=',mVirt
          write(u6,*) 'nVirt=',nVirt
          call Abend()
        end if
        mVirt = mVirt+1
        call DCopy_(nBas,PNew,1,C_Virt(1,mVirt),1)

        ! Update the P-matrix.

        do iBas=1,nBas
          do jBas=1,nBas
            P(iBas,jBas) = P(iBas,jBas)-PNew(iBas)*PNew(jBas)
          end do
        end do
      end if
    end if
    if (Polished) exit do_Polished
  end do do_Polished

  if (mVirt == nVirt) exit kBas_1_nBas

end do kBas_1_nBas

call mma_deallocate(C_tmp)

if (mVirt /= nVirt) then
  write(u6,*) 'mVirt.ne.nVirt'
  write(u6,*) 'mVirt,nVirt=',mVirt,nVirt
  call Abend()
end if

call mma_deallocate(P)
call mma_deallocate(PNew)

! Form S^(-1/2)

do iBas=1,nBas
  do jBas=1,nBas
    tmp = Zero
    do kBas=1,nBas
      tmp = tmp+(EVe(iBas,kBas)*EVe(jBas,kBas))/sqrt(EVa(kBas))
    end do
    Ovrlp_Sq(iBas,jBas) = tmp
  end do
end do
call mma_deallocate(EVe)
call mma_deallocate(EVa)

call mma_allocate(C_tmp,nBas,nVirt,Label='C_tmp')
call DCopy_(nBas*nVirt,C_Virt,1,C_tmp,1)
call DGEMM_('N','N',nBas,nVirt,nBas,One,Ovrlp_Sq,nBas,C_tmp,nBas,Zero,C_Virt,nBas)
call mma_deallocate(C_tmp)
call mma_deallocate(Ovrlp_Sq)

#ifdef _DEBUGPRINT_
call RecPrt('C_Virt(New)',' ',C_Virt,nBas,nVirt)
#endif

return

#ifdef _DEBUGPRINT_
contains

function iTri(i,j)
  integer(kind=iwp) :: iTri
  integer(kind=iwp), intent(in) :: i, j
  iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
end function iTri
#endif

end subroutine Virt_Space
