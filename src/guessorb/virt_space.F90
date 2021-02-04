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

implicit none
#include "stdalloc.fh"
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer i, j
integer ij, iTri
#endif
integer nOcc, nVirt, nBas
integer iBas, jBas, iOcc, kBas, iVirt, lBas
integer mVirt
logical Polished
real*8 C_Occ(nBas,nOcc), C_Virt(nBas,nVirt), Ovrlp(nBas*(nBas+1)/2)
real*8 tmp
real*8, dimension(:,:), allocatable :: P, Ovrlp_Sq, EVe, C_tmp
real*8, dimension(:), allocatable :: PNew, EVa
#ifdef _DEBUGPRINT_
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
write(6,*) 'nBas,nOcc,nVirt=',nBas,nOcc,nVirt
call RecPrt('C_Occ',' ',C_Occ,nBas,nOcc)
call RecPrt('C_Virt',' ',C_Virt,nBas,nVirt)
call TriPrt('Ovrlp',' ',Ovrlp,nBas)
do iOcc=1,nOcc
  tmp = 0.0d0
  do iBas=1,nBas
    do jBas=1,nBas
      ij = iTri(iBas,jBas)
      tmp = tmp+C_Occ(iBas,iOcc)*Ovrlp(ij)*C_Occ(jBas,iOcc)
    end do
  end do
  write(6,*) 'iOcc,tmp=',iOcc,tmp
end do
do iVirt=1,nVirt
  tmp = 0.0d0
  do iBas=1,nBas
    do jBas=1,nBas
      ij = iTri(iBas,jBas)
      tmp = tmp+C_Virt(iBas,iVirt)*Ovrlp(ij)*C_Virt(jBas,iVirt)
    end do
  end do
  write(6,*) 'iVirt,tmp=',iOcc,tmp
end do
#endif

if (nVirt == 0) call Abend()

! Compute S^{1/2}

call mma_allocate(Ovrlp_Sq,nBas,nBas,Label='Ovrlp_Sq')
call mma_allocate(EVa,nBas*(nBas+1)/2,Label='EVa')
call mma_allocate(EVe,nBas,nBas,Label='EVe')
call FZero(EVe,nBas**2)
call DCopy_(nBas,[1.0d0],0,EVe,nBas+1)
call DCopy_(nBas*(nBas+1)/2,Ovrlp,1,EVa,1)
call NIdiag(EVa,EVe,nBas,nBas,0)

do iBas=2,nBas
  EVa(iBas) = EVa(iBas*(iBas+1)/2)
end do

do iBas=1,nBas
  do jBas=1,nBas
    tmp = 0.0d0
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
call DGEMM_('N','N',nBas,nOcc,nBas,1.0d0,Ovrlp_Sq,nBas,C_Occ,nBas,0.0d0,C_tmp,nBas)

call mma_allocate(P,nBas,nBas,Label='P')
call mma_allocate(PNew,nBas,Label='PNew')

! Form the P matrix = 1 - |C(occ)><C(occ)|

do iBas=1,nBas
  do jBas=1,nBas
    tmp = 0.0d0
    if (iBas == jBas) tmp = 1.0d0
    do iOcc=1,nOcc
      tmp = tmp-C_tmp(iBas,iOcc)*C_tmp(jBas,iOcc)
    end do
    P(iBas,jBas) = tmp
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('P-mat',' ',P,nBas,nBas)
do iBas=1,nBas
  write(6,*) 'iBas,P(iBas,iBas)=',iBas,P(iBas,iBas)
end do
#endif

! Now, use Cholesky Decomposition to reduce
! the size of the P-matrix from a nBas x nBas size
! down to something that is nBas x nVirt

mVirt = 0
do kBas=1,nBas

  tmp = 0.0d0
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

  tmp = 0.0d0
  do iBas=1,nBas
    tmp = tmp+PNew(iBas)**2
  end do

  ! Skip if this is a null vector.

  if (tmp < 1.0d-14) cycle

  tmp = 1.0d0/sqrt(tmp)
  call DScal_(nBas,tmp,PNew,1)
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'New Trial vector'
  write(6,*) '================'
  write(6,*) 'kBas,lBas,tmp=',kBas,lBas,tmp
  call RecPrt('Normalized PNew',' ',PNew,nBas,1)
# endif

  !Polished = .false.
  Polished = .true.
666 continue

  do iOcc=1,nOcc

    ! From the trial vector eliminate the occupied space

    tmp = 0.0d0
    do iBas=1,nBas
      tmp = tmp+PNew(iBas)*C_tmp(iBas,iOcc)
    end do
#   ifdef _DEBUGPRINT_
    write(6,*) 'iOcc,tmp=',iOcc,tmp
#   endif
    ! Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)

    call DaXpY_(nBas,-tmp,C_Occ(1,iOcc),1,PNew,1)
  end do
  do iVirt=1,mVirt

    ! From the trial vector eliminate parts which already expressed.

    tmp = 0.0d0
    do iBas=1,nBas
      tmp = tmp+PNew(iBas)*C_Virt(iBas,iVirt)
    end do
#   ifdef _DEBUGPRINT_
    write(6,*) 'iVirt,tmp=',iVirt,tmp
#   endif

    ! Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)

    call DaXpY_(nBas,-tmp,C_Virt(1,iVirt),1,PNew,1)
  end do

  ! Test that it is not a null vector!

  tmp = 0.0d0
  do iBas=1,nBas
    tmp = tmp+PNew(iBas)**2
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'Norm after projection:',tmp
  write(6,*)
# endif
  if (tmp > 1.0d-14) then
    tmp = 1.0d0/sqrt(tmp)
    call DScal_(nBas,1.0d0/sqrt(tmp),PNew,1)
    if (.not. Polished) then
      Polished = .true.
      Go To 666
    end if

    if (tmp > 1.0d-14) then
      if (mVirt+1 > nVirt) then
        write(6,*) 'mVirt.gt.nVirt'
        write(6,*) 'mVirt=',mVirt
        write(6,*) 'nVirt=',nVirt
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

  if (mVirt == nVirt) Go To 667

end do

667 continue
call mma_deallocate(C_tmp)

if (mVirt /= nVirt) then
  write(6,*) 'mVirt.ne.nVirt'
  write(6,*) 'mVirt,nVirt=',mVirt,nVirt
  call Abend()
end if

call mma_deallocate(P)
call mma_deallocate(PNew)

! Form S^(-1/2)

do iBas=1,nBas
  do jBas=1,nBas
    tmp = 0.0d0
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
call DGEMM_('N','N',nBas,nVirt,nBas,1.0d0,Ovrlp_Sq,nBas,C_tmp,nBas,0.0d0,C_Virt,nBas)
call mma_deallocate(C_tmp)
call mma_deallocate(Ovrlp_Sq)

#ifdef _DEBUGPRINT_
call RecPrt('C_Virt(New)',' ',C_Virt,nBas,nVirt)
#endif

return

end subroutine Virt_Space
