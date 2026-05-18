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
! Copyright (C) 2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Cho_Amatrix(XMAT,nXMAT,CMO,NCMO,DDTR,NATR)
! Calculation of the "exchange" matrix for the G1,G2,G3 Fock operators
! from Cholesky vectors

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use CHOVEC_IO, only: NVLOC_CHOBATCH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use caspt2_module, only: nSym, nIsh, nAsh, nSsh, nOSqT, nOrb, nBtch, nBtches
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: nXMAT, nCMO, NATR
real(kind=wp), intent(inout) :: XMAT(nXMAT)
real(kind=wp), intent(in) :: CMO(nCMO), DDTR(NATR)
integer(kind=iwp) :: I, IB1, IB2, IBGRP, ISYM, J, JSYM, MXBGRP, MXCHOBUF, MXINT, MXPIQK, NADDBUFF, NBUF, NCHOBUF, NINTS, NLB, NLK, &
                     NV
integer(kind=iwp), allocatable :: BGRP(:,:,:), ICA(:), ICI(:), ICV(:), IXMAT(:), NBGRP(:), NVEC(:,:)
real(kind=wp), allocatable :: BRABUF(:), KETBUF(:)
real(kind=wp), allocatable, target :: INTBUF(:)
type(DSBA_Type) :: HDSQ
integer(kind=iwp), parameter :: Inac = 1, Acti = 2, Virt = 3

! Transform Cholesky vectors, this will have to be redone after
! the modified Fock matrix is diagonalized
call TRACHO3(CMO,NCMO)

! Square (Dd) to simplify multiplication
call Allocate_DT(HDSQ,NASH,NASH,NSYM,Label='HDSQ')
I = 1
do ISYM=1,NSYM
  call Square(DDTR(I),HDSQ%Sb(ISYM)%A1,1,NASH(ISYM),NASH(ISYM))
  I = I+nTri_Elem(NASH(ISYM))
end do
! To include a 1/2 factor in the final matrix, just scale (Dd)
HDSQ%A0(:) = Half*HDSQ%A0(:)

! Get offsets and sizes
call mma_allocate(IXMAT,NSYM,Label='IXMAT')
call mma_allocate(NBGRP,NSYM,Label='NBGRP')
MXINT = 0
MXBGRP = 0
do ISYM=1,NSYM
  ! Offsets of symmetry blocks in XMAT
  if (ISYM == 1) then
    IXMAT(ISYM) = 1
  else
    IXMAT(ISYM) = IXMAT(ISYM-1)+NORB(ISYM-1)**2
  end if
  ! Max size for integral matrix
  do JSYM=1,NSYM
    I = Mul(ISYM,JSYM)
    NINTS = max(NISH(I),NASH(I),NSSH(I))
    MXINT = max(MXINT,(NINTS*NASH(JSYM))**2)
  end do
  MXBGRP = max(MXBGRP,NBTCH(ISYM))
end do
call mma_allocate(INTBUF,MXINT,Label='INTBUF')
call mma_allocate(BGRP,2,MXBGRP,NSYM,Label='BGRP')
MXBGRP = 0
MXCHOBUF = 0
do ISYM=1,NSYM
  ! Batch groups per symmetry
  NBGRP(ISYM) = NBTCH(ISYM)
  do I=1,NBGRP(ISYM)
    BGRP(:,I,ISYM) = NBTCHES(ISYM)+I
  end do
  ! Max size for Cholesky vectors
  call MEMORY_ESTIMATE(ISYM,BGRP(:,:,ISYM),NBGRP(ISYM),NCHOBUF,MXPIQK,NADDBUFF)
  MXBGRP = max(MXBGRP,NBGRP(ISYM))
  MXCHOBUF = max(MXCHOBUF,NCHOBUF)
end do
! Number of Cholesky vectors per group and symmetry
call mma_allocate(NVEC,MXBGRP,NSYM,Label='NVEC')
NVEC(:,:) = 0
do ISYM=1,NSYM
  do IBGRP=1,NBGRP(ISYM)
    do I=BGRP(1,IBGRP,ISYM),BGRP(2,IBGRP,ISYM)
      NVEC(IBGRP,ISYM) = NVEC(IBGRP,ISYM)+NVLOC_CHOBATCH(I)
    end do
  end do
end do

call mma_allocate(BRABUF,MXCHOBUF,Label='BRABUF')
call mma_allocate(KETBUF,MXCHOBUF,Label='KETBUF')

! Now fetch Cholesky vectors and accumulate the A_pq matrix
call mma_allocate(ICI,NSYM,Label='ICI')
call mma_allocate(ICA,NSYM,Label='ICA')
call mma_allocate(ICV,NSYM,Label='ICV')
do ISYM=1,NSYM
  ! Index of symmetry blocks in Cholesky vectors,
  ! according to the active orbital symmetry
  J = 0
  ICI(Mul(ISYM,1)) = J
  ICA(1) = 0
  ICV(1) = 0
  do JSYM=1,NSYM-1
    I = Mul(ISYM,JSYM)
    J = J+NISH(JSYM)*NASH(I)
    ICI(Mul(ISYM,JSYM+1)) = J
    ICA(JSYM+1) = ICA(JSYM)+NASH(JSYM)*NASH(I)
    ICV(JSYM+1) = ICV(JSYM)+NASH(JSYM)*NSSH(I)
  end do
  do IBGRP=1,NBGRP(ISYM)
    IB1 = BGRP(1,IBGRP,ISYM)
    IB2 = BGRP(2,IBGRP,ISYM)
    NV = NVEC(IBGRP,ISYM)
    if (NV == 0) exit
    ! Inactive-Inactive
    call Get_Cholesky_Vectors(Inac,Acti,ISYM,BRABUF,size(BRABUF),NBUF,IB1,IB2)
    call Accum(Inac,Inac,BRABUF,size(BRABUF),BRABUF,size(BRABUF),ICI,ICI)
    ! Inactive-Active
    call Get_Cholesky_Vectors(Acti,Acti,ISYM,KETBUF,size(KETBUF),NBUF,IB1,IB2)
    call Accum(Inac,Acti,BRABUF,size(BRABUF),KETBUF,size(KETBUF),ICI,ICA)
    ! Active-Active
    call Accum(Acti,Acti,KETBUF,size(KETBUF),KETBUF,size(KETBUF),ICA,ICA)
    ! Inactive-Virtual
    call Get_Cholesky_Vectors(Virt,Acti,ISYM,KETBUF,size(KETBUF),NBUF,IB1,IB2)
    call Accum(Inac,Virt,BRABUF,size(BRABUF),KETBUF,size(KETBUF),ICI,ICV)
    ! Virtual-Virtual
    call Accum(Virt,Virt,KETBUF,size(KETBUF),KETBUF,size(KETBUF),ICV,ICV)
    ! Active-Virtual
    ! We could have saved these Cholesky vectors,
    ! but there's joy in repetition
    call Get_Cholesky_Vectors(Acti,Acti,ISYM,BRABUF,size(BRABUF),NBUF,IB1,IB2)
    call Accum(Acti,Virt,BRABUF,size(BRABUF),KETBUF,size(KETBUF),ICA,ICV)
  end do
end do

call GADGOp(XMAT,NOSQT,'+')

call Deallocate_DT(HDSQ)
call mma_deallocate(IXMAT)
call mma_deallocate(NBGRP)
call mma_deallocate(BGRP)
call mma_deallocate(NVEC)
call mma_deallocate(BRABUF)
call mma_deallocate(KETBUF)
call mma_deallocate(INTBUF)
call mma_deallocate(ICI)
call mma_deallocate(ICA)
call mma_deallocate(ICV)

contains

subroutine Accum(bBlock,kBlock,bBuf,nbBuf,kBuf,nkBuf,IB,IK)

  integer(kind=iwp), intent(in) :: nbBuf, nkBuf
  integer(kind=iwp) :: bBlock, kBlock, IB(NSYM), IK(NSYM)
  real(kind=wp) :: bBuf(nbBuf), kBuf(nkBuf)
  integer(kind=iwp) :: B1, BS, BSWCH, bOff(NSYM), I, II, IJ, IJT, J, JA, JJ, K1, KS, KSWCH, kOff(NSYM), NA, NB(NSYM), NK(NSYM), &
                       PQSYM, TUSYM
  logical(kind=iwp) :: diag
  real(kind=wp), pointer, contiguous :: INT2(:,:)
  real(kind=wp), external :: dDot_

  ! NB,NK = number of orbitals in bra/ket (not including NA factor)
  ! bOff,kOff = offset or starting orbital in bra/ket
  ! QB,QK = index function for bra/ket
  ! BSWCH,KSWCH = aux switch for generalizing integral access
  !               (1 if inactive, which come before active)
  BSWCH = 0
  select case (bBlock)
    case (Inac)
      NB(:) = NISH(1:NSYM)
      bOff(:) = 0
      BSWCH = 1
    case (Acti)
      NB(:) = NASH(1:NSYM)
      bOff(:) = NISH(1:NSYM)
    case (Virt)
      NB(:) = NSSH(1:NSYM)
      bOff(:) = NISH(1:NSYM)+NASH(1:NSYM)
    case default ! Nothing compares 2 U
      ! (just to keep compilers happy)
      NB(:) = 0
      bOff(:) = 0
      call Abend()
  end select
  KSWCH = 0
  select case (kBlock)
    case (Inac)
      NK(:) = NISH(1:NSYM)
      kOff(:) = 0
      KSWCH = 1
    case (Acti)
      NK(:) = NASH(1:NSYM)
      kOff(:) = NISH(1:NSYM)
    case (Virt)
      NK(:) = NSSH(1:NSYM)
      kOff(:) = NISH(1:NSYM)+NASH(1:NSYM)
    case default
      NK(:) = 0
      kOff(:) = 0
      call Abend()
  end select
  ! Is this a diagonal block?
  diag = bBlock == kBlock

  ! We want (pt|qu) integrals of symmetry ISYM, with:
  !   t,u of symmetry TUSYM
  !   p,q of symmetry PQSYM
  do TUSYM=1,NSYM
    PQSYM = Mul(ISYM,TUSYM)
    NA = NASH(TUSYM)
    ! Reconstruct the (bBlock,Active|kBlock,Active) integrals
    NLB = NB(PQSYM)*NA
    NLK = NK(PQSYM)*NA
    call dgemm_('N','T',NLB,NLK,NV,One,bBUF(IB(TUSYM)*NV+1),NLB,kBUF(IK(TUSYM)*NV+1),NLK,Zero,INTBUF,NLB)
    if (NLB*NLK == 0) cycle
    INT2(1:NLB,1:NLK) => INTBUF(1:NLB*NLK)

    ! Accumulate A_pq = sum_tu (pt|qu)*(Dd)_tu
    ! BS,KS = step sizes for traversing active orbitals in INT2
    !         (= 1 if inactive, = NB,NK otherwise)
    BS = (NB(PQSYM)-1)*(1-BSWCH)+1
    KS = (NK(PQSYM)-1)*(1-KSWCH)+1
    ! B1,K1 = for getting the index of the 1st active orbital
    !         (= NA if inactive, = 1 otherwise)
    B1 = (NA-1)*BSWCH+1
    K1 = (NA-1)*KSWCH+1
    do J=1,NK(PQSYM)
      do I=1,NB(PQSYM)
        ! Index of this element in XMAT
        IJ = IXMAT(PQSYM)+(kOff(PQSYM)+J-1)*NORB(PQSYM)+bOff(PQSYM)+I-1
        ! Index of (p1|**) and (**|q1)
        II = (I-1)*B1+1
        JJ = (J-1)*K1+1
        do JA=1,NA
          XMAT(IJ) = XMAT(IJ)+dDot_(NA,INT2(II:,JJ),BS,HDSQ%SB(TUSYM)%A2(:,JA),1)
          JJ = JJ+KS
        end do
        ! Symmetric element
        if (diag .and. (I == J)) exit
        IJT = IXMAT(PQSYM)+(bOff(PQSYM)+I-1)*NORB(PQSYM)+kOff(PQSYM)+J-1
        XMAT(IJT) = XMAT(IJ)
      end do
    end do
    nullify(INT2)
  end do

end subroutine Accum

end subroutine Cho_Amatrix
