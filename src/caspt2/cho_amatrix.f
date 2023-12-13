************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2023, Ignacio Fdez. Galvan                             *
************************************************************************

      SUBROUTINE Cho_Amatrix(XMAT,CMO,DDTR,NATR)
! Calculation of the "exchange" matrix for the G1,G2,G3 Fock operators
! from Cholesky vectors

      USE Index_Functions, ONLY: nTri_Elem
      USE Data_Structures, ONLY: Allocate_DT, Deallocate_DT, DSBA_Type
      USE CHOVEC_IO, ONLY: NVLOC_CHOBATCH
      USE stdalloc, ONLY: mma_allocate, mma_deallocate
      USE Constants, ONLY: Zero, One, Half

      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER :: NATR
      REAL*8 :: XMAT(NOSQT), CMO(NCMO), DDTR(NATR)
      INTEGER :: I, IB1, IB2, IBGRP, ISYM, J, JSYM, MXBGRP, MXCHOBUF,
     &           MXINT, MXPIQK, NADDBUFF, NBUF, NCHOBUF, NINTS, NLB,
     &           NLK, NV
      INTEGER, ALLOCATABLE :: BGRP(:,:,:), ICA(:), ICI(:), ICV(:),
     &                        IXMAT(:), NBGRP(:), NVEC(:,:)
      REAL*8, ALLOCATABLE :: BRABUF(:), KETBUF(:)
      REAL*8, ALLOCATABLE, TARGET :: INTBUF(:)
      TYPE(DSBA_Type) :: HDSQ
      INTEGER, PARAMETER :: Inac=1, Acti=2, Virt=3

      ! Transform Cholesky vectors, this will have to be redone after
      ! the modified Fock matrix is diagonalized
      CALL TRACHO3(CMO)

      ! Square (Dd) to simplify multiplication
      CALL Allocate_DT(HDSQ,NASH,NASH,NSYM,Label='HDSQ')
      I = 1
      DO ISYM=1,NSYM
        CALL Square(DDTR(I),HDSQ%Sb(ISYM)%A1,1,NASH(ISYM),NASH(ISYM))
        I = I+nTri_Elem(NASH(ISYM))
      END DO
      ! To include a 1/2 factor in the final matrix, just scale (Dd)
      HDSQ%A0(:) = Half*HDSQ%A0(:)

      ! Get offsets and sizes
      CALL mma_allocate(IXMAT,NSYM,Label='IXMAT')
      CALL mma_allocate(NBGRP,NSYM,Label='NBGRP')
      MXINT = 0
      MXBGRP = 0
      DO ISYM=1,NSYM
        ! Offsets of symmetry blocks in XMAT
        IF (ISYM == 1) THEN
          IXMAT(ISYM) = 1
        ELSE
          IXMAT(ISYM) = IXMAT(ISYM-1)+NORB(ISYM-1)**2
        END IF
        ! Max size for integral matrix
        DO JSYM=1,NSYM
          I = MUL(ISYM,JSYM)
          NINTS = MAX(NISH(I),NASH(I),NSSH(I))
          MXINT = MAX(MXINT,(NINTS*NASH(JSYM))**2)
        END DO
        MXBGRP = MAX(MXBGRP,NBTCH(ISYM))
      END DO
      CALL mma_allocate(INTBUF,MXINT,Label='INTBUF')
      CALL mma_allocate(BGRP,2,MXBGRP,NSYM,Label='BGRP')
      MXBGRP = 0
      MXCHOBUF = 0
      DO ISYM=1,NSYM
        ! Batch groups per symmetry
        NBGRP(ISYM) = NBTCH(ISYM)
        DO I=1,NBGRP(ISYM)
          BGRP(:,I,ISYM) = NBTCHES(ISYM)+I
        END DO
        ! Max size for Cholesky vectors
        CALL MEMORY_ESTIMATE(ISYM,BGRP(:,:,ISYM),NBGRP(ISYM),NCHOBUF,
     &                       MXPIQK,NADDBUFF)
        MXBGRP = MAX(MXBGRP,NBGRP(ISYM))
        MXCHOBUF = MAX(MXCHOBUF,NCHOBUF)
      END DO
      ! Number of Cholesky vectors per group and symmetry
      CALL mma_allocate(NVEC,MXBGRP,NSYM,Label='NVEC')
      NVEC(:,:) = 0
      DO ISYM=1,NSYM
        DO IBGRP=1,NBGRP(ISYM)
          DO I=BGRP(1,IBGRP,ISYM),BGRP(2,IBGRP,ISYM)
            NVEC(IBGRP,ISYM) = NVEC(IBGRP,ISYM)+NVLOC_CHOBATCH(I)
          END DO
        END DO
      END DO

      CALL mma_allocate(BRABUF,MXCHOBUF,Label='BRABUF')
      CALL mma_allocate(KETBUF,MXCHOBUF,Label='KETBUF')

      ! Now fetch Cholesky vectors and accumulate the A_pq matrix
      CALL mma_allocate(ICI,NSYM,Label='ICI')
      CALL mma_allocate(ICA,NSYM,Label='ICA')
      CALL mma_allocate(ICV,NSYM,Label='ICV')
      DO ISYM=1,NSYM
        ! Index of symmetry blocks in Cholesky vectors,
        ! according to the active orbital symmetry
        J = 0
        ICI(MUL(ISYM,1)) = J
        ICA(1) = 0
        ICV(1) = 0
        DO JSYM=1,NSYM-1
          I = MUL(ISYM,JSYM)
          J = J+NISH(JSYM)*NASH(I)
          ICI(MUL(ISYM,JSYM+1)) = J
          ICA(JSYM+1) = ICA(JSYM)+NASH(JSYM)*NASH(I)
          ICV(JSYM+1) = ICV(JSYM)+NASH(JSYM)*NSSH(I)
        END DO
        DO IBGRP=1,NBGRP(ISYM)
          IB1 = BGRP(1,IBGRP,ISYM)
          IB2 = BGRP(2,IBGRP,ISYM)
          NV = NVEC(IBGRP,ISYM)
          IF (NV == 0) EXIT
          ! Inactive-Inactive
          CALL Get_Cholesky_Vectors(Inac,Acti,ISYM,BRABUF,NBUF,IB1,IB2)
          CALL Accum(Inac,Inac,BRABUF,BRABUF,ICI,ICI)
          ! Inactive-Active
          CALL Get_Cholesky_Vectors(Acti,Acti,ISYM,KETBUF,NBUF,IB1,IB2)
          CALL Accum(Inac,Acti,BRABUF,KETBUF,ICI,ICA)
          ! Active-Active
          CALL Accum(Acti,Acti,KETBUF,KETBUF,ICA,ICA)
          ! Inactive-Virtual
          CALL Get_Cholesky_Vectors(Virt,Acti,ISYM,KETBUF,NBUF,IB1,IB2)
          CALL Accum(Inac,Virt,BRABUF,KETBUF,ICI,ICV)
          ! Virtual-Virtual
          CALL Accum(Virt,Virt,KETBUF,KETBUF,ICV,ICV)
          ! Active-Virtual
          ! We could have saved these Cholesky vectors,
          ! but there's joy in repetition
          CALL Get_Cholesky_Vectors(Acti,Acti,ISYM,BRABUF,NBUF,IB1,IB2)
          CALL Accum(Acti,Virt,BRABUF,KETBUF,ICA,ICV)
        END DO
      END DO

      CALL GADSum(XMAT,NOSQT)

      CALL Deallocate_DT(HDSQ)
      CALL mma_deallocate(IXMAT)
      CALL mma_deallocate(NBGRP)
      CALL mma_deallocate(BGRP)
      CALL mma_deallocate(NVEC)
      CALL mma_deallocate(BRABUF)
      CALL mma_deallocate(KETBUF)
      CALL mma_deallocate(INTBUF)
      CALL mma_deallocate(ICI)
      CALL mma_deallocate(ICA)
      CALL mma_deallocate(ICV)

      RETURN

      CONTAINS

      SUBROUTINE Accum(bBlock,kBlock,bBuf,kBuf,IB,IK)

      INTEGER :: bBlock, kBlock, IB(NSYM), IK(NSYM)
      REAL*8 :: bBuf(*), kBuf(*)
      INTEGER :: B1, BS, BSWCH, bOff(NSYM), I, II, IJ, IJT, J, JA, JJ,
     &           K1, KS, KSWCH, kOff(NSYM), NA, NB(NSYM), NK(NSYM),
     &           PQSYM, TUSYM
      LOGICAL :: diag
      REAL*8, POINTER, CONTIGUOUS :: INT2(:,:)
      REAL*8, EXTERNAL :: dDot_

      ! NB,NK = number of orbitals in bra/ket (not including NA factor)
      ! bOff,kOff = offset or starting orbital in bra/ket
      ! QB,QK = index function for bra/ket
      ! BSWCH,KSWCH = aux switch for generalizing integral access
      !               (1 if inactive, which come before active)
      BSWCH = 0
      SELECT CASE (bBlock)
        CASE (Inac)
          NB(:) = NISH(1:NSYM)
          bOff(:) = 0
          BSWCH = 1
        CASE (Acti)
          NB(:) = NASH(1:NSYM)
          bOff(:) = NISH(1:NSYM)
        CASE (Virt)
          NB(:) = NSSH(1:NSYM)
          bOff(:) = NISH(1:NSYM)+NASH(1:NSYM)
        CASE DEFAULT ! Nothing compares 2 U
          ! (just to keep compilers happy)
          NB(:) = 0
          bOff(:) = 0
          CALL Abend()
      END SELECT
      KSWCH = 0
      SELECT CASE (kBlock)
        CASE (Inac)
          NK(:) = NISH(1:NSYM)
          kOff(:) = 0
          KSWCH = 1
        CASE (Acti)
          NK(:) = NASH(1:NSYM)
          kOff(:) = NISH(1:NSYM)
        CASE (Virt)
          NK(:) = NSSH(1:NSYM)
          kOff(:) = NISH(1:NSYM)+NASH(1:NSYM)
        CASE DEFAULT
          NK(:) = 0
          kOff(:) = 0
          CALL Abend()
      END SELECT
      ! Is this a diagonal block?
      diag = bBlock == kBlock

      ! We want (pt|qu) integrals of symmetry ISYM, with:
      !   t,u of symmetry TUSYM
      !   p,q of symmetry PQSYM
      DO TUSYM=1,NSYM
        PQSYM = MUL(ISYM,TUSYM)
        NA = NASH(TUSYM)
        ! Reconstruct the (bBlock,Active|kBlock,Active) integrals
        NLB = NB(PQSYM)*NA
        NLK = NK(PQSYM)*NA
        CALL dgemm_('N','T',NLB,NLK,NV,One,bBUF(IB(TUSYM)*NV+1),NLB,
     &              kBUF(IK(TUSYM)*NV+1),NLK,Zero,INTBUF,NLB)
        IF (NLB*NLK == 0) CYCLE
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
        DO J=1,NK(PQSYM)
          DO I=1,NB(PQSYM)
            ! Index of this element in XMAT
            IJ = IXMAT(PQSYM)+
     &           (kOff(PQSYM)+J-1)*NORB(PQSYM)+bOff(PQSYM)+I-1
            ! Index of (p1|**) and (**|q1)
            II = (I-1)*B1+1
            JJ = (J-1)*K1+1
            DO JA=1,NA
              XMAT(IJ) = XMAT(IJ)+dDot_(NA,INT2(II:,JJ),BS,
     &                   HDSQ%SB(TUSYM)%A2(:,JA),1)
              JJ = JJ+KS
            END DO
            ! Symmetric element
            IF (diag .AND. (I == J)) EXIT
            IJT = IXMAT(PQSYM)+
     &            (bOff(PQSYM)+I-1)*NORB(PQSYM)+kOff(PQSYM)+J-1
            XMAT(IJT) = XMAT(IJ)
          END DO
        END DO
        NULLIFY(INT2)
      END DO

      END SUBROUTINE Accum

      END SUBROUTINE Cho_Amatrix
