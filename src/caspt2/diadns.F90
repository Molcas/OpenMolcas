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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DIADNS(ISYM,ICASE,VEC1,nVec1,VEC2,nVec2,DPT2,nDPT2,LIST,mList)
! Compute diagonal-block contribs to a trans density matrix.
! Vector blocks are in spectral resolution basis (ON).
! Each square matrix block of density matrix elements is
! computed and stored in full, even if VEC1=VEC2.
! Present implementation does not compute active-active
! contributions. This should be added in a separate routine,
! since it requires transformation to standard (Non-ON) basis.

use Symmetry_Info, only: Mul
use caspt2_global, only: do_grad
use EQSOLV, only: LLIST, NLIST
use Sigma_data, only: IFTEST, INCX1, INCX2, INCX3, INCY1, INCY2, LEN1, NLST1, VAL1
use caspt2_module, only: NAGEB, NAGTB, NASH, NASUP, NIGEJ, NIGTJ, NIMX, NINDEP, NISH, NISUP, NORB, NORB, NSMX, NSSH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISYM, ICASE, nVEC1, nVec2, nDPT2, mList, LIST(mList)
real(kind=wp), intent(in) :: VEC1(nVec1), VEC2(nVec2)
real(kind=wp), intent(inout) :: DPT2(nDPT2)
integer(kind=iwp) :: IA, IB, ICD, ICEM, ICEP, ICGM, ICGP, IDAB, IDII, IDIJ, IDTU, II, III, IJ, IJS, INCA, IOFCD(8,8), IOFDAB(8), &
                     IOFDIJ(8), IS, ISYMA, ISYMC, ISYMCA, ISYMI, ISYMK, ISYMKI, IV, IV1, IV11, IV2, IV22, IY, IY1, IY2, IYOFF, JS, &
                     LLST1, MU, NA, NAKI, NAS, NC, NCA, NCAY, NI, NICA, NIN, NIS, NK, NKI, NKIY, NO, NOA, NOI, NS, NVEC, NX
real(kind=wp) :: OVL, rSUM
real(kind=wp), allocatable :: X1(:), X2(:)
real(kind=wp), parameter :: SQR2 = sqrt(Two)
real(kind=wp), external :: DDOT_

NIN = NINDEP(ISYM,ICASE)
if (NIN == 0) return
NIS = NISUP(ISYM,ICASE)
if (NIS == 0) return
NAS = NASUP(ISYM,ICASE)
NVEC = NIN*NIS

IFTEST = 0
! Set up various offset arrays:
IDIJ = 0
do IS=1,NSYM
  NI = NISH(IS)
  NA = NASH(IS)
  NO = NORB(IS)
  IDTU = IDIJ+NO*NI+NI
  IDAB = IDTU+NO*NA+NA
  IOFDIJ(IS) = IDIJ
  IOFDAB(IS) = IDAB
  IDIJ = IDIJ+NO*NO
end do
do IS=1,NSYM
  ICD = 0
  ICEP = 0
  ICEM = 0
  ICGP = 0
  ICGM = 0
  do JS=1,NSYM
    IJS = Mul(IS,JS)
    IOFCD(IS,JS) = ICD
    ICD = ICD+NSSH(JS)*NISH(IJS)
    ICEP = ICEP+NSSH(JS)*NIGEJ(IJS)
    ICEM = ICEM+NSSH(JS)*NIGTJ(IJS)
    ICGP = ICGP+NISH(JS)*NAGEB(IJS)
    ICGM = ICGM+NISH(JS)*NAGTB(IJS)
  end do
end do

! Core contribution:
if (.not. do_grad) then
  OVL = DDOT_(NVEC,VEC1,1,VEC2,1)
  do IS=1,NSYM
    NI = NISH(IS)
    NO = NORB(IS)
    IDII = IOFDIJ(IS)+1
    do III=1,NI
      DPT2(IDII) = DPT2(IDII)+Two*OVL
      IDII = IDII+NO+1
    end do
  end do
end if

LLST1 = 0 ! dummy initialize
NLST1 = 0 ! dummy initialize

select case (ICASE)
  ! -----------------------------------------------
  case (1)
    ! Case A
    NI = NISH(ISYM)
    NO = NORB(ISYM)
    do II=1,NI
      IV2 = 1+NIN*(II-1)
      do IJ=1,NI
        IDIJ = IOFDIJ(ISYM)+II+NO*(IJ-1)
        IV1 = 1+NIN*(IJ-1)
        DPT2(IDIJ) = DPT2(IDIJ)-DDOT_(NIN,VEC1(IV1),1,VEC2(IV2),1)
      end do
    end do
    ! -----------------------------------------------
  case (2,3)
    ! Case BP
    ! Case BM
    ! Unfold VEC1 and VEC2 into X1(MU,K,I), X2(MU,K,I):
    NX = NIN*NIMX**2
    call mma_allocate(X1,NX,LABEL='X1')
    call mma_allocate(X2,NX,LABEL='X2')
    do ISYMK=1,NSYM
      ISYMI = Mul(ISYMK,ISYM)
      NK = NISH(ISYMK)
      NI = NISH(ISYMI)
      NKI = NK*NI
      if (NKI == 0) cycle
      X1(1:NIN*NKI) = Zero
      X2(1:NIN*NKI) = Zero
      if (ICASE == 2) then
        LLST1 = LLIST(ISYMK,ISYM,14)
        NLST1 = NLIST(ISYMK,ISYM,14)
        VAL1(1) = One
        VAL1(2) = SQR2
      else if (ICASE == 3) then
        LLST1 = LLIST(ISYMK,ISYM,15)
        NLST1 = NLIST(ISYMK,ISYM,15)
        VAL1(1) = One
        VAL1(2) = -One
      end if
      if (NLST1 == 0) cycle
      INCX1 = 1
      INCX2 = NIN
      INCX3 = NIN*NK
      INCY1 = 1
      INCY2 = NIN
      LEN1 = NIN
      call MLTUNF(LIST(LLST1),NLST1,X1,nX,VEC1,nVec1)
      call MLTUNF(LIST(LLST1),NLST1,X2,nX,VEC2,nVec2)
      ! D(I,J) := Add contraction -X2(MU,K,I)*X1(MU,K,J):
      IDIJ = 1+IOFDIJ(ISYMI)
      NO = NORB(ISYMI)
      call DGEMM_('T','N',NI,NI,NIN*NK,-One,X2,NIN*NK,X1,NIN*NK,One,DPT2(IDIJ),NO)
    end do
    call mma_deallocate(X1)
    call mma_deallocate(X2)
    ! -----------------------------------------------
  case (4)
    ! Case C
    NS = NSSH(ISYM)
    NO = NORB(ISYM)
    do IA=1,NS
      IV1 = 1+NIN*(IA-1)
      do IB=1,NS
        IDAB = IOFDAB(ISYM)+IA+NO*(IB-1)
        IV2 = 1+NIN*(IB-1)
        DPT2(IDAB) = DPT2(IDAB)+DDOT_(NIN,VEC1(IV1),1,VEC2(IV2),1)
      end do
    end do
    ! -----------------------------------------------
  case (5)
    ! Case D
    do ISYMA=1,NSYM
      NS = NSSH(ISYMA)
      NOA = NORB(ISYMA)
      ISYMI = Mul(ISYMA,ISYM)
      NI = NISH(ISYMI)
      NOI = NORB(ISYMI)
      IV = 1+NIN*IOFCD(ISYM,ISYMA)
      INCA = NIN*NI
      if (NI == 0) cycle
      do II=1,NI
        IV2 = IV+NIN*(II-1)
        do IJ=1,NI
          IDIJ = IOFDIJ(ISYMI)+II+NOI*(IJ-1)
          rSUM = DPT2(IDIJ)
          IV1 = IV+NIN*(IJ-1)
          do IA=1,NS
            IV11 = IV1+INCA*(IA-1)
            IV22 = IV2+INCA*(IA-1)
            rSUM = rSUM-DDOT_(NIN,VEC1(IV11),1,VEC2(IV22),1)
          end do
          DPT2(IDIJ) = rSUM
        end do
      end do
      do IA=1,NS
        IV1 = IV+INCA*(IA-1)
        do IB=1,NS
          IDAB = IOFDAB(ISYMA)+IA+NOA*(IB-1)
          IV2 = IV+INCA*(IB-1)
          DPT2(IDAB) = DPT2(IDAB)+DDOT_(INCA,VEC1(IV1),1,VEC2(IV2),1)
        end do
      end do
    end do
    ! -----------------------------------------------
  case (6,7)
    ! Case EP
    ! Case EM
    NX = NIN*NSMX*NIMX**2
    call mma_allocate(X1,NX,LABEL='X1')
    call mma_allocate(X2,NX,LABEL='X2')
    IYOFF = 0
    do ISYMA=1,NSYM
      ISYMKI = Mul(ISYMA,ISYM)
      NA = NSSH(ISYMA)
      NOA = NORB(ISYMA)
      if (ICASE == 6) NKIY = NIGEJ(ISYMKI)
      if (ICASE == 7) NKIY = NIGTJ(ISYMKI)
      IY = 1+IYOFF
      ! First, contributions to DIJ.
      ! Unfold VEC1 and VEC2 into X1(MU,A;K,I), X2(MU,A;K,I):
      do ISYMK=1,NSYM
        ISYMI = Mul(ISYMK,ISYMKI)
        NK = NISH(ISYMK)
        NI = NISH(ISYMI)
        NAKI = NA*NK*NI
        if (NAKI == 0) cycle
        X1(1:NIN*NAKI) = Zero
        X2(1:NIN*NAKI) = Zero
        if (ICASE == 6) then
          LLST1 = LLIST(ISYMK,ISYMKI,14)
          NLST1 = NLIST(ISYMK,ISYMKI,14)
          VAL1(1) = One
          VAL1(2) = SQR2
        else if (ICASE == 7) then
          LLST1 = LLIST(ISYMK,ISYMKI,15)
          NLST1 = NLIST(ISYMK,ISYMKI,15)
          VAL1(1) = One
          VAL1(2) = -One
        end if
        if (NLST1 == 0) cycle
        INCX1 = 1
        INCX2 = NIN*NA
        INCX3 = NIN*NA*NK
        INCY1 = 1
        INCY2 = NIN*NA
        LEN1 = NIN*NA
        call MLTUNF(LIST(LLST1),NLST1,X1,nX,VEC1(IY),nVEC1-IY+1)
        call MLTUNF(LIST(LLST1),NLST1,X2,nX,VEC2(IY),nVEC2-IY+1)
        ! D(I,J) := Add contraction -X2(MU,A,K,I)*X1(MU,A,K,J):
        IDIJ = 1+IOFDIJ(ISYMI)
        NOI = NORB(ISYMI)
        call DGEMM_('T','N',NI,NI,NIN*NA*NK,-One,X2,NIN*NA*NK,X1,NIN*NA*NK,One,DPT2(IDIJ),NOI)
      end do
      ! Second, contributions to DAB.
      if (NKIY == 0) cycle
      do IA=1,NA
        do IB=1,NA
          IDAB = IOFDAB(ISYMA)+IA+NOA*(IB-1)
          rSUM = DPT2(IDAB)
          do MU=1,NIN
            IY1 = IYOFF+MU+NIN*(IA-1)
            IY2 = IYOFF+MU+NIN*(IB-1)
            rSUM = rSUM+DDOT_(NKIY,VEC1(IY1),INCY2,VEC2(IY2),INCY2)
          end do
          DPT2(IDAB) = rSUM
        end do
      end do
      IYOFF = IYOFF+NIN*NA*NKIY
    end do
    call mma_deallocate(X1)
    call mma_deallocate(X2)
    ! -----------------------------------------------
  case (8,9)
    ! Case FP
    ! Case FM
    ! Unfold VEC1 and VEC2 into X1(MU,C,A), X2(MU,C,B):
    NX = NIN*NSMX**2
    call mma_allocate(X1,NX,LABEL='X1')
    call mma_allocate(X2,NX,LABEL='X2')
    do ISYMC=1,NSYM
      ISYMA = Mul(ISYMC,ISYM)
      NC = NSSH(ISYMC)
      NA = NSSH(ISYMA)
      NCA = NC*NA
      if (NCA == 0) cycle
      X1(1:NIN*NCA) = Zero
      X2(1:NIN*NCA) = Zero
      if (ICASE == 8) then
        LLST1 = LLIST(ISYMC,ISYM,16)
        NLST1 = NLIST(ISYMC,ISYM,16)
        VAL1(1) = One
        VAL1(2) = SQR2
      else if (ICASE == 9) then
        LLST1 = LLIST(ISYMC,ISYM,17)
        NLST1 = NLIST(ISYMC,ISYM,17)
        VAL1(1) = One
        VAL1(2) = -One
      end if
      if (NLST1 == 0) cycle
      INCX1 = 1
      INCX2 = NIN
      INCX3 = NIN*NC
      INCY1 = 1
      INCY2 = NIN
      LEN1 = NIN
      call MLTUNF(LIST(LLST1),NLST1,X1,nX,VEC1,nVec1)
      call MLTUNF(LIST(LLST1),NLST1,X2,nX,VEC2,nVec2)
      ! D(A,B) := Add contraction  X1(MU,C,A)*X2(MU,C,B):
      IDAB = 1+IOFDAB(ISYMA)
      NOA = NORB(ISYMA)
      call DGEMM_('T','N',NA,NA,NIN*NC,+One,X1,NIN*NC,X2,NIN*NC,One,DPT2(IDAB),NOA)
    end do
    call mma_deallocate(X1)
    call mma_deallocate(X2)
    ! -----------------------------------------------
  case (10,11)
    ! Case GP
    ! Case GM
    NX = NIN*NIMX*NSMX**2
    call mma_allocate(X1,NX,LABEL='X1')
    call mma_allocate(X2,NX,LABEL='X2')
    IYOFF = 0
    do ISYMI=1,NSYM
      ISYMCA = Mul(ISYMI,ISYM)
      NI = NISH(ISYMI)
      NOI = NORB(ISYMI)
      if (ICASE == 10) NCAY = NAGEB(ISYMCA)
      if (ICASE == 11) NCAY = NAGTB(ISYMCA)
      IY = 1+IYOFF
      ! First, contributions to DAB.
      ! Unfold VEC1 and VEC2 into X1(MU,I;C,A), X2(MU,I;C,A):
      do ISYMC=1,NSYM
        ISYMA = Mul(ISYMC,ISYMCA)
        NC = NSSH(ISYMC)
        NA = NSSH(ISYMA)
        NICA = NI*NC*NA
        if (NICA == 0) cycle
        X1(1:NIN*NICA) = Zero
        X2(1:NIN*NICA) = Zero
        if (ICASE == 10) then
          LLST1 = LLIST(ISYMC,ISYMCA,16)
          NLST1 = NLIST(ISYMC,ISYMCA,16)
          VAL1(1) = One
          VAL1(2) = SQR2
        else if (ICASE == 11) then
          LLST1 = LLIST(ISYMC,ISYMCA,17)
          NLST1 = NLIST(ISYMC,ISYMCA,17)
          VAL1(1) = One
          VAL1(2) = -One
        end if
        if (NLST1 == 0) cycle
        INCX1 = 1
        INCX2 = NIN*NI
        INCX3 = NIN*NI*NC
        INCY1 = 1
        INCY2 = NIN*NI
        LEN1 = NIN*NI
        call MLTUNF(LIST(LLST1),NLST1,X1,nX,VEC1(IY),nVec1-IY+1)
        call MLTUNF(LIST(LLST1),NLST1,X2,nX,VEC2(IY),nVec2-IY+1)
        ! D(A,B) := Add contraction +X1(MU,I,C,A)*X2(MU,I,C,B):
        IDAB = 1+IOFDAB(ISYMA)
        NOA = NORB(ISYMA)
        call DGEMM_('T','N',NA,NA,NIN*NI*NC,+One,X1,NIN*NI*NC,X2,NIN*NI*NC,One,DPT2(IDAB),NOA)
      end do
      ! Second, contributions to DIJ.
      if (NCAY == 0) cycle
      do II=1,NI
        do IJ=1,NI
          IDIJ = IOFDIJ(ISYMI)+II+NOI*(IJ-1)
          rSUM = DPT2(IDIJ)
          do MU=1,NIN
            IY1 = IYOFF+MU+NIN*(IJ-1)
            IY2 = IYOFF+MU+NIN*(II-1)
            rSUM = rSUM-DDOT_(NCAY,VEC1(IY1),INCY2,VEC2(IY2),INCY2)
          end do
          DPT2(IDIJ) = rSUM
        end do
      end do
      IYOFF = IYOFF+NIN*NI*NCAY
    end do
    call mma_deallocate(X1)
    call mma_deallocate(X2)
    ! -----------------------------------------------
  case (12,13)
    ! Case HP
    ! Case HM
    ! Unfold VEC1 and VEC2 into X1(MU,K,I), X2(MU,K,I):
    NX = NAS*NIMX**2
    call mma_allocate(X1,NX,LABEL='X1')
    call mma_allocate(X2,NX,LABEL='X2')
    do ISYMK=1,NSYM
      ISYMI = Mul(ISYMK,ISYM)
      NK = NISH(ISYMK)
      NI = NISH(ISYMI)
      NKI = NK*NI
      if (NKI == 0) cycle
      X1(1:NAS*NKI) = Zero
      X2(1:NAS*NKI) = Zero
      if (ICASE == 12) then
        LLST1 = LLIST(ISYMK,ISYM,14)
        NLST1 = NLIST(ISYMK,ISYM,14)
        VAL1(1) = One
        VAL1(2) = SQR2
      else if (ICASE == 13) then
        LLST1 = LLIST(ISYMK,ISYM,15)
        NLST1 = NLIST(ISYMK,ISYM,15)
        VAL1(1) = One
        VAL1(2) = -One
      end if
      if (NLST1 == 0) cycle
      INCX1 = 1
      INCX2 = NAS
      INCX3 = NAS*NK
      INCY1 = 1
      INCY2 = NAS
      LEN1 = NAS
      call MLTUNF(LIST(LLST1),NLST1,X1,nX,VEC1,nVec1)
      call MLTUNF(LIST(LLST1),NLST1,X2,nX,VEC2,nVec2)
      ! D(I,J) := Add contraction -X2(MU,K,I)*X1(MU,K,J):
      IDIJ = 1+IOFDIJ(ISYMI)
      NOI = NORB(ISYMI)
      call DGEMM_('T','N',NI,NI,NAS*NK,-One,X2,NAS*NK,X1,NAS*NK,One,DPT2(IDIJ),NOI)
    end do
    call mma_deallocate(X1)
    call mma_deallocate(X2)
    ! Unfold VEC1 and VEC2 into X1(A,C,IJ), X2(A,C,IJ):
    NX = NIS*NSMX**2
    call mma_allocate(X1,NX,LABEL='X1')
    call mma_allocate(X2,NX,LABEL='X2')
    do ISYMC=1,NSYM
      ISYMA = Mul(ISYMC,ISYM)
      NC = NSSH(ISYMC)
      NA = NSSH(ISYMA)
      NCA = NC*NA
      if (NCA == 0) cycle
      X1(1:NIS*NCA) = Zero
      X2(1:NIS*NCA) = Zero
      if (ICASE == 12) then
        LLST1 = LLIST(ISYMA,ISYM,16)
        NLST1 = NLIST(ISYMA,ISYM,16)
        VAL1(1) = One
        VAL1(2) = SQR2
      else if (ICASE == 13) then
        LLST1 = LLIST(ISYMA,ISYM,17)
        NLST1 = NLIST(ISYMA,ISYM,17)
        VAL1(1) = One
        VAL1(2) = -One
      end if
      if (NLST1 == 0) cycle
      INCX1 = NCA
      INCX2 = 1
      INCX3 = NA
      INCY1 = NAS
      INCY2 = 1
      LEN1 = NIS
      call MLTUNF(LIST(LLST1),NLST1,X1,nX,VEC1,nVec1)
      call MLTUNF(LIST(LLST1),NLST1,X2,nX,VEC2,nVec2)
      ! D(A,B) := Add contraction  X1(A,C,IJ)*X2(B,C,IJ):
      IDAB = 1+IOFDAB(ISYMA)
      NOA = NORB(ISYMA)
      call DGEMM_('N','T',NA,NA,NIS*NC,+One,X1,NA,X2,NA,One,DPT2(IDAB),NOA)
    end do
    call mma_deallocate(X1)
    call mma_deallocate(X2)
  case default
    call ABEND()
end select

contains

subroutine MLTUNF(LST,nLST,X,nX,Y,nY)
  ! Given a list with entries LST(4,ITEM), ITEM=1,NLST,
  ! with entries called L1,L2,L3,L4 for given ITEM, and
  ! an array of the form Y(p,q), compute the matrix
  !    X(p,L1,L2) := Add V*Y(p,L3), p=1..LEN1
  ! where V=VAL1(L4), looped over ITEM=1,NLST.
  ! Note: Arrays are addressed by strides given in common.

  integer(kind=iwp), intent(in) :: nLST, LST(4,nLST), nX, nY
  real(kind=wp), intent(inout) :: X(nX)
  real(kind=wp), intent(in) :: Y(nY)
  integer(kind=iwp) :: ILST, IX, IY, L1, L2, L3, L4
  real(kind=wp) :: V

  do ILST=1,NLST
    L1 = LST(1,ILST)
    L2 = LST(2,ILST)
    L3 = LST(3,ILST)
    L4 = LST(4,ILST)
    V = VAL1(L4)
    IX = 1+INCX2*(L1-1)+INCX3*(L2-1)
    IY = 1+INCY2*(L3-1)
    call DAXPY_(LEN1,V,Y(IY),INCY1,X(IX),INCX1)
  end do

end subroutine MLTUNF

end subroutine DIADNS
