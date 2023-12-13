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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine DELTAB(NREF,IOCR,L0,L1,L2,L3,INTNUM,LV,IFCORE,ICOR,NONE_,JONE,K00,K11,K22,K33,L00,L11,L22,L33)

use guga_global, only: IB, ILIM, IV0, IWAY, J2, K0, K1, K2, K3, LN, MXVERT, NIORB, NSM
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NREF, IOCR(*), INTNUM, LV, IFCORE, ICOR(*), NONE_, JONE(*), K00(*), K11(*), K22(*), K33(*), &
                                 L00(*), L11(*), L22(*), L33(*)
integer(kind=iwp), intent(_OUT_) :: L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp) :: I, IBS, IDIF, IEL, IIJ, IJJ, INHOLE, IOC(55), IPART, IREF, IRR, ISP(55), ISTA, JHOLE, JJ1, JPART, K, KM, KM1, &
                     LNS, LSYM, NCORR, NSJ
logical(kind=iwp) :: first, last
integer(kind=iwp), allocatable :: K0M(:), K1M(:), K2M(:), K3M(:), L0M(:), L1M(:), L2M(:), L3M(:)

K0(:) = 0
K1(:) = 0
K2(:) = 0
K3(:) = 0
do I=1,4*MXVERT
  L0(I) = 0
  L1(I) = 0
  L2(I) = 0
  L3(I) = 0
end do
IBS = 0
IEL = 2
LSYM = 1
LNS = NIORB+LV+1
IRR = 0
do I=LNS,LN
  IRR = IRR+1
  if (IOCR(IRR) == 1) LSYM = Mul(LSYM,NSM(I))
end do
call mma_allocate(K0M,MXVERT,label='K0M')
call mma_allocate(K1M,MXVERT,label='K1M')
call mma_allocate(K2M,MXVERT,label='K2M')
call mma_allocate(K3M,MXVERT,label='K3M')
call mma_allocate(L0M,MXVERT,label='L0M')
call mma_allocate(L1M,MXVERT,label='L1M')
call mma_allocate(L2M,MXVERT,label='L2M')
call mma_allocate(L3M,MXVERT,label='L3M')
do IIJ=1,ILIM
  ISTA = (IIJ-1)*MXVERT
  K0M(1:IV0) = 0
  K1M(1:IV0) = 0
  K2M(1:IV0) = 0
  K3M(1:IV0) = 0
  L0M(1:IV0) = 0
  L1M(1:IV0) = 0
  L2M(1:IV0) = 0
  L3M(1:IV0) = 0
  IJJ = IV0+1-IIJ
  KM = 1
  J2(KM) = IJJ
  first = .true.
  do
    if (first) then
      KM = KM+1
      IWAY(KM) = 0
      first = .false.
    end if
    KM1 = KM-1
    if ((L00(J2(KM1)) /= 0) .and. (IWAY(KM) < 1)) then
      J2(KM) = L00(J2(KM1))
      IWAY(KM) = 1
      IOC(KM1) = 0
      ISP(KM1) = 0
    else if ((L11(J2(KM1)) /= 0) .and. (IWAY(KM) < 2)) then
      J2(KM) = L11(J2(KM1))
      IWAY(KM) = 2
      IOC(KM1) = 1
      ISP(KM1) = 1
    else if ((L22(J2(KM1)) /= 0) .and. (IWAY(KM) < 3)) then
      J2(KM) = L22(J2(KM1))
      IWAY(KM) = 3
      IOC(KM1) = 1
      ISP(KM1) = 2
    else if ((L33(J2(KM1)) /= 0) .and. (IWAY(KM) < 4)) then
      J2(KM) = L33(J2(KM1))
      IWAY(KM) = 4
      IOC(KM1) = 2
      ISP(KM1) = 3
    else
      KM = KM-1
      if (KM == 1) then
        K0(ISTA+1:ISTA+IV0) = K0M(1:IV0)
        K1(ISTA+1:ISTA+IV0) = K1M(1:IV0)
        K2(ISTA+1:ISTA+IV0) = K2M(1:IV0)
        K3(ISTA+1:ISTA+IV0) = K3M(1:IV0)
        L0(ISTA+1:ISTA+IV0) = L0M(1:IV0)
        L1(ISTA+1:ISTA+IV0) = L1M(1:IV0)
        L2(ISTA+1:ISTA+IV0) = L2M(1:IV0)
        L3(ISTA+1:ISTA+IV0) = L3M(1:IV0)
        exit
      end if
      cycle
    end if
    if (KM1 == NIORB+LV) IBS = IB(J2(KM))
    if (KM /= LN+1) then
      first = .true.
      cycle
    end if
    NSJ = 1
    INHOLE = 0
    do I=1,LN
      if (IOC(I) == 1) NSJ = Mul(NSJ,NSM(I))
      if ((I <= NIORB+LV) .and. (I > LV)) INHOLE = INHOLE+2-IOC(I)
    end do
    ! STRIKE OUT INTERNAL CONFIGURATIONS
    IPART = 0
    if (IIJ > 1) IPART = IPART+1
    if (IIJ > 2) IPART = IPART+1
    JJ1 = 0
    last = .true.
    do IREF=1,NREF
      JHOLE = 0
      JPART = IPART
      do I=1,LN
        if (I <= LV) then
          IDIF = IOC(I)
        else if (I <= NIORB+LV) then
          IDIF = IOC(I)-2
        else
          JJ1 = JJ1+1
          if (IOC(I) == IOCR(JJ1)) cycle
          IDIF = IOC(I)-IOCR(JJ1)
        end if
        if (IDIF <= 0) then
          JHOLE = JHOLE-IDIF
        else
          JPART = JPART+IDIF
        end if
      end do
      if (JPART /= JHOLE) then
        write(u6,*) 'DeltaB: JPART /= JHOLE'
        write(u6,*) 'JPART,JHOLE=',JPART,JHOLE
        write(u6,*) 'iREF=',iREF
        call Abend()
      end if
      if (JPART <= IEL) then
        last = .false.
        exit
      end if
    end do
    if (last) cycle
    if ((IPART == 0) .and. (NSJ /= LSYM)) cycle
    if ((IPART == 2) .and. (INTNUM /= 0)) then
      ! INTERACTING SPACE
      if ((INHOLE == 2) .and. (IBS /= 0)) cycle
    end if
    ! NO CORE-CORE CORRELATION
    if (IFCORE /= 0) then
      NCORR = 0
      do I=1,LN
        if (ICOR(I) == 0) cycle
        NCORR = NCORR+2-IOC(I)
      end do
      if (NCORR > 1) cycle
    end if
    ! SINGLY OCCUPIED ORBITALS
    do I=1,NONE_
      if (IOC(JONE(I)) /= 1) cycle
    end do
    do K=1,LN
      if (ISP(K)-1 < 0) then
        K0M(J2(K+1)) = K00(J2(K+1))
        L0M(J2(K)) = L00(J2(K))
      else if (ISP(K)-1 == 0) then
        K1M(J2(K+1)) = K11(J2(K+1))
        L1M(J2(K)) = L11(J2(K))
      else if (ISP(K) == 3) then
        K3M(J2(K+1)) = K33(J2(K+1))
        L3M(J2(K)) = L33(J2(K))
      else
        K2M(J2(K+1)) = K22(J2(K+1))
        L2M(J2(K)) = L22(J2(K))
      end if
    end do
  end do
end do
call mma_deallocate(K0M)
call mma_deallocate(K1M)
call mma_deallocate(K2M)
call mma_deallocate(K3M)
call mma_deallocate(L0M)
call mma_deallocate(L1M)
call mma_deallocate(L2M)
call mma_deallocate(L3M)

return

end subroutine DELTAB
