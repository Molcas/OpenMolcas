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
! Copyright (C) 2000, Per Ake Malmqvist                                *
!***********************************************************************
!****************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST 2000-06-30
!  SUBROUTINE TRAINT
!  TRANSFORM TWO-ELECTRON INTEGRALS TO A MIXED ACTIVE ORBITAL
!  BASIS. THE RESULTING INTEGRALS ARE STORED IN ARRAY TUVX IN
!  THE SAME FORMAT AS THE 2-EL TRANSITION DENSITY MATRICES.
!****************************************************************

subroutine TRAINT(CMO1,CMO2,NGAM2,TUVX)

use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate
use TRNSFRM, only: NAP, NBP, ISP, LMOP1, NAQ, NBQ, ISQ, LMOQ1, NAR, NBR, LMOR1, ISR, NAS, NBS, LMOS1, ISS, NBPQ, NBRS, NAVX, &
                   NX1MX, NX2MX, NX3MX, NVXPQ, IAPR
use Symmetry_Info, only: nSym => nIrrep, MUL
use rassi_data, only: NCMO, NBMX, NASH, NBASF, NISH, NOSH

implicit none
integer NGAM2
real*8 CMO1(NCMO), CMO2(NCMO), TUVX(NGAM2)
integer KEEP(8), NBSX(8)
logical ISQARX
real*8, allocatable :: X1(:), X2(:), X3(:), VXPQ(:)
integer IRC, IA, IS, INTBUF, LMOP, NSP, LMOQ, NSQ, NSPQ, ISPQ, LMOR, NSRM, NSR, NSPQR, LMOS, NSSM, NSS, ISRS, NACT, NSYMX

! CLEAR THE ARRAY OF TRANSFORMED INTEGRALS.
TUVX(:) = Zero
! RETRIEVE STRUCTURE DATA FOR THE ORDERED INTEGRAL FILE.
! KEEP(IS) IS POSITIVE TO INDICATE THAT BASIS FUNCTIONS WITH
! SYMMETRY LABEL IS SHOULD BE SKIPPED. KEEP(2)=0 OR 1 TO SHOW
! IF SYMMETRY LABELS IJKL OF EACH SYMMETRY BLOCK ARE CANONICAL
! (ONLY CASES I>=J,K>=L,IJ>=KL PRESENT) OR IF LEFT AND RIGHT HAND
! PAIR ARE INDIVIDUALLY TREATED (I>=J,K>=L), RESPECTIVELY.
! THE IDATA ARRAY GIVES DISK ADDRESS, AND NUMBER OF INTEGRAL MATRICES
! PER BUFFER, FOR EACH SYMMETRY BLOCK. THERE ARE AT MOST 176 SUCH
! BLOCKS, AND THEIR ORDERING IS DETERMINED BY THE SYMMETRY
! LOOPS, WHICH MUST BE THE SAME AS IN THE ORDERING PROGRAM.
! RETRIEVE BASE DATA FROM UNIT LUORD:
! RETRIEVE BASE DATA FROM UNIT LUORD:
IRC = 0
call GETORD(IRC,ISQARX,NSYMX,NBSX,KEEP)
! SET UP IAPR(IS)=NR OF ACTIVES WITH PREVIOUS SYMMETRY LABEL.
IA = 0
do IS=1,NSYM
  IAPR(IS) = IA
  IA = IA+NASH(IS)
end do
! LOOP OVER QUADRUPLES OF SYMMETRIES (NSP,NSQ,NSR,NSS) NSR>=NSS
! IN THE SAME ORDER AS IN THE INTORD PROGRAM.
INTBUF = max(NBMX**2,256*256)
LMOP = 1
do NSP=1,NSYM
  if (NSP /= 1) LMOP = LMOP+NBASF(NSP-1)*NOSH(NSP-1)
  NAP = NASH(NSP)
  !KEEPP = KEEP(NSP)
  NBP = NBASF(NSP)
  ISP = NSP
  LMOP1 = LMOP+NISH(NSP)*NBP
  LMOQ = 1
  do NSQ=1,NSP
    if (NSQ /= 1) LMOQ = LMOQ+NBASF(NSQ-1)*NOSH(NSQ-1)
    NAQ = NASH(NSQ)
    !KEEPQ = KEEP(NSQ)
    NBQ = NBASF(NSQ)
    NSPQ = MUL(NSP,NSQ)
    ISQ = NSQ
    ISPQ = (ISP**2-ISP)/2+ISQ
    LMOQ1 = LMOQ+NISH(NSQ)*NBQ
    LMOR = 1
    NSRM = NSYM
    if (ISQARX) NSRM = NSP
    do NSR=1,NSRM
      if (NSR /= 1) LMOR = LMOR+NBASF(NSR-1)*NOSH(NSR-1)
      NAR = NASH(NSR)
      !KEEPR = KEEP(NSR)
      NBR = NBASF(NSR)
      LMOR1 = LMOR+NBR*NISH(NSR)
      NSPQR = MUL(NSPQ,NSR)
      ISR = NSR
      LMOS = 1
      NSSM = NSR
      do NSS=1,NSSM
        if (NSS /= 1) LMOS = LMOS+NBASF(NSS-1)*NOSH(NSS-1)
        if (NSPQR /= NSS) GO TO 101
        NAS = NASH(NSS)
        !KEEPS = KEEP(NSS)
        NBS = NBASF(NSS)
        LMOS1 = LMOS+NBS*NISH(NSS)
        ISS = NSS
        ISRS = (ISR**2-ISR)/2+ISS
        ! SHOULD THIS SYMMETRY BLOCK BE USED...?
        !KEEPT = KEEPP+KEEPQ+KEEPR+KEEPS
        NACT = NAP*NAQ*NAR*NAS
        ! ...WELL, IT IS PRESENT ON THE FILE
        if (ISPQ < ISRS) GO TO 101
        if (NACT == 0) GO TO 101
        ! ALLOCATE WORK AREAS FOR IN-CORE TRANSFORMATION ROUTINE TRACR:
        if (ISR == ISS) then
          NBPQ = (NBP+NBP**2)/2
          NBRS = (NBR+NBR**2)/2
        else
          NBPQ = NBP*NBQ
          NBRS = NBR*NBS
        end if
        NAVX = NAR*NAS
        NX1MX = max(INTBUF,NBP*NAQ,NBQ*NAP)
        NX2MX = max(NBR*NBS,NAP*NAQ)
        NX3MX = max(NBR*NAS,NAR*NBS,NBP*NBQ)
        NVXPQ = NAVX*NBPQ
        call mma_allocate(X1,NX1MX,Label='X1')
        call mma_allocate(X2,NX2MX,Label='X2')
        call mma_allocate(X3,NX3MX,Label='X3')
        call mma_allocate(VXPQ,NVXPQ,Label='VXPQ')
        call TRACR(INTBUF,CMO1,CMO2,NGAM2,TUVX,X1,X2,X3,VXPQ)
        call mma_deallocate(X1)
        call mma_deallocate(X2)
        call mma_deallocate(X3)
        call mma_deallocate(VXPQ)
101     continue
      end do
    end do
  end do
end do
call GADGOp(TUVX,NGAM2,'+')

end subroutine TRAINT
