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
! Copyright (C) 1989, Per Ake Malmqvist                                *
!***********************************************************************
!****************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST
!  SUBROUTINE TRACR     IBM-3090 RELEASE 89 01 30
!  TRANSFORM ONE SPECIFIC SYMMETRY BLOCK OF TWO-ELECTRON
!  INTEGRALS. THIS ROUTINE IS CALLED FROM TRAINT.
!****************************************************************

subroutine TRACR(LBUF,CMO1,CMO2,NGAM2,TUVX,X1,X2,X3,VXPQ)

use Index_Functions, only: iTri
use TRNSFRM, only: IAPR, ISP, ISQ, ISR, ISS, LMOP1, LMOQ1, LMOR1, LMOS1, NAP, NAQ, NAR, NAS, NAVX, NBP, NBPQ, NBQ, NBR, NBRS, NBS, &
                   NVXPQ, NX1MX, NX2MX, NX3MX
use rassi_data, only: NASHT, NCMO
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: LBUF, NGAM2
real(kind=wp) :: CMO1(NCMO), CMO2(NCMO), TUVX(NGAM2), X1(NX1MX), X2(NX2MX), X3(NX3MX), VXPQ(NVXPQ)
integer(kind=iwp) :: II, IOPT, IPQ, IPQST, IRC, IRSST, IT, ITF, ITU, IU, IUF, IUM, IV, IVF, IVX, IX, IXF, LPQ, NP, NPQ, NQ, NQM

! START LOOP OVER ORDERED AO-INTEGRALS: NPQ PQ-PAIRS IN EACH BUFFER.
! FOR EACH PQ PAIR, THERE IS A MATRIX CONTAINING THE (PQ,RS)
! INTEGRALS. INDEX S RUNS FASTEST, SO BY FORTRAN RULES, IF REGARDED
! AS A MATRIX JPQ(R,S), THEN THE MATRIX IS TRANSPOSED. IF ISR=ISS,
! THE MATRIX IS STORED IN ROW-MAJOR UNDERTRIANGULAR FORMAT.
! IN THE COMMENTS, NOTE THAT T,U, ETC DENOTES ACTUAL ORBITAL
! INDICES, WHILE IT,IU ETC ARE COUNTERS WITHIN THE SYMMETRY BLOCKS.
IRC = 0
IOPT = 1
IPQ = 0
LPQ = 0
NPQ = 0
IRSST = 1-NBRS
do NP=1,NBP
  NQM = NBQ
  if (ISP == ISQ) NQM = NP
  do NQ=1,NQM
    IPQ = IPQ+1
    ! IF NECESSARY, READ IN A FRESH INTEGRAL BUFFER OF NPQ MATRICES:
    if (LPQ == NPQ) then
      call RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,X1,LBUF,NPQ)
      IOPT = 2
      LPQ = 0
      IRSST = 1-NBRS
    end if
    LPQ = LPQ+1
    IRSST = IRSST+NBRS
    ! START TRANSFORMATION OF THIS PQ=PAIR
    if (ISR == ISS) then
      call SQUARE(X1(IRSST),X2,1,NBR,NBR)
    else
      X2(1:NBR*NBS) = X1(IRSST:IRSST+NBR*NBS-1)
    end if
    ! X2 CONTAINS THE MATRIX JPQ TRANSPOSED, X2(IS,IR)=JPQ(R,S)=(PQ/RS).
    call DGEMM_('T','N',NBR,NAS,NBS,One,X2,NBS,CMO2(LMOS1),NBS,Zero,X3,NBR)
    call DGEMM_('T','N',NAS,NAR,NBR,One,X3,NBR,CMO1(LMOR1),NBR,Zero,X2,NAS)
    ! X2 IS TRANSFORMED JPQ MATRIX, TRANSPOSED: X2(IX,IV)=(PQ/VX).
    ! SORT THE MATRIX X2 INTO VXPQ (SORT AFTER PQ INSTEAD OF VX).
    call DCOPY_(NAVX,X2,1,VXPQ(IPQ),NBPQ)
  end do
end do
! FIRST HALF TRANSFORMATION IS NOW DONE.VXPQ CONTAINS HALF TRANS-
! FORMED INTEGRALS FOR THIS SYMMETRY BLOCK: VXPQ(IPQ,IX,IV)=(PQ/VX).
! NOTE: V IS FIRST ORBITAL SET, X SECOND
! NOW TRANSFORM INDICES PQ TO TU FOR ALL (TU) >= (VX)
! FORTRAN ORDER IMPLIES THEN THAT U>=X AND, IF U=X, THEN T>=V.
IPQST = 1-NBPQ
do IV=1,NAR
  IVF = IV+IAPR(ISR)
  do IX=1,NAS
    IXF = IX+IAPR(ISS)
    IVX = IVF+NASHT*(IXF-1)
    IPQST = IPQST+NBPQ
    if (ISP == ISQ) then
      call SQUARE(VXPQ(IPQST),X3,1,NBP,NBP)
    else
      X3(1:NBPQ) = VXPQ(IPQST:IPQST+NBPQ-1)
    end if
    ! X3 IS HALF-TRANSFORMED JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
    ! WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
    IUM = 0
    if (ISQ == ISS) IUM = IX-1
    call DGEMM_('T','N',NBP,NAQ-IUM,NBQ,One,X3,NBQ,CMO2(LMOQ1+IUM*NBQ),NBQ,Zero,X1,NBP)
    ! X1(IU,IP) = (PU/VX) FOR THE GIVEN VX, ALL P, AND U>=X, WHERE P
    ! IS BASIS FUNCTIONS IN SYMMETRY ISP, U AND X ARE ACTIVE ORBITALS
    ! OF STATE 2 IN SYMMETRIES ISQ AND ISS, RESP., AND V IS ACTIVE ORB
    ! OF STATE 1 AND HAS SYMMETRY ISR.
    call DGEMM_('T','N',NAQ-IUM,NAP,NBP,One,X1,NBP,CMO1(LMOP1),NBP,Zero,X2,NAQ-IUM)
    ! X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
    ! X2(IU-IUM,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAQ,  IT=1,..,NAP.
    II = 0
    do IT=1,NAP
      ITF = IT+IAPR(ISP)
      do IU=IUM+1,NAQ
        IUF = IU+IAPR(ISQ)
        II = II+1
        ITU = ITF+NASHT*(IUF-1)
        TUVX(iTri(ITU,IVX)) = X2(II)
      end do
    end do
    if (ISP == ISQ) cycle
    ! X3 STILL CONTAINS THE JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
    ! PERM SYMMETRY OF P AND Q MEANS WE CAN LET THEM SWITCH IDENTITY.
    ! FOR REMAINDER OF THE LOOP, IF ISP /= ISQ, THEN REGARD X3 AS
    ! CONTAINING X3(IQ,IP)=(QP/VX) AND TRANSFORM AS BEFORE.
    ! WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
    IUM = 0
    if (ISP == ISS) IUM = IX-1
    call DGEMM_('N','N',NBQ,NAP-IUM,NBP,One,X3,NBQ,CMO2(LMOP1+IUM*NBP),NBP,Zero,X1,NBQ)
    ! X1(IQ,IU) = (QU/VX) FOR THE GIVEN VX, ALL Q, AND U>=X, WHERE Q
    ! IS BASIS FUNCTIONS IN SYMMETRY ISQ, U AND X ARE ACTIVE ORBITALS
    ! OF STATE 2 IN SYMMETRIES ISP AND ISS, RESP., AND V IS ACTIVE ORB
    ! OF STATE 1 AND HAS SYMMETRY ISR.
    call DGEMM_('T','N',NAP-IUM,NAQ,NBQ,One,X1,NBQ,CMO1(LMOQ1),NBQ,Zero,X2,NAP-IUM)
    ! X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
    ! X2(IU-IUM,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAP,  IT=1,..,NAQ.
    II = 0
    do IT=1,NAQ
      ITF = IT+IAPR(ISQ)
      do IU=IUM+1,NAP
        IUF = IU+IAPR(ISP)
        II = II+1
        ITU = ITF+NASHT*(IUF-1)
        TUVX(iTri(ITU,IVX)) = X2(II)
      end do
    end do
  end do
end do
if (ISR /= ISS) then
  ! NOW REPEAT IT ALL OVER AGAIN. THIS TIME, USE PERMUTATION
  ! SYMMETRY TO SWITCH IDENTITY OF BASIS INDICES R AND S.
  IRC = 0
  IOPT = 1
  IPQ = 0
  LPQ = 0
  NPQ = 0
  IRSST = 1-NBRS
  do NP=1,NBP
    NQM = NBQ
    if (ISP == ISQ) NQM = NP
    do NQ=1,NQM
      IPQ = IPQ+1
      ! IF NECESSARY, READ IN A FRESH INTEGRAL BUFFER OF NPQ MATRICES:
      if (LPQ == NPQ) then
        call RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,X1,LBUF,NPQ)
        IOPT = 2
        LPQ = 0
        IRSST = 1-NBRS
      end if
      LPQ = LPQ+1
      IRSST = IRSST+NBRS
      X2(1:NBR*NBS) = X1(IRSST:IRSST+NBR*NBS-1)
      ! X2 CONTAINS THE MATRIX JPQ , X2(IS,IR)=JPQ(S,R)=(PQ/SR).
      ! NOTE THAT SYMMETRY OF R IS ISS AND SYMMETRY OF S IS ISR NOW.
      call DGEMM_('N','N',NBS,NAR,NBR,One,X2,NBS,CMO2(LMOR1),NBR,Zero,X3,NBS)
      call DGEMM_('T','N',NAR,NAS,NBS,One,X3,NBS,CMO1(LMOS1),NBS,Zero,X2,NAR)
      ! X2 IS TRANSFORMED JPQ MATRIX, TRANSPOSED: X2(IX,IV)=(PQ/VX).
      ! SORT THE MATRIX X2 INTO VXPQ (SORT AFTER PQ INSTEAD OF VX).
      call DCOPY_(NAVX,X2,1,VXPQ(IPQ),NBPQ)
    end do
  end do
  ! AS BEFORE, EXCEPT THAT SYMMETRY OF V IS ISS, SYMMETRY OF X IS ISR.
  ! NOW TRANSFORM INDICES PQ TO TU FOR ALL (TU) >= (VX)
  ! FORTRAN ORDER IMPLIES THEN THAT U>=X AND, IF U=X, THEN T>=V.
  IPQST = 1-NBPQ
  do IV=1,NAS
    IVF = IV+IAPR(ISS)
    do IX=1,NAR
      IXF = IX+IAPR(ISR)
      IVX = IVF+NASHT*(IXF-1)
      IPQST = IPQST+NBPQ
      if (ISP == ISQ) then
        call SQUARE(VXPQ(IPQST),X3,1,NBP,NBP)
      else
        X3(1:NBPQ) = VXPQ(IPQST:IPQST+NBPQ-1)
      end if
      ! X3 IS HALF-TRANSFORMED JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
      ! WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
      IUM = 0
      if (ISQ == ISR) IUM = IX-1
      call DGEMM_('T','N',NBP,NAQ-IUM,NBQ,One,X3,NBQ,CMO2(LMOQ1+IUM*NBQ),NBQ,Zero,X1,NBP)
      ! X1(IU,IP) = (PU/VX) FOR THE GIVEN VX, ALL P, AND U>=X, WHERE P
      ! IS BASIS FUNCTIONS IN SYMMETRY ISP, U AND X ARE ACTIVE ORBITALS
      ! OF STATE 2 IN SYMMETRIES ISQ AND ISR, RESP., AND V IS ACTIVE ORB
      ! OF STATE 1 AND HAS SYMMETRY ISS.
      call DGEMM_('T','N',NAQ-IUM,NAP,NBP,One,X1,NBP,CMO1(LMOP1),NBP,Zero,X2,NAQ-IUM)
      ! X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
      ! X2(IU,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAQ,  IT=1,..,NAP.
      II = 0
      do IT=1,NAP
        ITF = IT+IAPR(ISP)
        do IU=IUM+1,NAQ
          IUF = IU+IAPR(ISQ)
          II = II+1
          ITU = ITF+NASHT*(IUF-1)
          TUVX(iTri(ITU,IVX)) = X2(II)
        end do
      end do
      if (ISP == ISQ) cycle
      ! X3 STILL CONTAINS THE JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
      ! PERM SYMMETRY OF P AND Q MEANS WE CAN LET THEM SWITCH IDENTITY.
      ! FOR REMAINDER OF THE LOOP, IF ISP /= ISQ, THEN REGARD X3 AS
      ! CONTAINING X3(IQ,IP)=(QP/VX) AND TRANSFORM AS BEFORE.
      ! WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
      IUM = 0
      if (ISP == ISR) IUM = IX-1
      call DGEMM_('N','N',NBQ,NAP-IUM,NBP,One,X3,NBQ,CMO2(LMOP1+IUM*NBP),NBP,Zero,X1,NBQ)
      ! X1(IQ,IU) = (QU/VX) FOR THE GIVEN VX, ALL Q, AND U>=X, WHERE Q
      ! IS BASIS FUNCTIONS IN SYMMETRY ISQ, U AND X ARE ACTIVE ORBITALS
      ! OF STATE 2 IN SYMMETRIES ISP AND ISR, RESP., AND V IS ACTIVE ORB
      ! OF STATE 1 AND HAS SYMMETRY ISS.
      call DGEMM_('T','N',NAP-IUM,NAQ,NBQ,One,X1,NBQ,CMO1(LMOQ1),NBQ,Zero,X2,NAP-IUM)
      ! X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
      ! X2(IU,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAP,  IT=1,..,NAQ.
      II = 0
      do IT=1,NAQ
        ITF = IT+IAPR(ISQ)
        do IU=IUM+1,NAP
          IUF = IU+IAPR(ISP)
          II = II+1
          ITU = ITF+NASHT*(IUF-1)
          TUVX(iTri(ITU,IVX)) = X2(II)
        end do
      end do
    end do
  end do
end if

end subroutine TRACR
