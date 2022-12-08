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
! Copyright (C) 1987, Bjorn O. Roos                                    *
!               1992, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1987  B. O. ROOS                           *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND, SWEDEN                 *
!--------------------------------------------*

subroutine TR2Sq(CMO,X1,X2,X3,URPQ,RUPQ,TUPQ,lBuf)
! SECOND ORDER TWO-ELECTRON TRANSFORMATION ROUTINE
!
! THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
! (ISP,ISQ,ISR,ISS) WITH ISP >= ISQ AND ISR >= ISS.
! P,Q,R,S are SO indices.
! A,B are MO indices, counting only non-frozen and non-deleted.
! T,U are occupied MO indices, only non-frozen and non-deleted.
! INTEGRALS (AB/TU) ARE ALWAYS GENERATED
! EXCHANGE INTEGRALS (AT/BU) ARE GENERATED AS FOLLOWS:
! (AT/BU) IF ISP >= ISR
! (AT/UB) IF ISP > ISS AND ISP /= ISQ
! (TA/BU) IF ISQ > ISR AND ISP /= ISQ
! (TA/UB) IF ISQ >= ISS AND ISP /= ISQ
!
! ********** IBM-3090 RELEASE 87 09 14 **********
! Replace MXMA with DGEMM. P-AA Malmqvist 1992-05-06.

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: lBuf
#include "rasdim.fh"
#include "caspt2.fh"
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(_OUT_) :: X1(*), X2(*), X3(*)
real(kind=wp), intent(inout) :: URPQ(*), RUPQ(*), TUPQ(*)
integer(kind=iwp) :: IAD1, IAD1S, IAD2, IAD2S, IAD3, IAD3S, iOpt, IOUT1, IOUT2, IOUT3, IPQ, IPQMX1, IPQMX2, IPQMX3, IPQST, IR, &
                     iRc, IRSST, IRU, ISPQRS, IST, ITU, IX2, KKTU, LAR, LPQ, LR, NA, NAT, NORU, NOTU, NOUR, NP, NQ, NQM, NR, &
                     NSYMP, NT, NTM, NTMAX, NU, NUM, NUMAX
#include "intgrl.fh"
#include "trafo.fh"

NSYMP = (NSYM**2+NSYM)/2
NORU = NBR*NOCS
NOUR = NBS*NOCR
NOTU = NOCR*NOCS
if (ISR == ISS) NOTU = (NOCR**2+NOCR)/2

! CHECK FOR IN CORE OR OUT OF CORE TRANSFORMATION

! 1. SORT OF PARTIALLY TRANSFORMED INTEGRALS (RU/PQ) ON UNIT LUHLF1
IPQMX1 = NBPQ
if (NBPQ*NORU > LRUPQ) then
  IPQMX1 = LRUPQ/NORU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (RU/PQ)',IPQMX1
  IAD1S = 0
  call dDAFILE(LUHLF1,0,RUPQ,IPQMX1,IAD1S)
end if
IAD1 = 0
IOUT1 = 0
! 2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (UR/PQ) ON UNIT LUHLF2
IPQMX2 = NBPQ
if (NBPQ*NOUR > LURPQ) then
  IPQMX2 = LURPQ/NOUR
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (UR/PQ)',IPQMX2
  IAD2S = 0
  call dDAFILE(LUHLF2,0,URPQ,IPQMX2,IAD2S)
end if
IAD2 = 0
IOUT2 = 0
! 3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (TU/PQ) ON UNIT LUHLF3
IPQMX3 = NBPQ
if (NBPQ*NOTU > LTUPQ) then
  IPQMX3 = LTUPQ/NOTU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (TU/PQ)',IPQMX3
  IAD3S = 0
  call dDAFILE(LUHLF3,0,TUPQ,IPQMX3,IAD3S)
end if
IAD3 = 0
IOUT3 = 0

! START LOOP OVER SORTED AO-INTEGRALS: NPQ PQ-PAIRS IN EACH BUFFER

IPQ = 0
LPQ = 0
NPQ = 0
iOpt = 1
iRc = 0
IRSST = 1-NBRS
do NP=1,NBP
  NQM = NBQ
  if (ISP == ISQ) NQM = NP
  do NQ=1,NQM
    IPQ = IPQ+1
    IOUT1 = IOUT1+1
    IOUT2 = IOUT2+1
    IOUT3 = IOUT3+1

    ! READ IN A BLOCK OF INTEGRALS FOR NPQ PQ-VALUES

    if (LPQ == NPQ) then
      call RdOrd(iRc,iOpt,isP,isQ,isR,isS,X1,lBuf,nPQ)
      if (IRC > 1) then
        write(u6,*) ' ERROR RETURN CODE IRC=',IRC
        write(u6,*) ' FROM RDORD, CALLED FROM TRA2.'
        call Abend()
      end if
      iOpt = 2
      LPQ = 0
      IRSST = 1-NBRS
    end if
    LPQ = LPQ+1
    IRSST = IRSST+NBRS

    ! START TRANFORMATION OF THIS PQ PAIR

    if (ISR == ISS) then
      call SQUARE(X1(IRSST),X2,1,NBS,NBS)
    else
      call DCOPY_(NBRS,X1(IRSST),1,X2,1)
    end if

    ! INTEGRALS (PQ/UR)

    if (((ISP /= ISQ) .and. (ISP > ISS)) .and. (NOCR /= 0)) then
      call DGEMM_('N','N',NBS,NOCR,NBR,One,X2,NBS,CMO(LMOR2),NBR,Zero,X3,NBS)

      ! SORT THESE INTEGRALS AS (UR/PQ)

      if (IOUT2 > IPQMX2) then
        IOUT2 = 1
        call dDAFILE(LUHLF2,1,URPQ,IPQMX2*NOUR,IAD2)
      end if
      call DCOPY_(NOUR,X3,1,URPQ(IOUT2),IPQMX2)
    end if

    ! INTEGRALS (PQ/RU)

    if (NOCS /= 0) then
      call DGEMM_('T','N',NBR,NOCS,NBS,One,X2,NBS,CMO(LMOS2),NBS,Zero,X3,NBR)

      ! SORT THESE INTEGRALS AS (RU/PQ)

      if (ISP >= ISR) then
        if (IOUT1 > IPQMX1) then
          IOUT1 = 1
          call dDAFILE(LUHLF1,1,RUPQ,IPQMX1*NORU,IAD1)
        end if
        call DCOPY_(NORU,X3,1,RUPQ(IOUT1),IPQMX1)
      end if
    end if

    ! INTEGRALS (PQ/TU)

    if (NOCR*NOCS /= 0) then
      if (ISR == ISS) then
        call DGEMM_Tri('T','N',NOCR,NOCR,NBR,One,X3,NBR,CMO(LMOR2),NBR,Zero,X2,NOCR)
      else
        call DGEMM_('T','N',NOCS,NOCR,NBR,One,X3,NBR,CMO(LMOR2),NBR,Zero,X2,NOCS)
      end if

      ! SORT INTEGRALS (PQ/TU) INTO TUPQ (SORT AFTER PQ INSTEAD OF TU)

      if (IOUT3 > IPQMX3) then
        IOUT3 = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF3,1,TUPQ(1+IPQMX3*(I-1)),IPQMX3,IAD3)
        !vv end do
        call dDAFILE(LUHLF3,1,TUPQ,IPQMX3*NOTU,IAD3)

      end if
      call DCOPY_(NOTU,X2,1,TUPQ(IOUT3),IPQMX3)
    end if

    ! WE NOW HAVE THREE SETS OF PARTIALLY TRANSFORMED INTEGRALS
    ! IN TUPQ: (TU/PQ)   TRIANGULAR FOR ISR == ISS
    ! IN RUPQ: (RU/PQ)   IF ISP >= ISR
    ! IN URPQ: (UR/PQ)   IF (ISP > ISQ) .AND. (ISP > ISS)

  end do
end do

! EMPTY LAST BUFFERS

if (IPQMX1 < NBPQ) then
  !vv do I=1,NORU
  !vv   call dDAFILE(LUHLF1,1,RUPQ(1+IPQMX1*(I-1)),IPQMX1,IAD1)
  !vv end do
  call dDAFILE(LUHLF1,1,RUPQ,IPQMX1*NORU,IAD1)
end if
if (IPQMX2 < NBPQ) then
  !vv do I=1,NOUR
  !vv   call dDAFILE(LUHLF2,1,URPQ(1+IPQMX2*(I-1)),IPQMX2,IAD2)
  !vv end do
  call dDAFILE(LUHLF2,1,URPQ,IPQMX2*NOUR,IAD2)
end if
if (IPQMX3 < NBPQ) then
  !vv do I=1,NOTU
  !vv   call dDAFILE(LUHLF3,1,TUPQ(1+IPQMX3*(I-1)),IPQMX3,IAD3)
  !vv end do
  call dDAFILE(LUHLF3,1,TUPQ,IPQMX3*NOTU,IAD3)

end if

! FIRST PARTIAL TRANSFORMATION FINISHED
! SORTED INTEGRALS ARE ON UNITS LUHLF1 (RUPQ), LUHLF2 (URPQ),
! AND LUHLF3 (TUPQ), CONTROLLED BY THE ADRESSES IAD1,IAD2, AND IAD3,
! OR IN CORE (RUPQ, URPQ, AND TUPQ)

! SECOND HALF TRANSFORMATION FOR INTEGRALS (PQ/TU)
! FIRST SAVE THE START ADDRESS ON LUINTM FOR THIS BLOCK OF INTEGRALS
! NOTE THAT THE SYMMETRY LABEL ISPQRS ASSUMES THAT SYMMETRY LOOPS
! IN THE ORDER T,U,A,B FOR ALL INTEGRAL TYPES.

if (NOCR*NOCS /= 0) then
  ISPQRS = ((ISR**2-ISR)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISQ
  IAD2M(1,ISPQRS) = IAD13
  ITU = 0
  do NT=1,NOCR
    NUM = NT
    if (ISS /= ISR) NUM = NOCS
    do NU=1,NUM
      IPQST = 1+NBPQ*ITU
      ITU = ITU+1

      ! READ ONE BUFFER OF INTEGRALS BACK INTO CORE

      if (IPQMX3 < NBPQ) then
        !*sta0830
        call RBuf_tra2(LUHLF3,TUPQ,NBPQ,IPQMX3,NOTU,ITU,IPQST,IAD3S)
        !*end0830
      end if

      if (ISP == ISQ) then
        call SQUARE(TUPQ(IPQST),X2,1,NBQ,NBQ)
        call DGEMM_('N','N',NBQ,NOP,NBP,One,X2,NBQ,CMO(LMOP),NBP,Zero,X1,NBP)
        call DGEMM_Tri('T','N',NOP,NOP,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,ZERO,X2,NOP)
        IX2 = (NOP+NOP**2)/2
      else
        call DGEMM_('N','N',NBQ,NOP,NBP,One,TUPQ(IPQST),NBQ,CMO(LMOP),NBP,Zero,X1,NBQ)
        call DGEMM_('T','N',NOQ,NOP,NBQ,One,CMO(LMOQ),NBQ,X1,NBQ,Zero,X2,NOQ)
        IX2 = NOP*NOQ
      end if

      ! WRITE INTEGRALS (AB/TU) ON OUTPUT UNIT LUINTM
      ! INTEGRALS FOR SYMMETRY BLOCK (ISP,ISQ,ISR,ISS) ARE STORED
      ! ONE BLOCK FOR EACH TU STARTING AT ADDRESS IAD2M(1,ISPQRS).
      ! TRIANGULAR IN AB AND TU IF ISP == ISQ ( AND ISR == ISS)

      call GADSum(X2,IX2)
      call dDAFILE(LUINTM,1,X2,IX2,IAD13)

      ! EXTRACT INTEGRALS WITH ALL INDICES ACTIVE INTO TUVX
      ! Not used with caspt2. If wanted, TUVX must be declared.

      !if ((ISP >= ISR) .and. (NASH(ISP)*NASH(ISQ) /= 0)) call GTUVX(X2,TUVX,NT,NU,ITP,ITQ,ITR,ITS,ISP,ISQ)

    end do
  end do
end if

! SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/RU)-> (AT/BU)
! IF ISP == ISR THEN T >= U BUT ALWAYS ALL A AND B

NOTU = NOCQ*NOCS
if (ISQ == ISS) NOTU = (NOCQ**2+NOCQ)/2
if ((ISP >= ISR) .and. (NOTU /= 0)) then
  LAR = LTUPQ/NOTU
  LR = LAR/NOP
  if (LR > NBR) LR = NBR
  LAR = NOP*LR
  IAD3S = 0
  call dDAFILE(LUHLF3,0,TUPQ,LAR,IAD3S)
  IAD3 = 0
  IR = 0
  do NR=1,NBR
    IR = IR+1
    do NU=1,NOCS
      IRU = NBR*(NU-1)+NR
      IPQST = 1+NBPQ*(IRU-1)

      ! READ ONE BUFFER OF INTEGRALS BACK INTO CORE

      if (IPQMX1 < NBPQ) then
        !*sta0830
        call RBuf_tra2(LUHLF1,RUPQ,NBPQ,IPQMX1,NORU,IRU,IPQST,IAD1S)
        !*end0830
      end if
      if (ISP == ISQ) then
        call SQUARE(RUPQ(IPQST),X2,1,NBQ,NBQ)
      else
        call DCOPY_(NBPQ,RUPQ(IPQST),1,X2,1)
      end if
      if (ISQ == ISS) then
        call DGEMM_('T','N',NBP,NOCQ-NU+1,NBQ,One,X2,NBQ,CMO(LMOQ2+NBQ*(NU-1)),NBQ,Zero,X1,NBP)
        call DGEMM_('T','N',NOCQ-NU+1,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOCQ-NU+1)
      else
        call DGEMM_('T','N',NBP,NOCQ,NBQ,One,X2,NBQ,CMO(LMOQ2),NBQ,Zero,X1,NBP)
        call DGEMM_('T','N',NOCQ,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOCQ)
      end if

      ! INTEGRALS (AT/RU) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
      ! VALUES OF A,T (T >= U IF ISQ == ISS)
      !
      ! SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
      ! OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
      ! INTEGRALS (TU/PQ)

      if (IR > LR) then
        IR = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
        !vv end do
        call dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)
      end if
      NAT = 0
      do NA=1,NOP
        NTM = 1
        if (ISQ == ISS) NTM = NU
        do NT=NTM,NOCQ
          ITU = NOCS*(NT-1)+NU-1
          if (ISQ < ISS) ITU = NOCQ*(NU-1)+NT-1
          if (ISQ == ISS) ITU = (NT**2-NT)/2+NU-1
          NAT = NAT+1
          TUPQ(LAR*ITU+NOP*(IR-1)+NA) = X2(NAT)
        end do
      end do
    end do
  end do

  ! EMPTY LAST BUFFER IF LR < NBR

  if (LR < NBR) then
    !vv do I=1,NOTU
    !vv   call dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
    !vv end do
    call dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)
  end if

  ! NOW TRANSFORM INDEX R TO B

  ! FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK

  if (ISQ >= ISS) then
    ISPQRS = ((ISQ**2-ISQ)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISR
    IAD2M(2,ISPQRS) = IAD13
    NTMAX = NOCQ
    NUMAX = NOCS
  else
    ISPQRS = ((ISS**2-ISS)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISR
    IAD2M(3,ISPQRS) = IAD13
    NTMAX = NOCS
    NUMAX = NOCQ
  end if
  KKTU = 0
  IST = 1-NOP*NBR
  do NT=1,NTMAX
    NUM = NUMAX
    if (ISQ == ISS) NUM = NT
    do NU=1,NUM
      KKTU = KKTU+1
      IST = IST+NOP*NBR
      if (LR < NBR) then
        !*sta0830
        call RBuf_tra2(LUHLF3,TUPQ,NBR*NOP,LAR,NOTU,KKTU,IST,IAD3S)
        !*end0830
      end if
      call DGEMM_('T','T',NOR,NOP,NBR,One,CMO(LMOR),NBR,TUPQ(IST),NOP,Zero,X2,NOR)

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOP*NOR)
      call dDAFILE(LUINTM,1,X2,NOP*NOR,IAD13)
    end do
  end do
end if

! SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/RU)-> (TA/BU)

NOTU = NOCP*NOCS
if (((ISP /= ISQ) .and. (ISQ > ISR)) .and. (NOTU /= 0)) then
  LAR = LTUPQ/NOTU
  LR = LAR/NOQ
  if (LR > NBR) LR = NBR
  LAR = NOQ*LR
  IAD3S = 0
  call dDAFILE(LUHLF3,0,TUPQ,LAR,IAD3S)
  IAD3 = 0
  IRU = 0
  IR = 0
  do NR=1,NBR
    IR = IR+1
    do NU=1,NOCS
      IRU = NBR*(NU-1)+NR
      IPQST = 1+NBPQ*(IRU-1)

      ! READ ONE BUFFER OF INTEGRALS BACK INTO CORE

      if (IPQMX1 < NBPQ) then
        !*sta0830
        call RBuf_tra2(LUHLF1,RUPQ,NBPQ,IPQMX1,NORU,IRU,IPQST,IAD1S)
        !*end0830
      end if
      call DGEMM_('N','N',NBQ,NOCP,NBP,One,RUPQ(IPQST),NBQ,CMO(LMOP2),NBP,Zero,X1,NBQ)
      call DGEMM_('T','N',NOCP,NOQ,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,Zero,X2,NOCP)

      ! INTEGRALS (TA/RU) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
      ! VALUES OF T,A . NOTE THAT T AND U HERE ALWAYS ARE OF DIFFERENT
      ! SYMMETRIES
      !
      ! SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
      ! OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
      ! INTEGRALS (TU/PQ)

      if (IR > LR) then
        IR = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
        !vv end do
        call dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)

      end if
      NAT = 0
      do NA=1,NOQ
        do NT=1,NOCP
          ITU = NOCS*(NT-1)+NU-1
          NAT = NAT+1
          TUPQ(LAR*ITU+NOQ*(IR-1)+NA) = X2(NAT)
        end do
      end do
    end do
  end do

  ! EMPTY LAST BUFFER IF LR < NBR

  if (LR < NBR) then
    !vv do I=1,NOTU
    !vv   call dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
    !vv end do
    call dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)
  end if

  ! NOW TRANSFORM INDEX R TO B

  ! FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK

  ISPQRS = ((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISQ**2-ISQ)/2+ISR
  IAD2M(2,ISPQRS) = IAD13
  KKTU = 0
  IST = 1-NOQ*NBR
  do NT=1,NOCP
    do NU=1,NOCS
      KKTU = KKTU+1
      IST = IST+NOQ*NBR
      if (LR < NBR) then
        !*sta0830
        call RBuf_tra2(LUHLF3,TUPQ,NBR*NOQ,LAR,NOTU,KKTU,IST,IAD3S)
        !*end0830
      end if
      call DGEMM_('T','T',NOR,NOQ,NBR,One,CMO(LMOR),NBR,TUPQ(IST),NOQ,Zero,X2,NOR)

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOR*NOQ)
      call dDAFILE(LUINTM,1,X2,NOR*NOQ,IAD13)
    end do
  end do
end if

! SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/UR)-> (AT/UB)

NOTU = NOCQ*NOCR
if (((ISP /= ISQ) .and. (ISP > ISS)) .and. (NOTU /= 0)) then
  LAR = (LRUPQ+LTUPQ)/NOTU
  LR = LAR/NOP
  if (LR > NBS) LR = NBS
  LAR = NOP*LR
  IAD3S = 0
  call dDAFILE(LUHLF3,0,RUPQ,LAR,IAD3S)
  IAD3 = 0
  IRU = 0
  IR = 0
  do NR=1,NBS
    IR = IR+1
    do NU=1,NOCR
      IRU = NBS*(NU-1)+NR
      IPQST = 1+NBPQ*(IRU-1)

      ! READ ONE BUFFER OF INTEGRALS BACK INTO CORE

      if (IPQMX2 < NBPQ) then
        !*sta0830
        call RBuf_tra2(LUHLF2,URPQ,NBPQ,IPQMX2,NOUR,IRU,IPQST,IAD2S)
        !*end0830
      end if
      call DGEMM_('T','N',NBP,NOCQ,NBQ,One,URPQ(IPQST),NBQ,CMO(LMOQ2),NBQ,Zero,X1,NBP)
      call DGEMM_('T','N',NOCQ,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOCQ)

      ! INTEGRALS (AT/UR) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
      ! VALUES OF A,T. NOTE THAT T AND U HAVE DIFFERENT SYMMETRIES.
      !
      ! SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
      ! OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
      ! INTEGRALS (TU/PQ)

      if (IR > LR) then
        IR = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
        !vv end do
        call dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
      end if
      NAT = 0
      do NA=1,NOP
        do NT=1,NOCQ
          ITU = NOCR*(NT-1)+NU-1
          if (ISQ < ISR) ITU = NOCQ*(NU-1)+NT-1
          NAT = NAT+1
          RUPQ(LAR*ITU+NOP*(IR-1)+NA) = X2(NAT)
        end do
      end do
    end do
  end do

  ! EMPTY LAST BUFFER IF LR < NBS

  if (LR < NBS) then
    !vv do I=1,NOTU
    !vv   call dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
    !vv end dO
    call dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
  end if

  ! NOW TRANSFORM INDEX R TO B

  ! FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK

  if (ISQ >= ISR) then
    ISPQRS = ((ISQ**2-ISQ)/2+ISR-1)*NSYMP+(ISP**2-ISP)/2+ISS
    IAD2M(2,ISPQRS) = IAD13
    NTMAX = NOCQ
    NUMAX = NOCR
  else
    ISPQRS = ((ISR**2-ISR)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISS
    IAD2M(3,ISPQRS) = IAD13
    NTMAX = NOCR
    NUMAX = NOCQ
  end if
  KKTU = 0
  IST = 1-NBS*NOP
  do NT=1,NTMAX
    do NU=1,NUMAX
      KKTU = KKTU+1
      IST = IST+NBS*NOP
      if (LR < NBS) then
        !*sta0830
        call RBuf_tra2(LUHLF3,RUPQ,NBS*NOP,LAR,NOTU,KKTU,IST,IAD3S)
        !*end0830
      end if
      call DGEMM_('T','T',NOS,NOP,NBS,One,CMO(LMOS),NBS,RUPQ(IST),NOP,Zero,X2,NOS)

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOS*NOP)
      call dDAFILE(LUINTM,1,X2,NOS*NOP,IAD13)
    end do
  end do
end if

! SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/UR)-> (TA/UB)

NOTU = NOCP*NOCR
if (((ISP /= ISQ) .and. (ISQ >= ISS)) .and. (NOTU /= 0)) then
  if (ISP == ISR) NOTU = (NOCP**2+NOCP)/2
  LAR = (LRUPQ+LTUPQ)/NOTU
  LR = LAR/NOQ
  if (LR > NBS) LR = NBS
  LAR = NOQ*LR
  IAD3S = 0
  call dDAFILE(LUHLF3,0,RUPQ,LAR,IAD3S)
  IAD3 = 0
  IRU = 0
  IR = 0
  do NR=1,NBS
    IR = IR+1
    do NU=1,NOCR
      IRU = NBS*(NU-1)+NR
      IPQST = 1+NBPQ*(IRU-1)

      ! READ ONE BUFFER OF INTEGRALS BACK INTO CORE

      if (IPQMX2 < NBPQ) then
        !*sta0830
        call RBuf_tra2(LUHLF2,URPQ,NBPQ,IPQMX2,NOUR,IRU,IPQST,IAD2S)
        !*end0830
      end if
      if (ISP == ISR) then
        call DGEMM_('N','N',NBQ,NOCP-NU+1,NBP,One,URPQ(IPQST),NBQ,CMO(LMOP2+NBP*(NU-1)),NBP,Zero,X1,NBQ)
        call DGEMM_('T','N',NOCP-NU+1,NOQ,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,Zero,X2,NOCP-NU+1)
      else
        call DGEMM_('N','N',NBQ,NOCP,NBP,One,URPQ(IPQST),NBQ,CMO(LMOP2),NBP,Zero,X1,NBQ)
        call DGEMM_('T','N',NOCP,NOQ,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,Zero,X2,NOCP)
      end if

      ! INTEGRALS (TA/UR) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
      ! VALUES OF A,T
      !
      ! SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
      ! OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
      ! INTEGRALS (TU/PQ)

      if (IR > LR) then
        IR = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
        !vv end do
        call dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
      end if
      NAT = 0
      do NA=1,NOQ
        NTM = 1
        if (ISP == ISR) NTM = NU
        do NT=NTM,NOCP
          ITU = NOCR*(NT-1)+NU-1
          if (ISP < ISR) ITU = NOCP*(NU-1)+NT-1
          if (ISP == ISR) ITU = (NT**2-NT)/2+NU-1
          NAT = NAT+1
          RUPQ(LAR*ITU+NOQ*(IR-1)+NA) = X2(NAT)
        end do
      end do
    end do
  end do

  ! EMPTY LAST BUFFER IF LR < NBS

  if (LR < NBS) then
    !vv do I=1,NOTU
    !vv   call dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
    !vv end do
    call dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
  end if

  !NOW TRANSFORM INDEX R TO B

  ! FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK

  if (ISP >= ISR) then
    ISPQRS = ((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
    IAD2M(2,ISPQRS) = IAD13
    NTMAX = NOCP
    NUMAX = NOCR
  else
    ISPQRS = ((ISR**2-ISR)/2+ISP-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
    IAD2M(3,ISPQRS) = IAD13
    NTMAX = NOCR
    NUMAX = NOCP
  end if
  KKTU = 0
  IST = 1-NOQ*NBS
  do NT=1,NTMAX
    NUM = NUMAX
    if (ISP == ISR) NUM = NT
    do NU=1,NUM
      KKTU = KKTU+1
      IST = IST+NOQ*NBS
      if (LR < NBS) then
        !*sta0830
        call RBuf_tra2(LUHLF3,RUPQ,NBS*NOQ,LAR,NOTU,KKTU,IST,IAD3S)
        !*end0830
      end if
      call DGEMM_('T','T',NOS,NOQ,NBS,One,CMO(LMOS),NBS,RUPQ(IST),NOQ,Zero,X2,NOS)

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOS*NOQ)
      call dDAFILE(LUINTM,1,X2,NOS*NOQ,IAD13)
    end do
  end do
end if

! END OF TRANSFORMATION FOR THIS SYMMETRY BLOCK
!
! IAD2M CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:
! IAD2M(1,ISPQRS)   COULOMB INTEGRALS (AB|TU)
! IAD2M(2,ISPQRS)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T > SYM U
! IAD2M(3,ISPQRS)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T < SYM U
! THE LAST ADRESS IS ZERO IF SYM T = SYM U
! TO SEE HOW THE INTEGRALS ARE USED LOOK IN RDINT2

return

end subroutine TR2Sq
