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
!***********************************************************************
!--------------------------------------------*
! 1987  B. O. ROOS                           *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND, SWEDEN                 *
!--------------------------------------------*

subroutine tr2NsA2(CMO,X1,nX1,X2,nX2,pqrU,npqrU,pqTU,npqTU)
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
! This and tr2NsB routines transform non-squared AO integrals. The
! transformed MO integrals are stored as the same as Tr2Sq
! subroutine does.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "rasdim.fh"
#include "caspt2.fh"
integer(kind=iwp), intent(in) :: nX1, nX2, npqrU, npqTU
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(out) :: X1(nX1), X2(nX2)
real(kind=wp), intent(inout) :: pqrU(npqrU), pqTU(npqTU)
integer(kind=iwp) :: IAD2S, IAD3, IAD3S, icc, icxc1, icxc5, IPQMX2, IPQMX3, IPQST, IPQTU, IR, IRU, ISPQRS, IST, ITU, IX2, KKTU, &
                     LAR, LR, NA, NAT, NORU, NOTU, NR, NSYMP, NT, NTM, NTMAX, NU, Num, NUMAX
#include "trafo.fh"
#include "intgrl.fh"

NSYMP = NSYM*(NSYM+1)/2
NOTU = NOCR*NOCS
if (ISR == ISS) NOTU = (NOCR**2+NOCR)/2
NORU = NBR*NOCS
icc = NOP*NOQ*NOCR*NOCS
icxc1 = NOP*NOCQ*NOR*NOCS
icxc5 = NOCP*NOQ*NOR*NOCS

! Check for in core or out of core transformation

! 2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
IPQMX2 = NBPQ
if (NBPQ*NORU > LRUPQ) then
  IPQMX2 = LRUPQ/NORU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
  IAD2S = 0
  call dDAFILE(LUHLF2,0,PQRU,IPQMX2,IAD2S)
end if
!4x     3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
IPQMX3 = NBPQ
if (NBPQ*NOTU > LTUPQ) then
  IPQMX3 = LTUPQ/NOTU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
  IAD3S = 0
  call dDAFILE(LUHLF3,0,PQTU,IPQMX3,IAD3S)
end if
IAD3 = 0

!====================================================
! Second half transformation
!====================================================
!-------------------------------------------
! Second half transformation (AB,TU) coulomb
!  Always calculated.
!-------------------------------------------
! Loop over t,u pair,
if (icc /= 0) then
  ISPQRS = ((ISR**2-ISR)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISQ
  IAD2M(1,ISPQRS) = IAD13
  IPQTU = 0
  ITU = 0
  do NT=1,NOCR
    Num = NOCS
    if (ISR == ISS) Num = NT
    do NU=1,Num
      IPQTU = 1+NBPQ*ITU
      ITU = ITU+1
      ! Read sorted integral from disk
      if (IPQMX3 < NBPQ) then
        call RBuf_tra2(LUHLF3,PQTU,NBPQ,IPQMX3,NOTU,ITU,IPQTU,IAD3S)
      end if
      if (ISP == ISQ) then
        ! Square if necessary
        call SQUARE(PQTU(IPQTU),X2,1,NBP,NBP)
        ! (pq,TU) -> (Aq,TU)
        call DGEMM_('N','N',NBQ,NOP,NBP,One,X2,NBQ,CMO(LMOP),NBP,Zero,X1,NBQ)
        ! (Aq,TU) -> (AB,TU)
        call DGEMM_Tri('T','N',NOP,NOP,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,Zero,X2,NOP)
        IX2 = (NOP+NOP**2)/2
      else
        ! (pq,TU) -> (Aq,TU)
        call DGEMM_('N','N',NBQ,NOP,NBP,One,PQTU(IPQTU),NBQ,CMO(LMOP),NBP,Zero,X1,NBQ)
        ! (Aq,TU) -> (AB,TU)
        call DGEMM_('T','N',NOQ,NOP,NBQ,One,CMO(LMOQ),NBQ,X1,NBQ,Zero,X2,NOQ)
        IX2 = NOP*NOQ
      end if
      ! Store (AB,TU) of this t,u pair
      call GADSum(X2,IX2)
      call dDAFILE(LUINTM,1,X2,IX2,IAD13)
    end do
  end do
  ! End of loop over t,u pair
end if
!-----------------------------------------------------------------------
! Second half transformation Case 1 (AT,BU)
!  Always calculated. Both type 1 and 2. Case 2 (BU,AT) can be abandoned
!  since always ISP > ISR(equality can be removed by symmetry)
!-----------------------------------------------------------------------
! Loop over r,u pair
NOTU = NOCQ*NOCS
if (ISQ == ISS) NOTU = (NOCQ**2+NOCQ)/2
if (icxc1 /= 0) then
  LAR = LTUPQ/NOTU
  LR = LAR/NOP
  if (LR > NBR) LR = NBR
  LAR = NOP*LR
  IAD3S = 0
  call dDAFILE(LUHLF3,0,PQTU,LAR,IAD3S)
  IAD3 = 0
  IR = 0
  do NR=1,NBR
    IR = IR+1
    do NU=1,NOCS
      ! Square if necessary
      IRU = NBR*(NU-1)+NR
      IPQST = 1+NBPQ*(IRU-1)
      if (IPQMX2 < NBPQ) then
        call RBuf_tra2(LUHLF2,PQRU,NBPQ,IPQMX2,NORU,IRU,IPQST,IAD2S)
      end if
      if (ISP == ISQ) then
        call Square(PQRU(IPQST),X2,1,NBP,NBP)
      else
        call dcopy_(NBPQ,PQRU(IPQST),1,X2,1)
      end if
      ! if (ISQ == ISS) then triangular
      if (ISQ == ISS) then
        ! (pq,rU) -> (pT,rU)
        call DGEMM_('T','N',NBP,NOCQ-NU+1,NBQ,One,X2,NBQ,CMO(LMOQ2+NBQ*(NU-1)),NBQ,Zero,X1,NBP)
        ! (pT,rU) -> (AT,rU)
        call DGEMM_('T','N',NOCQ-NU+1,NOP,NBP,One,X1,NBP,CMO(LMOP2),NBP,Zero,X2,NOCQ-NU+1)
      else
        ! (pq,rU) -> (pT,rU)
        call DGEMM_('T','N',NBP,NOCQ,NBQ,One,X2,NBQ,CMO(LMOQ2),NBQ,Zero,X1,NBP)
        ! (pT,rU) -> (AT,rU)
        call DGEMM_('T','N',NOCQ,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOCQ)
      end if
      ! Store buffer
      if (IR > LR) then
        IR = 1
        call dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
      end if
      ! Sort
      NAT = 0
      do NA=1,NOP
        NTM = 1
        if (ISQ == ISS) NTM = NU
        do NT=NTM,NOCQ
          ITU = NOCS*(NT-1)+NU-1
          if (ISQ < ISS) ITU = NOCQ*(NU-1)+NT-1
          if (ISQ == ISS) ITU = (NT**2-NT)/2+NU-1
          NAT = NAT+1
          PQTU(LAR*ITU+NOP*(IR-1)+NA) = X2(NAT)
          !if ((isp == 7) .and. (isq == 1) .and. (isr == 6)) write(u6,'(f13.6)') pqtu(1)
        end do
      end do
      ! End of loop over r,u pair
    end do
  end do
  ! Store last buffer
  if (LR < NBR) then
    call dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
  end if
  ! Transform fourth index
  if (ISQ >= ISS) then
    ! Exchange type 1
    ISPQRS = ((ISQ**2-ISQ)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISR
    IAD2M(2,ISPQRS) = IAD13
    NTMAX = NOCQ
    NUMAX = NOCS
  else
    ! Exchange type 2
    ISPQRS = ((ISS**2-ISS)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISR
    IAD2M(3,ISPQRS) = IAD13
    NTMAX = NOCS
    NUMAX = NOCQ
  end if
  ! Loop over t,u pair, If (ISQ == ISS) loop should be triangle
  IST = 1-NOP*NBR
  KKTU = 0
  do NT=1,NTMAX
    NUM = NUMAX
    if (ISQ == ISS) NUM = NT
    do NU=1,NUM
      IST = IST+NOP*NBR
      KKTU = KKTU+1
      if (LR < NBR) then
        call RBuf_tra2(LUHLF3,PQTU,NBR*NOP,LAR,NOTU,KKTU,IST,IAD3S)
      end if
      ! (AT,rU) -> (AT,BU)
      call DGEMM_('T','T',NOR,NOP,NBR,One,CMO(LMOR2),NBR,PQTU(IST),NOP,Zero,X2,NOR)

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOP*NOR)
      call dDAFILE(LUINTM,1,X2,NOP*NOR,IAD13)
    end do
  end do
  ! End of loop over t,u pair
end if
!-----------------------------------------------------------------------
! Case 5 (TA,BU) and 6 (UB,AT)
!  Calculated if (ISP /= ISQ) .or. ((ISQ < ISR) .and. (ISP == ISR))
!  If ISQ >= ISR then Case 5, which alwats gives type 1 integrals
!  If (ISQ < ISR) .and. (ISP /= ISR)
!                then Case 6, which alwats gives type 2 integrals
!-----------------------------------------------------------------------
NOTU = NOCP*NOCS
! ISP == ISR in Case 6 should be skipped.
if (((ISQ >= ISR) .or. (ISP /= ISR)) .and. (ISP /= ISQ) .and. (icxc5 /= 0)) then
  LAR = LTUPQ/NOTU
  LR = LAR/NOQ
  if (LR > NBR) LR = NBR
  LAR = NOQ*LR
  IAD3S = 0
  call dDAFILE(LUHLF3,0,PQTU,LAR,IAD3S)
  IAD3 = 0
  IRU = 0
  IR = 0
  ! Loop over r,u pair
  do NR=1,NBR
    IR = IR+1
    do NU=1,NOCS
      ! Square is unnecessary
      IRU = NBR*(NU-1)+NR
      IPQST = 1+NBPQ*(IRU-1)
      if (IPQMX2 < NBPQ) then
        call RBuf_tra2(LUHLF2,PQRU,NBPQ,IPQMX2,NORU,IRU,IPQST,IAD2S)
      end if
      ! Always ISP > ISS i.e. s(T) > s(U)
      ! (pq,rU) -> (Tq,rU)
      call DGEMM_('N','N',NBQ,NOCP,NBP,One,PQRU(IPQST),NBQ,CMO(LMOP2),NBP,Zero,X1,NBQ)
      ! (Tq,rU) -> (TA,rU)
      call DGEMM_('T','N',NOCP,NOQ,NBQ,One,X1,NBQ,CMO(LMOQ2),NBQ,Zero,X2,NOCP)
      ! Store buffer
      if (IR > LR) then
        IR = 1
        call dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
      end if
      ! Sorting
      NAT = 0
      do NA=1,NOQ
        do NT=1,NOCP
          ITU = NOCS*(NT-1)+NU-1
          NAT = NAT+1
          PQTU(LAR*ITU+NOQ*(IR-1)+NA) = X2(NAT)
        end do
      end do
      ! End of loop over r,u pair
    end do
  end do
  ! Store last buffer
  if (LR < NBR) then
    call dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
  end if
  if (ISQ >= ISR) then
    ! Store(Only type1)
    ISPQRS = ((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISQ**2-ISQ)/2+ISR
    IAD2M(2,ISPQRS) = IAD13
  else if ((ISP /= ISR) .and. (ISR > ISQ)) then
    ! Store(Only type2)
    ISPQRS = ((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISR**2-ISR)/2+ISQ
    IAD2M(3,ISPQRS) = IAD13
  end if
  ! Loop over t,u pair
  IST = 1-NOQ*NBR
  KKTU = 0
  do NT=1,NOCP
    do NU=1,NOCS
      IST = IST+NOQ*NBR
      KKTU = KKTU+1
      if (LR < NBR) then
        call RBuf_tra2(LUHLF3,PQTU,NBR*NOQ,LAR,NOTU,KKTU,IST,IAD3S)
      end if
      ! (TA,rU) -> (TA,BU)
      if (ISQ >= ISR) then
        call DGEMM_('T','T',NOR,NOQ,NBR,One,CMO(LMOR),NBR,PQTU(IST),NOQ,Zero,X2,NOR)
      else if ((ISP /= ISR) .and. (ISR > ISQ)) then
        call DGEMM_('N','N',NOQ,NOR,NBR,One,PQTU(IST),NOQ,CMO(LMOR),NBR,Zero,X2,NOQ)
      end if

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOQ*NOR)
      call dDAFILE(LUINTM,1,X2,NOQ*NOR,IAD13)
    end do
  end do
  ! End of loop over t,u pair
end if

return

end subroutine tr2NsA2
