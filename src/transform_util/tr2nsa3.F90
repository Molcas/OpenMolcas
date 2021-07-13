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

subroutine tr2nsa3(CMO,X1,nX1,X2,nX2,pqUs,npqUS,pqrU,npqrU)
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
integer(kind=iwp), intent(in) :: nX1, nX2, npqUS, npqrU
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(out) :: X1(nX1), X2(nX2)
real(kind=wp), intent(inout) :: pqUS(npqUS), pqrU(npqrU)
integer(kind=iwp) :: IAD1S, IAD2, IAD2S, icxc3, icxc7, IPQMX1, IPQMX2, IPQST, IRU, IS, ISPQRS, IST, ITU, IUS, KKTU, LAS, LS, NA, &
                     NAT, NORU, NOTU, NOUS, NS, NSYMP, NT, NTM, NTMAX, NU, Num, NUMAX
#include "trafo.fh"
#include "intgrl.fh"

NSYMP = NSYM*(NSYM+1)/2
NOTU = NOCR*NOCS
if (ISR == ISS) NOTU = (NOCR**2+NOCR)/2
NOUS = NOCR*NBS
NORU = NBR*NOCS
icxc3 = NOP*NOCQ*NOCR*NOS
icxc7 = NOCP*NOQ*NOCR*NOS

! Check for in core or out of core transformation

! 1. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|US) ON UNIT LUHLF1
IPQMX1 = NBPQ
if (NBPQ*NOUS > LURPQ) then
  IPQMX1 = LURPQ/NOUS
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|US)',IPQMX1
  IAD1S = 0
  call dDAFILE(LUHLF1,0,PQUS,IPQMX1,IAD1S)
end if
! 2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
IPQMX2 = NBPQ
if (NBPQ*NORU > LRUPQ) then
  IPQMX2 = LRUPQ/NORU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
  IAD2S = 0
  call dDAFILE(LUHLF2,0,PQRU,IPQMX2,IAD2S)
end if
IAD2 = 0
! 3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
!IPQMX3 = NBPQ
!if (NBPQ*NOTU > LTUPQ) then
!  IPQMX3 = LTUPQ/NOTU
!  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
!  IAD3S = 0
!  call dDAFILE(LUHLF3,0,PQTU,IPQMX3,IAD3S)
!end if
!IAD3 = 0
!IOUT3 = 0

!-----------------------------------------------------------------------
! Second half transformation Case 3 (AT,UB)
!  Calculated if ISR /= ISS. Both type 1 and 2. ISQ /= ISR, since
!  equality makes integral to be zero symmetrically.
!  Case 4 (BU,TA) need not be calculated, since always ISP > ISS.
!-----------------------------------------------------------------------
NOTU = NOCQ*NOCR
if ((ISR /= ISS) .and. (icxc3 /= 0)) then
  LAS = (LRUPQ+LTUPQ)/NOTU
  LS = LAS/NOP
  if (LS > NBS) LS = NBS
  LAS = NOP*LS
  IAD2S = 0
  call dDAFILE(LUHLF2,0,PQRU,LAS,IAD2S)
  IAD2 = 0
  ! Loop over u,s pair
  IS = 0
  do NS=1,NBS
    IS = IS+1
    do NU=1,NOCR
      ! Square if necessary
      IUS = NBS*(NU-1)+NS
      IPQST = 1+NBPQ*(IUS-1)
      if (IPQMX1 < NBPQ) then
        call RBuf_tra2(LUHLF1,PQUS,NBPQ,IPQMX1,NOUS,IUS,IPQST,IAD1S)
      end if
      if (ISP == ISQ) then
        call Square(PQUS(IPQST),X2,1,NBP,NBP)
      else
        call dcopy_(NBPQ,PQUS(IPQST),1,X2,1)
      end if
      ! Always ISQ /= ISR  i.e. s(T) /= s(U)
      ! (pq,Us) -> (pT,Us)
      call DGEMM_('T','N',NBP,NOCQ,NBQ,One,X2,NBQ,CMO(LMOQ2),NBQ,Zero,X1,NBP)
      ! (pT,Us) -> (AT,Us)
      call DGEMM_('T','N',NOCQ,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOCQ)
      ! Store buffer
      if (IS > LS) then
        IS = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
        !vv end do
        call dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
      end if
      ! Sort
      NAT = 0
      do NA=1,NOP
        do NT=1,NOCQ
          ITU = NOCR*(NT-1)+NU-1
          if (ISQ < ISR) ITU = NOCQ*(NU-1)+NT-1
          NAT = NAT+1
          PQRU(LAS*ITU+NOP*(IS-1)+NA) = X2(NAT)
        end do
      end do
      ! End of loop over u,s pair
    end do
  end do
  ! Store last buffer
  if (LS < NBS) then
    !vv do I=1,NOTU
    !vv   call dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
    !vv end do
    call dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
  end if

  ! Trasform last index

  if (ISQ >= ISR) then
    ! Store(type 1)
    ISPQRS = ((ISQ**2-ISQ)/2+ISR-1)*NSYMP+(ISP**2-ISP)/2+ISS
    IAD2M(2,ISPQRS) = IAD13
    NTMAX = NOCQ
    NUMAX = NOCR
  else
    ! Store(type 2)
    ISPQRS = ((ISR**2-ISR)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISS
    IAD2M(3,ISPQRS) = IAD13
    NTMAX = NOCR
    NUMAX = NOCQ
  end if
  ! Loop over t,u pair, always ISQ /= ISR
  IST = 1-NBS*NOP
  KKTU = 0
  do NT=1,NTMAX
    do NU=1,NUMAX
      IST = IST+NBS*NOP
      KKTU = KKTU+1
      if (LS < NBS) then
        call RBuf_tra2(LUHLF2,PQRU,NBS*NOP,LAS,NOTU,KKTU,IST,IAD2S)
      end if
      ! (AT,Us) -> (AT,UB)
      call DGEMM_('T','T',NOS,NOP,NBS,One,CMO(LMOS2),NBS,PQRU(IST),NOP,Zero,X2,NOS)

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOP*NOS)
      call dDAFILE(LUINTM,1,X2,NOP*NOS,IAD13)
    end do
  end do
  ! End of loop over t,u pair
end if
!-----------------------------------------------------------------------
! Case 7 (TA,UB) and 8 (UB,TA)
!  Calculated if (ISP /= ISQ) .and. (ISR /= ISS)
!   Case 7, if ISQ >= ISS, always type 1
!   Case 8, if (ISQ < ISS) .and. (ISP /= ISR), always type 2
!-----------------------------------------------------------------------
NOTU = NOCP*NOCR
if (((ISS <= ISQ) .or. (ISP /= ISR)) .and. (ISP /= ISQ) .and. (ISR /= ISS) .and. (icxc7 /= 0)) then
  LAS = (LRUPQ+LTUPQ)/NOTU
  LS = LAS/NOQ
  if (LS > NBS) LS = NBS
  LAS = NOQ*LS
  IAD2S = 0
  call dDAFILE(LUHLF2,0,PQRU,LAS,IAD2S)
  IAD2 = 0
  IRU = 0
  ! Loop over u,s pair
  IS = 0
  do NS=1,NBS
    IS = IS+1
    do NU=1,NOCR
      IRU = NBS*(NU-1)+NS
      IPQST = 1+NBPQ*(IRU-1)
      if (IPQMX1 < NBPQ) then
        call RBuf_tra2(LUHLF1,PQUS,NBPQ,IPQMX1,NOUS,IRU,IPQST,IAD1S)
      end if
      ! Always ISP > ISQ
      ! Square unnecessary
      if (ISP == ISR) then
        ! (pq,Us) -> (Tq,Us)
        call DGEMM_('N','N',NBQ,NOCP-NU+1,NBP,One,PQUS(IPQST),NBQ,CMO(LMOP2+NBP*(NU-1)),NBP,Zero,X1,NBQ)
        ! (Tq,Us) -> (TA,Us)
        call DGEMM_('T','N',NOCP-NU+1,NOQ,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,Zero,X2,NOCP-NU+1)
      else
        ! (pq,Us) -> (Tq,Us)
        call DGEMM_('N','N',NBQ,NOCP,NBP,One,PQUS(IPQST),NBQ,CMO(LMOP2),NBP,Zero,X1,NBQ)
        ! (Tq,Us) -> (TA,Us)
        call DGEMM_('T','N',NOCP,NOQ,NBQ,One,X1,NBQ,CMO(LMOQ),NBQ,Zero,X2,NOCP)
      end if
      ! Store buffer
      if (IS > LS) then
        IS = 1
        !vv do I=1,NOTU
        !vv   call dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
        !vv end do
        call dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
      end if
      ! Sorting
      ! Note: LAS is supposed to be small enough that (NOCP*NOCR)*LAS is
      ! at most = (LRUPQ+LTUPQ). NU can be as large as NOCR, so ITU
      ! can become NOCR*NOCP-1.
      ! Thus (LAS*ITU+NOQ*(IS-1)+NA can become as large as:
      ! LAS*(NOCP*NOCR-1)+NOQ*(NBS-1)+NOQ=LAS*(NOCP*NOCR-1)+NOQ*NBS
      ! which may become nearly as large as (LRUPQ+LTUPQ)+NOQ*NBS.
      ! But when called, only the size LRUPQ has been reserved for
      ! use (tractl:LW6=LW5+LRUPQ) by the array PQRU!
      ! This may be deliberate, if it is intended that PQRU is
      ! overlaying i.e. extended above the LW6 address which is
      ! then assumed not to be used any longer....
      NAT = 0
      do NA=1,NOQ
        NTM = 1
        if (ISP == ISR) NTM = NU
        do NT=NTM,NOCP
          ITU = NOCR*(NT-1)+NU-1
          if (ISP == ISR) ITU = (NT*NT-NT)/2+NU-1
          NAT = NAT+1
          PQRU(LAS*ITU+NOQ*(IS-1)+NA) = X2(NAT)
        end do
      end do
      ! End of Loop over u,s pair
    end do
  end do
  ! Store the last buffer
  if (LS < NBS) then
    !vv do I=1,NOTU
    !vv   call dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
    !vv end do
    call dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
  end if
  if (ISQ >= ISS) then
    ! Store(Only type1)
    ISPQRS = ((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
    IAD2M(2,ISPQRS) = IAD13
  else if ((ISP /= ISR) .and. (ISQ < ISS)) then
    ! Store(Only type2)
    ISPQRS = ((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISS**2-ISS)/2+ISQ
    IAD2M(3,ISPQRS) = IAD13
  end if
  IST = 1-NOQ*NBS
  KKTU = 0
  ! Loop over t,u pair, If (ISP == ISR) then loop should be triangle
  do NT=1,NOCP
    Num = NOCR
    if (ISP == ISR) Num = NT
    do NU=1,Num
      ! (TA,Us) -> (TA,UB)
      IST = IST+NOQ*NBS
      KKTU = KKTU+1
      if (LS < NBS) then
        call RBuf_tra2(LUHLF2,PQRU,NBS*NOQ,LAS,NOTU,KKTU,IST,IAD2S)
      end if
      if (ISQ >= ISS) then
        call DGEMM_('T','T',NOS,NOQ,NBS,One,CMO(LMOS),NBS,PQRU(IST),NOQ,Zero,X2,NOS)
      else if ((ISP /= ISR) .and. (ISS > ISQ)) then
        call DGEMM_('N','N',NOQ,NOS,NBS,One,PQRU(IST),NOQ,CMO(LMOS),NBS,Zero,X2,NOQ)
      end if

      ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

      call GADSum(X2,NOQ*NOS)
      call dDAFILE(LUINTM,1,X2,NOQ*NOS,IAD13)
    end do
  end do
  ! End of Loop over t,u pair
end if

return

end subroutine tr2nsa3
