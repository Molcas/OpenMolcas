!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SORTB(BUFS,INDS,ACBDS,ACBDT,ISAB,BFACBD)
! SORTS INTEGRALS (AB/CD)
! FOR FIXED A,C ALL B,D

use mrci_global, only: ICH, IPASS, IRC, IROW, JJS, KBUFF1, LASTAD, LN, LSYM, Lu_60, Lu_80, LUTRA, MCHAIN, NBITM2, NCHN2, NORB, &
                       NSM, NSYM, NTIBUF, NVIR, NVIRP, NVIRT, TIBUF
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: BUFS(NBITM2,NCHN2)
integer(kind=iwp), intent(out) :: INDS(NBITM2+2,NCHN2)
real(kind=wp), intent(_OUT_) :: ACBDS(*), ACBDT(*), BFACBD(*)
integer(kind=iwp), intent(in) :: ISAB(*)
#include "tratoc.fh"
#include "warnings.h"
integer(kind=iwp) :: I, IAC, IACMAX, IACMIN, IAD16, IAD50, IADR, IBDS, IDISK, IFIN1, IFIN2, ILOOP, IN1, INB, INND, INPS, INPT, &
                     INS, INSB, INSOUT, INUMB, IOUT, IPOS, IST, IST1, IST2, ISTEP, ISYM, ITAIL, ITURN, JDISK, KK, LENGTH, M1, M2, &
                     M3, M4, N1, N2, N3, N4, NA, NAC, NB, NC, ND, NDMAX, NI, NJ, NK, NL, NOP, NOQ, NOR, NORB0(9), NORBP, NOS, &
                     NOVM, NOVST, NSAC, NSACL, NSC, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NSSM, NT, NTM, NU, NUMAX, NUMIN, NV, NVT, NX, &
                     NXM
real(kind=wp) :: FINI

NVT = IROW(NVIRT+1)
NOVST = LN*NVIRT+1
IAD16 = 0

NORB0(1) = 0
do I=1,NSYM
  NORB0(I+1) = NORB0(I)+NORB(I)
end do

INSOUT = 0
IACMAX = 0
do ISTEP=1,IPASS
  IAD50 = 0
  call iDAFILE(LUTRA,2,iTraToc,nTraToc,IAD50)
  IDISK = 0
  IACMIN = IACMAX+1
  IACMAX = IACMAX+NCHN2
  if (IACMAX > NVT) IACMAX = NVT
  if (IACMIN > IACMAX) cycle

  ! Initialize Buffer Counts and BackChain Links.
  INDS(NBITM2+1,:) = 0
  INDS(NBITM2+2,:) = -1

  ! Loop over symmetry blocks of all-virtual integrals.
  do NSP=1,NSYM
    NOP = NORB(NSP)
    do NSQ=1,NSP
      NSPQ = MUL(NSP,NSQ)
      NOQ = NORB(NSQ)
      do NSR=1,NSP
        NSPQR = MUL(NSPQ,NSR)
        NOR = NORB(NSR)
        NSSM = NSR
        if (NSR == NSP) NSSM = NSQ
        do NSS=1,NSSM
          if (NSS /= NSPQR) cycle
          NOS = NORB(NSS)
          NORBP = NOP*NOQ*NOR*NOS
          if (NORBP == 0) cycle

          call dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)

          ! Loop over index quadruples in this symm block
          IOUT = 0
          do NV=1,NOR
            NXM = NOS
            if (NSR == NSS) NXM = NV
            do NX=1,NXM
              NTM = 1
              if (NSP == NSR) NTM = NV
              do NT=NTM,NOP
                NUMIN = 1
                if ((NSP == NSR) .and. (NT == NV)) NUMIN = NX
                NUMAX = NOQ
                if (NSP == NSQ) NUMAX = NT
                do NU=NUMIN,NUMAX

                  ! MO integral value is made accessable at TIBUF(IOUT)
                  IOUT = IOUT+1
                  if (IOUT > NTIBUF) then
                    call dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)
                    IOUT = 1
                  end if

                  ! M1..M4 seqential number of orbitals.
                  M1 = ICH(NORB0(NSP)+NT)
                  M2 = ICH(NORB0(NSQ)+NU)
                  M3 = ICH(NORB0(NSR)+NV)
                  M4 = ICH(NORB0(NSS)+NX)
                  if ((M1 <= LN) .or. (M2 <= LN)) cycle
                  if ((M3 <= LN) .or. (M4 <= LN)) cycle

                  ! Permute orbital indices to canonical order
                  ! and put integral value in FINI
                  N1 = max(M1,M2)
                  N2 = min(M1,M2)
                  N3 = max(M3,M4)
                  N4 = min(M3,M4)
                  NI = N1
                  NJ = N2
                  NK = N3
                  NL = N4
                  if (NI <= NK) then
                    if (NI /= NK) then
                      NI = N3
                      NJ = N4
                      NK = N1
                      NL = N2
                    else if (NJ > NL) then
                      NL = N2
                      NJ = N4
                    end if
                  end if
                  FINI = TIBUF(IOUT)

                  ! Compute virtual indices.
                  NA = NI-LN
                  NB = NJ-LN
                  NC = NK-LN
                  ND = NL-LN
                  ITURN = 0
                  if ((NA == NB) .and. (NC == ND)) cycle
                  do
                    IAC = IROW(NA)+NC
                    if ((IAC >= IACMIN) .and. (IAC <= IACMAX)) then
                      if ((NA == NC) .and. (NB == ND)) FINI = FINI/2
                      NAC = IAC-IACMIN+1
                      IPOS = INDS(NBITM2+1,NAC)+1
                      INDS(NBITM2+1,NAC) = IPOS
                      INDS(IPOS,NAC) = NB+2**8*ND
                      BUFS(IPOS,NAC) = FINI
                      if (IPOS >= NBITM2) then
                        ! Save this buffer if filled up.
                        JDISK = IDISK
                        call iDAFILE(Lu_60,1,INDS(1,NAC),NBITM2+2,IDISK)
                        call dDAFILE(Lu_60,1,BUFS(1,NAC),NBITM2,IDISK)
                        INDS(NBITM2+1,NAC) = 0
                        INDS(NBITM2+2,NAC) = JDISK
                      end if
                    end if

                    if (ITURN == 1) exit
                    if ((NA == NC) .and. (NB == ND)) exit
                    if ((NA == NB) .or. (NC == ND)) exit
                    ITURN = 1
                    NC = NL-LN
                    ND = NK-LN
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  ! EMPTY LAST BUFFERS
  NOVM = IACMAX-IACMIN+1
  if ((NOVST+IACMIN-1+NOVM) > MCHAIN) then
    write(u6,*) 'SORTB Error: NOVST+IACMIN-1+NOVM > MCHAIN'
    write(u6,*) 'NOVST =',NOVST
    write(u6,*) 'IACMIN=',IACMIN
    write(u6,*) 'NOVM  =',NOVM
    write(u6,*) 'MCHAIN=',MCHAIN
    write(u6,*) '  (See code).'
    call QUIT(_RC_GENERAL_ERROR_)
  end if
  do I=1,NOVM
    JDISK = IDISK
    call iDAFILE(Lu_60,1,INDS(1,I),NBITM2+2,IDISK)
    call dDAFILE(Lu_60,1,BUFS(1,I),NBITM2,IDISK)
    LASTAD(NOVST+IACMIN-1+I) = JDISK
  end do
  do ISYM=1,NSYM
    IST1 = IRC(3)+JJS(ISYM+9)+1
    IFIN1 = IRC(3)+JJS(ISYM+10)
    INPS = IFIN1-IST1+1
    IST2 = IRC(2)+JJS(ISYM)+1
    IFIN2 = IRC(2)+JJS(ISYM+1)
    INPT = IFIN2-IST2+1
    ITAIL = INPS+INPT
    if (ITAIL == 0) cycle
    IN1 = -NVIRT
    do NA=1,NVIRT
      IN1 = IN1+NVIRT
      do NC=1,NA
        IAC = IROW(NA)+NC
        if (IAC < IACMIN) cycle
        if (IAC > IACMAX) cycle
        if (NA == 1) cycle
        NSAC = MUL(NSM(LN+NA),NSM(LN+NC))
        NSACL = MUL(NSAC,LSYM)
        if (NSACL /= ISYM) cycle
        NSC = NSM(LN+NC)
        NDMAX = NVIRP(NSC)+NVIR(NSC)
        if (NDMAX > NA) NDMAX = NA
        INS = ISAB(NA+(NDMAX-1)*NVIRT)
        ACBDS(1:INS) = Zero
        ACBDT(1:INS) = Zero
        IADR = LASTAD(NOVST+IAC)
        do
          call iDAFILE(Lu_60,2,INDS,NBITM2+2,IADR)
          call dDAFILE(Lu_60,2,BUFS,NBITM2,IADR)
          LENGTH = INDS(NBITM2+1,1)
          IADR = INDS(NBITM2+2,1)
          do KK=1,LENGTH
            INND = INDS(KK,1)
            NB = ibits(INND,0,8)
            ND = ibits(INND,8,8)

            IBDS = ISAB(NB+(ND-1)*NVIRT)
            ACBDS(IBDS) = ACBDS(IBDS)+BUFS(KK,1)
            if (NB > ND) ACBDT(IBDS) = ACBDT(IBDS)+BUFS(KK,1)
            if (NB < ND) ACBDT(IBDS) = ACBDT(IBDS)-BUFS(KK,1)
          end do
          if (IADR == -1) exit
        end do
        ILOOP = 0
        do
          INSB = INS
          do
            INB = KBUFF1-INSOUT
            INUMB = INSB
            if (INSB > INB) INUMB = INB
            IST = INS-INSB+1
            if (ILOOP == 0) call DCOPY_(INUMB,ACBDS(IST),1,BFACBD(INSOUT+1),1)
            if (ILOOP == 1) call DCOPY_(INUMB,ACBDT(IST),1,BFACBD(INSOUT+1),1)
            INSOUT = INSOUT+INUMB
            if (INSOUT >= KBUFF1) then
              call dDAFILE(Lu_80,1,BFACBD,KBUFF1,IAD16)
              INSOUT = 0
            end if
            INSB = INSB-INUMB
            if (INSB <= 0) exit
          end do
          ILOOP = ILOOP+1
          if (ILOOP /= 1) exit
        end do
      end do
    end do
  end do
end do
! EMPTY LAST BUFFER
if (INSOUT /= 0) call dDAFILE(Lu_80,1,BFACBD,KBUFF1,IAD16)

return

end subroutine SORTB
