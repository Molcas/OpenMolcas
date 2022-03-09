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
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine SORTA_CPF(BUFOUT,INDOUT,ICAD,IBUFL,TIBUF,ISAB,BUFBI,INDBI,BIAC,BICA,NINTGR)
! SORTS INTEGRALS (AB/CI)
! FOR FIXED B,I ALL A,C
! FIRST CHAIN FOR IJKL

use cpf_global, only: IADABCI, ICH, IFIRST, IPRINT, IROW, KBUF, KBUFF1, LASTAD, LN, Lu_CIGuga, Lu_TiABCI, Lu_TiABIJ, Lu_TraInt, &
                      MADR, NNS, NORB, NSM, NSYM, NTIBUF, NVIRT
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoI

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: BUFOUT(*), TIBUF(NTIBUF), BUFBI(*), BIAC(*), BICA(*)
integer(kind=iwp), intent(_OUT_) :: INDOUT(*), ICAD(*), IBUFL(*), INDBI(*), NINTGR
integer(kind=iwp), intent(in) :: ISAB(*)
#include "tratoc.fh"
integer(kind=iwp) :: I, IACS, IAD15, IAD50, IADD10, IADR, IBUFIJ, ICHK, ICP, ICPP, ICQ, ID, IDISK, IDIV, IIJ, IIN, IJ, IJKL, ILEN, &
                     ILOOP, INND, INS, INSOUT, IOUT, IREC, ITURN, JDISK, KBUF0, KBUF1, KBUF2, KK, KKBUF0, KKBUF1, KKBUF2, KL, &
                     LENGTH, M1, M2, M3, M4, N1, N2, N3, N4, NA, NAC, NAT, NB, NC, NI, NIB, NJ, NK, NL, NOP, NOQ, NOR, NORB0(9), &
                     NORBP, NOS, NOV, NOVST, NSAVE, NSIB, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NSSM, NT, NTM, NU, NUMAX, NUMIN, NV, &
                     NX, NXM
real(kind=wp) :: FINI
logical(kind=iwp) :: Skip

call COUNT_CPF(NINTGR,NSYM,NORB,MUL)
if (IPRINT >= 2) then
  write(6,*) ' NUMBER OF TWO-ELECTRON INTEGRALS:',NINTGR
end if
IAD50 = 0
call iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
KKBUF0 = (RTOI*(KBUFF1+2)-2)/(RTOI+1)
KKBUF1 = RTOI*KKBUF0+KKBUF0+1
KKBUF2 = KKBUF1+1
NOV = LN*NVIRT+1
if (IFIRST /= 0) NOV = 1
IDISK = 0
KBUF0 = RTOI*KBUF
KBUF1 = KBUF0+KBUF+1
KBUF2 = KBUF1+1
IDIV = RTOI
ID = 0
do IREC=1,NOV
  IBUFL(IREC) = 0
  ICAD(IREC) = ID
  INDOUT(ID+KBUF2) = -1
  ID = ID+KBUF2
end do
NORB0(1) = 0
do I=1,NSYM
  NORB0(I+1) = NORB0(I)+NORB(I)
end do

! TWO-ELECTRON INTEGRALS

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
        call dDAFILE(Lu_TraInt,2,TIBUF,NTIBUF,IAD50)
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
              loop1: do NU=NUMIN,NUMAX
                IOUT = IOUT+1
                if (IOUT > NTIBUF) then
                  call dDAFILE(Lu_TraInt,2,TIBUF,NTIBUF,IAD50)
                  IOUT = 1
                end if
                M1 = ICH(NORB0(NSP)+NT)
                M2 = ICH(NORB0(NSQ)+NU)
                M3 = ICH(NORB0(NSR)+NV)
                M4 = ICH(NORB0(NSS)+NX)
                if ((M1 <= 0) .or. (M2 <= 0) .or. (M3 <= 0) .or. (M4 <= 0)) cycle loop1
                ! ORDER THESE INDICES CANONICALLY
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
                  else if (NJ <= NL) then
                    NL = N2
                    NJ = N4
                  end if
                end if
                FINI = TIBUF(IOUT)
                if (abs(FINI) < 1.0e-9_wp) cycle loop1
                if (NI > LN) then
                  if (NK <= LN) cycle loop1
                  if (NJ > LN) then
                    if (NL > LN) cycle loop1
                    if (IFIRST /= 0) cycle loop1

                    ! ABCI

                    NA = NI-LN
                    NB = NJ-LN
                    NC = NK-LN
                    NI = NL
                    Skip = .false.
                  else
                    Skip = .true.
                  end if
                  do
                    if (Skip) then
                      Skip = .false.
                    else
                      ITURN = 0
                      do
                        NIB = (NI-1)*NVIRT+NB+1
                        IBUFL(NIB) = IBUFL(NIB)+1
                        ICQ = ICAD(NIB)
                        ICP = ICQ/IDIV+IBUFL(NIB)
                        BUFOUT(ICP) = FINI
                        ICPP = ICQ+KBUF0+IBUFL(NIB)
                        INDOUT(ICPP) = (NA-1)*NVIRT+NC
                        if (IBUFL(NIB) >= KBUF) then
                          INDOUT(ICQ+KBUF1) = KBUF
                          JDISK = IDISK
                          call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,IDISK)
                          INDOUT(ICQ+KBUF2) = JDISK
                          IBUFL(NIB) = 0
                        end if
                        if ((ITURN == 1) .or. (NA == NB)) cycle loop1
                        ITURN = 1
                        NAT = NA
                        NA = NB
                        NB = NAT
                      end do
                    end if
                    if ((NJ <= 0) .or. (NL <= LN)) cycle loop1

                    ! CIAB

                    if (IFIRST /= 0) cycle loop1
                    NA = NK-LN
                    NB = NL-LN
                    NC = NI-LN
                    NI = NJ
                  end do
                end if

                ! IJKL

                IIJ = IROW(NI)+NJ
                KL = IROW(NK)+NL
                IJKL = IIJ*(IIJ-1)/2+KL
                IJ = 1
                IBUFL(IJ) = IBUFL(IJ)+1
                ICQ = ICAD(IJ)
                ICP = ICQ/IDIV+IBUFL(IJ)
                BUFOUT(ICP) = FINI
                ICPP = ICQ+KBUF0+IBUFL(IJ)
                INDOUT(ICPP) = IJKL
                if (IBUFL(IJ) >= KBUF) then
                  INDOUT(ICQ+KBUF1) = KBUF
                  JDISK = IDISK
                  call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,IDISK)
                  INDOUT(ICQ+KBUF2) = JDISK
                  IBUFL(IJ) = 0
                end if
              end do loop1
            end do
          end do
        end do
      end do
    end do
  end do
end do
! EMPTY LAST BUFFERS
!FUE Start of insertion
if (NOV > MADR) then
  write(u6,*) 'SORTA_CPF Error: NOV > MADR (See code).'
  call Abend()
end if
!FUE End of insertion
do I=1,NOV
  ICQ = ICAD(I)
  INDOUT(ICQ+KBUF1) = IBUFL(I)
  JDISK = IDISK
  call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,IDISK)
  LASTAD(I) = JDISK
end do

! IJKL

IDISK = 0
IBUFIJ = 0
INDBI(KKBUF2) = -1
IADR = LASTAD(1)
do
  call iDAFILE(Lu_TiABIJ,2,INDOUT,KBUF2,IADR)
  LENGTH = INDOUT(KBUF1)
  IADR = INDOUT(KBUF2)
  do I=1,LENGTH
    IBUFIJ = IBUFIJ+1
    BUFBI(IBUFIJ) = BUFOUT(I)
    INDBI(RTOI*KKBUF0+IBUFIJ) = INDOUT(KBUF0+I)
    if (IBUFIJ >= KKBUF0) then
      INDBI(KKBUF1) = KKBUF0
      JDISK = IDISK
      call iDAFILE(Lu_TiABCI,1,INDBI,KKBUF2,IDISK)
      INDBI(KKBUF2) = JDISK
      IBUFIJ = 0
    end if
  end do
  if (IADR == -1) exit
end do
! EMPTY LAST BUFFER
INDBI(KKBUF1) = IBUFIJ
JDISK = IDISK
call iDAFILE(Lu_TiABCI,1,INDBI,KKBUF2,IDISK)
LASTAD(1) = JDISK

! ABCI

ICHK = 0
IAD15 = IDISK
IADABCI = IAD15
INSOUT = 0
NOVST = 1
IADD10 = IAD10(4)
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
IIN = 2
NSAVE = ICOP1(IIN)
do
  NI = NSAVE
  IOUT = 0
  Skip = .false.
  do
    IIN = IIN+1
    if (IIN > ILEN) then
      call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      ILEN = ICOP1(nCOP+1)
      if (ILEN <= 0) then
        Skip = .true.
        exit
      end if
      IIN = 1
    end if
    if (ICHK /= 0) exit
    if (ICOP1(IIN) /= 0) then
      IOUT = IOUT+1
    else
      ICHK = 1
    end if
  end do
  if (.not. Skip) then
    ICHK = 0
    NSAVE = ICOP1(IIN)
  end if
  NIB = (NI-1)*NVIRT+NOVST
  do NB=1,NVIRT
    NSIB = MUL(NSM(LN+NB),NSM(NI))
    INS = NNS(NSIB)
    BIAC(1:INS) = Zero
    BICA(1:INS) = Zero
    NIB = NIB+1
    IADR = LASTAD(NIB)
    do
      call iDAFILE(Lu_TiABIJ,2,INDOUT,KBUF2,IADR)
      LENGTH = INDOUT(KBUF1)
      IADR = INDOUT(KBUF2)
      do KK=1,LENGTH
        INND = INDOUT(KBUF0+KK)
        NA = (INND-1)/NVIRT+1
        NC = INND-(NA-1)*NVIRT
        NAC = (NA-1)*NVIRT+NC
        IACS = ISAB(NAC)
        BIAC(IACS) = BIAC(IACS)+BUFOUT(KK)
        if (NA > NC) BICA(IACS) = BICA(IACS)-BUFOUT(KK)
        if (NA < NC) BICA(IACS) = BICA(IACS)+BUFOUT(KK)
      end do
      if (IADR == -1) exit
    end do
    ILOOP = 0
    do
      do I=1,INS
        INSOUT = INSOUT+1
        if (ILOOP == 0) BUFBI(INSOUT) = BIAC(I)
        if (ILOOP == 1) BUFBI(INSOUT) = BICA(I)
        if (INSOUT >= KBUFF1) then
          call dDAFILE(Lu_TiABCI,1,BUFBI,KBUFF1,IAD15)
          INSOUT = 0
        end if
      end do
      ILOOP = ILOOP+1
      if (ILOOP /= 1) exit
    end do
  end do
  if (ILEN < 0) exit
end do
! EMPTY LAST BUFFER
if (INSOUT /= 0) call dDAFILE(Lu_TiABCI,1,BUFBI,KBUFF1,IAD15)

return

end subroutine SORTA_CPF
