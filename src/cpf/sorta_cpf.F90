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

use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoI

implicit none
#include "cpfmcpf.fh"
real(kind=wp) :: BUFOUT(*), TIBUF(NTIBUF), BUFBI(*), BIAC(*), BICA(*)
integer(kind=iwp) :: INDOUT(*), ICAD(*), IBUFL(*), ISAB(*), INDBI(*), NINTGR
#include "files_cpf.fh"
integer(kind=iwp) :: I, IACS, IAD15, IAD50, IADR, IBUFIJ, ICHK, ICP, ICPP, ICQ, ID, IDISK, IDIV, IIJ, IIN, IJ, IJKL, ILEN, ILOOP, &
                     INND, INS, INSOUT, IOUT, IREC, ITURN, JDISK, KBUF0, KBUF1, KBUF2, KK, KKBUF0, KKBUF1, KKBUF2, KL, LENGTH, M1, &
                     M2, M3, M4, N1, N2, N3, N4, NA, NAC, NAT, NB, NC, NI, NIB, NJ, NK, NL, NOP, NOQ, NOR, NORB0(9), NORBP, NOS, &
                     NOV, NOVST, NSAVE, NSIB, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NSSM, NT, NTM, NU, NUMAX, NUMIN, NV, NX, NXM
real(kind=wp) :: FINI

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
        if (NSS /= NSPQR) GO TO 310
        NOS = NORB(NSS)
        NORBP = NOP*NOQ*NOR*NOS
        if (NORBP == 0) GO TO 310
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
              do NU=NUMIN,NUMAX
                IOUT = IOUT+1
                if (IOUT > NTIBUF) then
                  call dDAFILE(Lu_TraInt,2,TIBUF,NTIBUF,IAD50)
                  IOUT = 1
                end if
                M1 = ICH(NORB0(NSP)+NT)
                M2 = ICH(NORB0(NSQ)+NU)
                M3 = ICH(NORB0(NSR)+NV)
                M4 = ICH(NORB0(NSS)+NX)
                if ((M1 <= 0) .or. (M2 <= 0)) GO TO 306
                if ((M3 <= 0) .or. (M4 <= 0)) GO TO 306
                ! ORDER THESE INDICES CANONICALLY
                N1 = M1
                N2 = M2
                if (M1 > M2) GO TO 11
                N1 = M2
                N2 = M1
11              N3 = M3
                N4 = M4
                if (M3 > M4) GO TO 12
                N3 = M4
                N4 = M3
12              NI = N1
                NJ = N2
                NK = N3
                NL = N4
                if (NI > NK) GO TO 502
                if (NI == NK) GO TO 14
                NI = N3
                NJ = N4
                NK = N1
                NL = N2
                GO TO 502
14              if (NJ > NL) GO TO 502
                NL = N2
                NJ = N4
502             FINI = TIBUF(IOUT)
                if (abs(FINI) < 1.0e-9_wp) GO TO 306
                if (NI <= LN) GO TO 109
                if (NK <= LN) GO TO 306
                if (NJ <= LN) GO TO 42
                if (NL > LN) GO TO 306
                if (IFIRST /= 0) GO TO 306

                ! ABCI

                NA = NI-LN
                NB = NJ-LN
                NC = NK-LN
                NI = NL
108             ITURN = 0
107             NIB = (NI-1)*NVIRT+NB+1
                IBUFL(NIB) = IBUFL(NIB)+1
                ICQ = ICAD(NIB)
                ICP = ICQ/IDIV+IBUFL(NIB)
                BUFOUT(ICP) = FINI
                ICPP = ICQ+KBUF0+IBUFL(NIB)
                INDOUT(ICPP) = (NA-1)*NVIRT+NC
                if (IBUFL(NIB) < KBUF) GO TO 106
                INDOUT(ICQ+KBUF1) = KBUF
                JDISK = IDISK
                call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,IDISK)
                INDOUT(ICQ+KBUF2) = JDISK
                IBUFL(NIB) = 0
106             if ((ITURN == 1) .or. (NA == NB)) GO TO 306
                ITURN = 1
                NAT = NA
                NA = NB
                NB = NAT
                GO TO 107
42              if ((NJ <= 0) .or. (NL <= LN)) GO TO 306

                ! CIAB

                if (IFIRST /= 0) GO TO 306
                NA = NK-LN
                NB = NL-LN
                NC = NI-LN
                NI = NJ
                GO TO 108

                ! IJKL

109             IIJ = IROW(NI)+NJ
                KL = IROW(NK)+NL
                IJKL = IIJ*(IIJ-1)/2+KL
                IJ = 1
                IBUFL(IJ) = IBUFL(IJ)+1
                ICQ = ICAD(IJ)
                ICP = ICQ/IDIV+IBUFL(IJ)
                BUFOUT(ICP) = FINI
                ICPP = ICQ+KBUF0+IBUFL(IJ)
                INDOUT(ICPP) = IJKL
                if (IBUFL(IJ) < KBUF) GO TO 306
                INDOUT(ICQ+KBUF1) = KBUF
                JDISK = IDISK
                call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,IDISK)
                INDOUT(ICQ+KBUF2) = JDISK
                IBUFL(IJ) = 0
306             continue
              end do
            end do
          end do
        end do
310     continue
      end do
    end do
  end do
end do
! EMPTY LAST BUFFERS
!FUE Start of insertion
if (NOV > mAdr) then
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
201 call iDAFILE(Lu_TiABIJ,2,INDOUT,KBUF2,IADR)
LENGTH = INDOUT(KBUF1)
IADR = INDOUT(KBUF2)
if (LENGTH == 0) GO TO 209
do I=1,LENGTH
  IBUFIJ = IBUFIJ+1
  BUFBI(IBUFIJ) = BUFOUT(I)
  INDBI(RTOI*KKBUF0+IBUFIJ) = INDOUT(KBUF0+I)
  if (IBUFIJ < KKBUF0) GO TO 202
  INDBI(KKBUF1) = KKBUF0
  JDISK = IDISK
  call iDAFILE(Lu_TiABCI,1,INDBI,KKBUF2,IDISK)
  INDBI(KKBUF2) = JDISK
  IBUFIJ = 0
202 continue
end do
209 if (IADR /= -1) GO TO 201
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
100 NI = NSAVE
IOUT = 0
110 IIN = IIN+1
if (IIN <= ILEN) GO TO 15
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
if (ILEN <= 0) GO TO 6
IIN = 1
15 if (ICHK /= 0) GO TO 460
if (ICOP1(IIN) == 0) GO TO 10
IOUT = IOUT+1
GO TO 110
10 ICHK = 1
GO TO 110
460 ICHK = 0
NSAVE = ICOP1(IIN)
6 continue
NIB = (NI-1)*NVIRT+NOVST
do NB=1,NVIRT
  NSIB = MUL(NSM(LN+NB),NSM(NI))
  INS = NNS(NSIB)
  if (INS == 0) GO TO 18
  do I=1,INS
    BIAC(I) = Zero
    BICA(I) = Zero
  end do
18 NIB = NIB+1
  IADR = LASTAD(NIB)
203 call iDAFILE(Lu_TiABIJ,2,INDOUT,KBUF2,IADR)
  LENGTH = INDOUT(KBUF1)
  IADR = INDOUT(KBUF2)
  if (LENGTH == 0) GO TO 210
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
210 if (IADR /= -1) GO TO 203
  ILOOP = 0
72 do I=1,INS
    INSOUT = INSOUT+1
    if (ILOOP == 0) BUFBI(INSOUT) = BIAC(I)
    if (ILOOP == 1) BUFBI(INSOUT) = BICA(I)
    if (INSOUT < KBUFF1) GO TO 75
    call dDAFILE(Lu_TiABCI,1,BUFBI,KBUFF1,IAD15)
    INSOUT = 0
75  continue
  end do
  ILOOP = ILOOP+1
  if (ILOOP == 1) GO TO 72
end do
if (ILEN >= 0) GO TO 100
! EMPTY LAST BUFFER
if (INSOUT == 0) return
call dDAFILE(Lu_TiABCI,1,BUFBI,KBUFF1,IAD15)

return

end subroutine SORTA_CPF
