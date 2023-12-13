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

subroutine SORTB_CPF(BUFOUT,INDOUT,ICAD,IBUFL,TIBUF,ACBDS,ACBDT,ISAB,BUFACBD)
! SORTS INTEGRALS (AB/CD) FOR FIXED A,C ALL B,D

use cpf_global, only: ICH, IPASS, IRC, IROW, JBUF, JJS, KBUFF1, LASTAD, LN, LSYM, Lu_TiABCD, Lu_TiABIJ, Lu_TraInt, MADR, NORB, &
                      NSM, NSYM, NSYS, NTIBUF, NVIRT
use Symmetry_Info, only: Mul
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6, RtoI

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: BUFOUT(*), TIBUF(NTIBUF), ACBDS(*), ACBDT(*), BUFACBD(*)
integer(kind=iwp), intent(_OUT_) :: INDOUT(*), ICAD(*), IBUFL(*)
integer(kind=iwp), intent(in) :: ISAB(*)
#include "tratoc.fh"
integer(kind=iwp) :: I, IAC, IACMAX, IACMIN, IAD16, IAD50, IADR, IBDS, ICP, ICPP, ICQ, ID, IDISK, IDIV, IFIN1, IFIN2, ILOOP, IN1, &
                     INND, INPS, INPT, INS, INSOUT, IOUT, IREC, IST1, IST2, ISTEP, ISYM, ITAIL, ITURN, JBUF0, JBUF1, JBUF2, JDISK, &
                     KK, LENGTH, M1, M2, M3, M4, N1, N2, N3, N4, NA, NAC, NB, NBD, NC, ND, NDMAX, NI, NJ, NK, NL, NOP, NOQ, NOR, &
                     NORB0(9), NORBP, NOS, NOV, NOVM, NOVST, NSAC, NSACL, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NSSM, NT, NTM, NU, &
                     NUMAX, NUMIN, NV, NVT, NX, NXM
real(kind=wp) :: FINI

NVT = IROW(NVIRT+1)
NOV = (NVT-1)/IPASS+1
NOVST = LN*NVIRT+1
IAD16 = 0
JBUF0 = RTOI*JBUF
JBUF1 = JBUF0+JBUF+1
JBUF2 = JBUF1+1
IDIV = RTOI
NORB0(1) = 0
do I=1,NSYM
  NORB0(I+1) = NORB0(I)+NORB(I)
end do
INSOUT = 0
IACMAX = 0
do ISTEP=1,IPASS
  IAD50 = 0
  call iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
  IDISK = 0
  IACMIN = IACMAX+1
  IACMAX = IACMAX+NOV
  if (IACMAX > NVT) IACMAX = NVT
  if (IACMIN > IACMAX) cycle
  ID = 0
  do IREC=1,NOV
    IBUFL(IREC) = 0
    ICAD(IREC) = ID
    INDOUT(ID+JBUF2) = -1
    ID = ID+JBUF2
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
                  if ((M1 <= LN) .or. (M2 <= LN) .or. (M3 <= LN) .or. (M4 <= LN)) cycle
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
                  if (abs(FINI) < 1.0e-9_wp) cycle
                  NA = NI-LN
                  NB = NJ-LN
                  NC = NK-LN
                  ND = NL-LN
                  ITURN = 0
                  if ((NA == NB) .and. (NC == ND)) cycle
                  do
                    IAC = IROW(NA)+NC
                    if ((IAC >= IACMIN) .and. (IAC <= IACMAX)) then
                      if ((NA == NC) .and. (NB == ND)) FINI = FINI*Half
                      NAC = IAC-IACMIN+1
                      IBUFL(NAC) = IBUFL(NAC)+1
                      ICQ = ICAD(NAC)
                      ICP = ICQ/IDIV+IBUFL(NAC)
                      BUFOUT(ICP) = FINI
                      ICPP = ICQ+JBUF0+IBUFL(NAC)
                      !PAM97 INDOUT(ICPP) = ior(NB,ishft(ND,8))
                      INDOUT(ICPP) = NB+ND*2**8
                      if (IBUFL(NAC) >= JBUF) then
                        INDOUT(ICQ+JBUF1) = JBUF
                        JDISK = IDISK
                        call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),JBUF2,IDISK)
                        INDOUT(ICQ+JBUF2) = JDISK
                        IBUFL(NAC) = 0
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
  if ((NOVST+IACMIN-1+NOVM) > MADR) then
    write(u6,*) 'SORTB_CPF Error: NOVST+IACMIN-1+NOVM > MADR'
    write(u6,*) '  (See code).'
    call Abend()
  end if
  do I=1,NOVM
    ICQ = ICAD(I)
    INDOUT(ICQ+JBUF1) = IBUFL(I)
    JDISK = IDISK
    call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),JBUF2,IDISK)
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
        NDMAX = NSYS(NSM(LN+NC)+1)
        if (NDMAX > NA) NDMAX = NA
        INS = ISAB(IN1+NDMAX)
        ACBDS(1:INS) = Zero
        ACBDT(1:INS) = Zero
        IADR = LASTAD(NOVST+IAC)
        do
          call iDAFILE(Lu_TiABIJ,2,INDOUT,JBUF2,IADR)
          LENGTH = INDOUT(JBUF1)
          IADR = INDOUT(JBUF2)
          do KK=1,LENGTH
            INND = INDOUT(JBUF0+KK)
            NB = ibits(INND,0,8)
            ND = ibits(INND,8,8)
            NBD = (NB-1)*NVIRT+ND
            IBDS = ISAB(NBD)
            ACBDS(IBDS) = ACBDS(IBDS)+BUFOUT(KK)
            if (NB > ND) ACBDT(IBDS) = ACBDT(IBDS)+BUFOUT(KK)
            if (NB < ND) ACBDT(IBDS) = ACBDT(IBDS)-BUFOUT(KK)
          end do
          if (IADR == -1) exit
        end do
        ILOOP = 0
        do
          do I=1,INS
            INSOUT = INSOUT+1
            if (ILOOP == 0) BUFACBD(INSOUT) = ACBDS(I)
            if (ILOOP == 1) BUFACBD(INSOUT) = ACBDT(I)
            if (INSOUT >= KBUFF1) then
              call dDAFILE(Lu_TiABCD,1,BUFACBD,KBUFF1,IAD16)
              INSOUT = 0
            end if
          end do
          ILOOP = ILOOP+1
          if (ILOOP /= 1) exit
        end do
      end do
    end do
  end do
end do
! EMPTY LAST BUFFER
if (INSOUT /= 0) call dDAFILE(Lu_TiABCD,1,BUFACBD,KBUFF1,IAD16)

return

end subroutine SORTB_CPF
