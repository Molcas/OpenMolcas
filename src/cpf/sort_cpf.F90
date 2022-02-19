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

subroutine SORT_CPF(BUFOUT,INDOUT,ICAD,IBUFL,FC,FIJ,FJI,TIBUF)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension BUFOUT(*), INDOUT(*)
dimension ICAD(*), FC(*), IBUFL(*), FIJ(*), FJI(*), TIBUF(*)
dimension IVEC(20), IPOF(65)
dimension NORB0(9)

IAD50 = 0
call iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
NVT = IROW(NVIRT+1)
do I=1,20
  IVEC(I) = 0
end do
IN = 1
do I=1,NSYM
  call IPO_CPF(IPOF(IN),NVIR,MUL,NSYM,I,-1)
  IN = IN+NSYM
end do
! ORDER OF RECORD-CHAINS IS
! 1. NOT2 CHAINS (AB/IJ)
! 2. NOT2 CHAINS (AI/BJ)
! 3. NOT2 CHAINS (AI/JK)
! RECORD STRUCTURE IS
! 1. LBUF INTEGRALS
! 2. LBUF INDICES
! 3. NUMBER OF INTEGRALS IN THIS RECORD
! 4. ADDRESS OF LAST RECORD
NOT2 = IROW(LN+1)
NOV = 3*NOT2
NOTT = 2*NOT2
NOVST = LN*NVIRT+1+NVT
IDISK = 0
LBUF0 = RTOI*LBUF
LBUF1 = LBUF0+LBUF+1
LBUF2 = LBUF1+1
IDIV = RTOI
ID = 0
do IREC=1,NOV
  IBUFL(IREC) = 0
  ICAD(IREC) = ID
  INDOUT(ID+LBUF2) = -1
  ID = ID+LBUF2
end do
NORB0(1) = 0
do I=1,NSYM
  NORB0(I+1) = NORB0(I)+NORB(I)
end do

! ONE-ELECTRON INTEGRALS

NORBTT = 0
do ISYM=1,nsym
  NORBTT = NORBTT+(NORB(ISYM)*(NORB(ISYM)+1))/2
end do
EMY = POTNUC
NOB2 = IROW(NORBT+1)
IADD17 = ITOC17(2)
call dDAFILE(Lu_TraOne,2,FIJ,NORBTT,IADD17)
call DCOPY_(NOB2,[0.0d0],0,FC,1)
IBUF = 0
KORBI = 0
do ISYM=1,NSYM
  do JORBI=KORBI+1,KORBI+NORB(ISYM)
    do IORBI=KORBI+1,JORBI
      IBUF = IBUF+1
      ONEHAM = FIJ(IBUF)
      NI = ICH(IORBI)
      NJ = ICH(JORBI)
      if ((NI == 0) .or. (NJ == 0)) GO TO 199
      if (NI < NJ) then
        NTMP = NI
        NI = NJ
        NJ = NTMP
      end if
      if (NJ > 0) then
        IJT = IROW(NI)+NJ
        FC(IJT) = FC(IJT)+ONEHAM
      else if (NI == NJ) then
        EMY = EMY+2.0d0*ONEHAM
      end if
199   continue
    end do
  end do
  KORBI = KORBI+NORB(ISYM)
end do
call DCOPY_(NOB2,[0.0d0],0,FIJ,1)
call DCOPY_(NOB2,[0.0d0],0,FJI,1)
if (IPRINT >= 20) then
  call TRIPRT('FC IN SORT BEFORE TWOEL',' ',FC,NORBT)
  write(6,'(A,F20.8)') ' EMY:',EMY
  call XFLUSH(6)
end if

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
                if ((M1 == 0) .or. (M2 == 0)) GO TO 306
                if ((M3 == 0) .or. (M4 == 0)) GO TO 306
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
                if (abs(FINI) < 1.D-09) GO TO 306
                if ((NI <= 0) .or. (NJ <= 0)) GO TO 41
                if ((NK <= 0) .or. (NL <= 0)) GO TO 41
                DFINI = abs(FINI)
                IEXP = int(-log10(DFINI+1.0D-20)+5)
                if (IEXP > 20) IEXP = 20
                if (IEXP < 1) IEXP = 1
                IVEC(IEXP) = IVEC(IEXP)+1
                if ((NI /= NJ) .or. (NK /= NL)) GO TO 42
                IJ = IROW(NI)+NK
                FIJ(IJ) = FINI
                ! SKIP (AA/II) INTEGRALS
                GO TO 306
42              if ((NI /= NK) .or. (NJ /= NL)) GO TO 43
                IJ = IROW(NI)+NJ
                FJI(IJ) = FINI
43              if (NI <= LN) GO TO 306
                if (NJ > LN) GO TO 102
                if (NK > LN) GO TO 103
                ! AIJK
                JK = NOTT+IROW(NK)+NL
                IBUFL(JK) = IBUFL(JK)+1
                ICQ = ICAD(JK)
                ICP = ICQ/IDIV+IBUFL(JK)
                BUFOUT(ICP) = FINI
                ICPP = ICQ+LBUF0+IBUFL(JK)
                INDOUT(ICPP) = IROW(NI)+NJ
                if (IBUFL(JK) < LBUF) GO TO 306
                INDOUT(ICQ+LBUF1) = LBUF
                JDISK = IDISK
                call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
                INDOUT(ICQ+LBUF2) = JDISK
                IBUFL(JK) = 0
                GO TO 306
103             if (NL > LN) GO TO 306
                ! AIBJ
                IIJ = NOT2+IROW(NJ)+NL
                if (NL > NJ) IIJ = NOT2+IROW(NL)+NJ
                IBUFL(IIJ) = IBUFL(IIJ)+1
                ICQ = ICAD(IIJ)
                ICP = ICQ/IDIV+IBUFL(IIJ)
                BUFOUT(ICP) = FINI
                NSA = NSM(NI)
                NAV = NI-LN-NSYS(NSA)
                NSB = NSM(NK)
                NBV = NK-LN-NSYS(NSB)
                NSIJT = (MUL(NSM(NJ),NSM(NL))-1)*NSYM
                if (NL > NJ) GO TO 105
                INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                GO TO 104
105             INAV = IPOF(NSIJT+NSB)+(NAV-1)*NVIR(NSB)+NBV
104             ICPP = ICQ+LBUF0+IBUFL(IIJ)
                INDOUT(ICPP) = INAV
                if (IBUFL(IIJ) < LBUF) GO TO 108
                INDOUT(ICQ+LBUF1) = LBUF
                JDISK = IDISK
                call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
                INDOUT(ICQ+LBUF2) = JDISK
                IBUFL(IIJ) = 0
108             if (NJ /= NL) GO TO 306
                if (NI == NK) GO TO 306
                JNAV = IROW(NI)+NK
                FC(JNAV) = FC(JNAV)-FINI
                GO TO 306
102             if (NK > LN) GO TO 306
                ! ABIJ
                IIJ = IROW(NK)+NL
                IBUFL(IIJ) = IBUFL(IIJ)+1
                ICQ = ICAD(IIJ)
                ICP = ICQ/IDIV+IBUFL(IIJ)
                BUFOUT(ICP) = FINI
                NSA = NSM(NI)
                NAV = NI-LN-NSYS(NSA)
                NSB = NSM(NJ)
                NBV = NJ-LN-NSYS(NSB)
                NSIJT = (MUL(NSM(NK),NSM(NL))-1)*NSYM
                INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                ICPP = ICQ+LBUF0+IBUFL(IIJ)
                INDOUT(ICPP) = INAV
                if (IBUFL(IIJ) < LBUF) GO TO 106
                INDOUT(ICQ+LBUF1) = LBUF
                JDISK = IDISK
                call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
                INDOUT(ICQ+LBUF2) = JDISK
                IBUFL(IIJ) = 0
106             if (NK /= NL) GO TO 306
                if (NI == NJ) GO TO 306
                JNAV = IROW(NI)+NJ
                FC(JNAV) = FC(JNAV)+D2*FINI
                GO TO 306
                ! CHECK FOR FOCK-MATRIX CONTRIBUTION
41              if (NI /= NJ) GO TO 51
                II = 1
                call IFOCK(FC,NI,NK,NL,FINI,II)
                if (NK /= NL) GO TO 52
                if ((NI > 0) .or. (NK > 0)) GO TO 57
                EMY = EMY+D2*FINI
                if (NI /= NK) EMY = EMY+D2*FINI
57              if (NI == NK) GO TO 52
                call IFOCK(FC,NK,NI,NJ,FINI,II)
                GO TO 52
51              if (NK /= NL) GO TO 52
                II = 1
                call IFOCK(FC,NK,NI,NJ,FINI,II)
52              II = 0
                if (NI /= NK) GO TO 53
                call IFOCK(FC,NI,NJ,NL,FINI,II)
                if (NJ /= NL) GO TO 306
                if ((NI > 0) .or. (NJ > 0)) GO TO 58
                EMY = EMY-FINI
                if (NI /= NJ) EMY = EMY-FINI
58              if (NI == NJ) GO TO 306
                call IFOCK(FC,NJ,NI,NK,FINI,II)
                GO TO 306
53              if (NI /= NL) GO TO 54
                call IFOCK(FC,NI,NJ,NK,FINI,II)
                GO TO 306
54              if (NJ /= NK) GO TO 55
                call IFOCK(FC,NJ,NI,NL,FINI,II)
                GO TO 306
55              if (NJ /= NL) GO TO 306
                call IFOCK(FC,NJ,NI,NK,FINI,II)
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
if ((NOVST+NOV) > mAdr) then
  write(6,*) 'SORT Error: NOVST+NOV>MADR (See code).'
  call Abend()
end if
do I=1,NOV
  ICQ = ICAD(I)
  INDOUT(ICQ+LBUF1) = IBUFL(I)
  JDISK = IDISK
  call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
  LASTAD(NOVST+I) = JDISK
end do
do J=1,NORBT
  IND = IROW(J+1)
  FC(IND) = FC(IND)+EMY/N
end do
IADD25 = 0
call dDAFILE(Lu_25,1,FC,NOB2,IADD25)
IAD25S = IADD25
write(6,154)
call XFLUSH(6)
154 format(//6X,'STATISTICS FOR INTEGRALS, FIRST ENTRY 10**3-10**4',/)
write(6,155) (IVEC(I),I=1,20)
call XFLUSH(6)
155 format(6X,5I10)

return

end subroutine SORT_CPF
