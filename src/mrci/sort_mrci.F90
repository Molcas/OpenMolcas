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

!PAM04 subroutine SORT(BUFOUT,INDOUT,FC,FIIJJ,FIJIJ,NINTGR)
subroutine SORT_MRCI(BUFS,INDS,FC,FIIJJ,FIJIJ,NINTGR)

use Definitions, only: wp, iwp, u6

implicit none
#include "mrci.fh"
real(kind=wp) :: BUFS(NBITM3,NCHN3), FC(NBTRI), FIIJJ(*), FIJIJ(*)
integer(kind=iwp) :: INDS(NBITM3+2,NCHN3), NINTGR
#include "warnings.h"
integer(kind=iwp) :: I, IAD50, IADD17, IADD25, IBUF, IDISK, IEXP, IIJ, IIN, IJ, IJT, IKT, INAV, IND, IORBI, IOUT, IPOF(65), IPOS, &
                     IREC, ISYM, IVEC(20), J, JDISK, JK, JORBI, KORBI, M1, M2, M3, M4, N1, N2, N3, N4, NAV, NBV, NI, NJ, NK, NL, &
                     NOP, NOQ, NOR, NORB0(9), NORBP, NORBTT, NOS, NOT2, NOTT, NOVST, NSA, NSB, NSIJT, NSP, NSPQ, NSPQR, NSQ, NSR, &
                     NSS, NSSM, NT, NTM, NTMP, NU, NUMAX, NUMIN, NV, NVT, NX, NXM
real(kind=wp) :: DFINI, EFROZ, FINI, ONEHAM

IAD50 = 0
call iDAFILE(LUTRA,2,iTraToc,nTraToc,IAD50)
NVT = IROW(NVIRT+1)
do I=1,20
  IVEC(I) = 0
end do
IIN = 1
do I=1,NSYM
  call IPO(IPOF(IIN),NVIR,MUL,NSYM,I,-1)
  IIN = IIN+NSYM
end do
! ORDER OF RECORD-CHAINS IS
! 1.  NOT2 CHAINS (AB/IJ)
! 2.  NOT2 CHAINS (AI/BJ)
! 3.  NOT2 CHAINS (AI/JK)
! RECORD STRUCTURE IS
! 1.  NBITM3 INTEGRALS
! 2.  NBITM3 INDICES
! 3.  NUMBER OF INTEGRALS IN THIS RECORD
! 4.  ADDRESS OF LAST RECORD
NOT2 = IROW(LN+1)
NOTT = 2*NOT2
NOVST = LN*NVIRT+1+NVT
IDISK = 0
!PAM97 The portable code should then be:
!PAM04 NBITM3 = (RTOI*NBSIZ3-2)/(RTOI+1)
!PAM04 IBOFF3 = RTOI*NBITM3
!PAM04 IBBC3 = IBOFF3+NBITM3+1
!PAM04 IBDA3 = IBBC3+1

do IREC=1,NCHN3
  !PAM04 INDOUT(IBBC3+(IREC-1)*RTOI*NBSIZ3) =  0
  !PAM04 INDOUT(IBDA3+(IREC-1)*RTOI*NBSIZ3) = -1
  INDS(NBITM3+1,IREC) = 0
  INDS(NBITM3+2,IREC) = -1
end do
NORB0(1) = 0
do I=1,NSYM
  NORB0(I+1) = NORB0(I)+NORB(I)
end do
! READ ONE-ELECTRON ORBITALS. USE FIIJJ TEMPORARILY AS READ BUFFER.
NORBTT = 0
do ISYM=1,nsym
  NORBTT = NORBTT+(NORB(ISYM)*(NORB(ISYM)+1))/2
end do
EFROZ = POTNUC
call FZERO(FC,NBTRI)
IADD17 = ITOC17(2)
call dDAFILE(LUONE,2,FIIJJ,NORBTT,IADD17)
IBUF = 0
KORBI = 0
do ISYM=1,NSYM
  do JORBI=KORBI+1,KORBI+NORB(ISYM)
    do IORBI=KORBI+1,JORBI
      IBUF = IBUF+1
      ONEHAM = FIIJJ(IBUF)
      NI = ICH(IORBI)
      NJ = ICH(JORBI)
      if ((NI == 0) .or. (NJ == 0)) GO TO 198
      if (NI < NJ) then
        NTMP = NI
        NI = NJ
        NJ = NTMP
      end if
      if (NJ > 0) then
        IJT = IROW(NI)+NJ
        FC(IJT) = FC(IJT)+ONEHAM
      else if (NI == NJ) then
        EFROZ = EFROZ+2*ONEHAM
      end if
198   continue
    end do
  end do
  KORBI = KORBI+NORB(ISYM)
end do
if (IPRINT >= 20) then
  call TRIPRT('FC IN SORT_MRCI BEFORE TWOEL',' ',FC,NORBT)
  write(u6,'(A,F20.8)') ' EFROZ:',EFROZ
  call XFLUSH(u6)
end if
call FZERO(FIIJJ,NBTRI)
call FZERO(FIJIJ,NBTRI)
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
        call dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)
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
                  call dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)
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
                if ((NI <= 0) .or. (NJ <= 0)) GO TO 41
                if ((NK <= 0) .or. (NL <= 0)) GO TO 41
                DFINI = abs(FINI)+1.0e-20_wp
                IEXP = int(-log10(DFINI))+5
                if (IEXP <= 20) IVEC(IEXP) = IVEC(IEXP)+1
                if ((NI /= NJ) .or. (NK /= NL)) GO TO 42
                IJ = IROW(NI)+NK
                FIIJJ(IJ) = FINI
                ! SKIP (AA/II) INTEGRALS
                GO TO 306
42              if ((NI /= NK) .or. (NJ /= NL)) GO TO 43
                IJ = IROW(NI)+NJ
                FIJIJ(IJ) = FINI
43              if (NI <= LN) GO TO 306
                if (NJ > LN) GO TO 102
                if (NK > LN) GO TO 103
                ! AIJK
                JK = NOTT+IROW(NK)+NL
                !PAM04 IPOS = INDOUT(IBBC3+(JK-1)*RTOI*NBSIZ3)+1
                !PAM04 INDOUT(IBBC3+(JK-1)*RTOI*NBSIZ3) = IPOS
                IPOS = INDS(NBITM3+1,JK)+1
                INDS(NBITM3+1,JK) = IPOS
                !BUFOUT(IPOS,JK) = FINI
                !PAM04 BUFOUT(IPOS+(JK-1)*NBSIZ3) = FINI
                BUFS(IPOS,JK) = FINI
                !PAM04 INDOUT(IBOFF3+IPOS+(JK-1)*RTOI*NBSIZ3) = IROW(NI)+NJ
                INDS(IPOS,JK) = IROW(NI)+NJ
                if (IPOS < NBITM3) GO TO 306
                JDISK = IDISK
                !PAM04 call dDAFILE(Lu_60,1,INDOUT(1+(JK-1)*RTOI*NBSIZ3),NBSIZ3,IDISK)
                call iDAFILE(Lu_60,1,INDS(1,JK),NBITM3+2,IDISK)
                call dDAFILE(Lu_60,1,BUFS(1,JK),NBITM3,IDISK)
                !PAM04 INDOUT(IBBC3+(JK-1)*RTOI*NBSIZ3) = 0
                !PAM04 INDOUT(IBDA3+(JK-1)*RTOI*NBSIZ3) = JDISK
                INDS(NBITM3+1,JK) = 0
                INDS(NBITM3+2,JK) = JDISK
                GO TO 306
103             if (NL > LN) GO TO 306
                ! AIBJ
                IIJ = NOT2+IROW(NJ)+NL
                if (NL > NJ) IIJ = NOT2+IROW(NL)+NJ
                !PAM04 IPOS = INDOUT(IBBC3+(IIJ-1)*RTOI*NBSIZ3)+1
                !PAM04 INDOUT(IBBC3+(IIJ-1)*RTOI*NBSIZ3) = IPOS
                IPOS = INDS(NBITM3+1,IIJ)+1
                INDS(NBITM3+1,IIJ) = IPOS
                !BUFOUT(IPOS,IIJ) = FINI
                !PAM04 BUFOUT(IPOS+(IIJ-1)*NBSIZ3) = FINI
                BUFS(IPOS,IIJ) = FINI
                NSA = NSM(NI)
                NAV = NI-LN-NVIRP(NSA)
                NSB = NSM(NK)
                NBV = NK-LN-NVIRP(NSB)
                NSIJT = (MUL(NSM(NJ),NSM(NL))-1)*NSYM
                if (NL > NJ) GO TO 105
                INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                GO TO 104
105             INAV = IPOF(NSIJT+NSB)+(NAV-1)*NVIR(NSB)+NBV
104             continue
                !PAM04 INDOUT(IPOS+IBOFF3+(IIJ-1)*RTOI*NBSIZ3) = INAV
                INDS(IPOS,IIJ) = INAV
                if (IPOS < NBITM3) GO TO 108
                JDISK = IDISK
                !PAM04 call dDAFILE(Lu_60,1,INDOUT(1+(IIJ-1)*RTOI*NBSIZ3),NBSIZ3,IDISK)
                call iDAFILE(Lu_60,1,INDS(1,IIJ),NBITM3+2,IDISK)
                call dDAFILE(Lu_60,1,BUFS(1,IIJ),NBITM3,IDISK)
                !PAM04 INDOUT(IBBC3+(IIJ-1)*RTOI*NBSIZ3) = 0
                !PAM04 INDOUT(IBDA3+(IIJ-1)*RTOI*NBSIZ3) = JDISK
                INDS(NBITM3+1,IIJ) = 0
                INDS(NBITM3+2,IIJ) = JDISK
108             if (NJ /= NL) GO TO 306
                if (NI == NK) GO TO 306
                IKT = IROW(NI)+NK
                FC(IKT) = FC(IKT)-FINI
                GO TO 306
102             if (NK > LN) GO TO 306
                ! ABIJ
                IIJ = IROW(NK)+NL
                !PAM04 IPOS = INDOUT(IBBC3+(IIJ-1)*RTOI*NBSIZ3)+1
                !PAM04 INDOUT(IBBC3+(IIJ-1)*RTOI*NBSIZ3) = IPOS
                IPOS = INDS(NBITM3+1,IIJ)+1
                INDS(NBITM3+1,IIJ) = IPOS
                !PAM04 BUFOUT(IPOS+(IIJ-1)*NBSIZ3) = FINI
                BUFS(IPOS,IIJ) = FINI
                NSA = NSM(NI)
                NAV = NI-LN-NVIRP(NSA)
                NSB = NSM(NJ)
                NBV = NJ-LN-NVIRP(NSB)
                NSIJT = (MUL(NSM(NK),NSM(NL))-1)*NSYM
                INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                !PAM04 INDOUT(IBOFF3+IPOS+(IIJ-1)*RTOI*NBSIZ3) = INAV
                INDS(IPOS,IIJ) = INAV
                if (IPOS < NBITM3) GO TO 106
                JDISK = IDISK
                !PAM04 call dDAFILE(Lu_60,1,INDOUT(1+(IIJ-1)*RTOI*NBSIZ3),NBSIZ3,IDISK)
                call iDAFILE(Lu_60,1,INDS(1,IIJ),NBITM3+2,IDISK)
                call dDAFILE(Lu_60,1,BUFS(1,IIJ),NBITM3,IDISK)
                !PAM04 INDOUT(IBBC3+(IIJ-1)*RTOI*NBSIZ3) = 0
                !PAM04 INDOUT(IBDA3+(IIJ-1)*RTOI*NBSIZ3) = JDISK
                INDS(NBITM3+1,IIJ) = 0
                INDS(NBITM3+2,IIJ) = JDISK
106             if (NK /= NL) GO TO 306
                if (NI == NJ) GO TO 306
                IJT = IROW(NI)+NJ
                FC(IJT) = FC(IJT)+2*FINI
                GO TO 306
                ! CHECK FOR FOCK-MATRIX, AND FROZEN ENERGY, CONTRIBUTIONS
41              continue
                if (NI < 0) then
                  if ((NI == NJ) .and. (NK == NL)) EFROZ = EFROZ+4*FINI
                  if ((NI == NK) .and. (NJ == NL)) EFROZ = EFROZ-2*FINI
                else if (NL < 0) then
                  if ((NK == NL) .and. (NJ > 0)) then
                    IJT = IROW(NI)+NJ
                    FC(IJT) = FC(IJT)+2*FINI
                  else if ((NJ == NL) .and. (NK > 0)) then
                    IKT = IROW(NI)+NK
                    FC(IKT) = FC(IKT)-FINI
                  end if
                end if
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
if ((NOVST+NCHN3) > mChain) then
  write(u6,*) 'SORT_MRCI Error: NOVST+NCHN3>MCHAIN (See code).'
  call QUIT(_RC_GENERAL_ERROR_)
end if
do I=1,NCHN3
  JDISK = IDISK
  !PAM04 call dDAFILE(Lu_60,1,INDOUT(1+(I-1)*RTOI*NBSIZ3),NBSIZ3,IDISK)
  call iDAFILE(Lu_60,1,INDS(1,I),NBITM3+2,IDISK)
  call dDAFILE(Lu_60,1,BUFS(1,I),NBITM3,IDISK)
  LASTAD(NOVST+I) = JDISK
end do
do J=1,NORBT
  IND = IROW(J+1)
  FC(IND) = FC(IND)+EFROZ/NELEC
end do
IADD25 = 0
call dDAFILE(Lu_25,1,FC,NBTRI,IADD25)
IAD25S = IADD25
!if (IPRINT >= 2) then
write(u6,154)
call XFLUSH(u6)
write(u6,155) (IVEC(I),I=1,20)
154 format(//6X,'STATISTICS FOR INTEGRALS, FIRST ENTRY 10**3-10**4',/)
call XFLUSH(u6)
155 format(6X,5I10)
!end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(NINTGR)

end subroutine SORT_MRCI
