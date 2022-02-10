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

subroutine SORTA(BUFS,INDS,ISAB,BUFBI,BIAC,BICA,NINTGR)
! SORTS INTEGRALS (AB/CI)
! FOR FIXED B,I ALL A,C
! FIRST CHAIN FOR IJKL

implicit real*8(A-H,O-Z)
external COUNT
#include "SysDef.fh"
#include "warnings.h"
#include "mrci.fh"
dimension BUFS(NBITM1,NCHN1)
dimension INDS(NBITM1+2,NCHN1)
dimension BUFBI(KBUFF1)
dimension BIAC(ISMAX), BICA(ISMAX)
dimension ISAB(*)
dimension NORB0(9)

call count(NINTGR,NSYM,NORB,MUL)
if (IPRINT >= 6) write(6,1234) NINTGR
call XFLUSH(6)
1234 format(//6X,'NUMBER OF TWO-ELECTRON INTEGRALS',I10)

IAD50 = 0
call iDAFILE(LUTRA,2,iTraToc,nTraToc,IAD50)
!IBOFF1 = RTOI*NBITM1
!IBBC1 = IBOFF1+NBITM1+1
!IBDA1 = IBBC1+1

IDISK = 0
ICHK = 0
do IREC=1,NCHN1
  !INDOUT(IBBC1+(IREC-1)*RTOI*NBSIZ1) = 0
  !INDOUT(IBDA1+(IREC-1)*RTOI*NBSIZ1) = -1
  INDS(NBITM1+1,IREC) = 0
  INDS(NBITM1+2,IREC) = -1
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
                FINI = TIBUF(IOUT)
                NI = ICH(NORB0(NSP)+NT)
                if (NI <= 0) goto 306
                NJ = ICH(NORB0(NSQ)+NU)
                if (NJ <= 0) goto 306
                NK = ICH(NORB0(NSR)+NV)
                if (NK <= 0) goto 306
                NL = ICH(NORB0(NSS)+NX)
                if (NL <= 0) goto 306
                ! ORDER THESE INDICES CANONICALLY
                if (NI < NJ) then
                  M = NI
                  NI = NJ
                  NJ = M
                end if
                if (NK < NL) then
                  M = NK
                  NK = NL
                  NL = M
                end if
                if (NI < NK) then
                  M = NK
                  NK = NI
                  NI = M
                  M = NL
                  NL = NJ
                  NJ = M
                else if ((NI == NK) .and. (NJ < NL)) then
                  M = NL
                  NL = NJ
                  NJ = M
                end if
                if (NI <= LN) GO TO 109
                if (NK <= LN) GO TO 306
                if (IFIRST /= 0) GO TO 306
                if (NJ <= LN) GO TO 42
                if (NL > LN) GO TO 306
                ! ABCI
                NA = NI-LN
                NB = NJ-LN
                NC = NK-LN
                NI = NL
                GO TO 108
42              continue
                if (NL <= LN) goto 306
                ! CIAB
                NA = NK-LN
                NB = NL-LN
                NC = NI-LN
                NI = NJ
108             NIB = (NI-1)*NVIRT+NB+1
                !IPOS = INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)+1
                !INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1) = IPOS
                IPOS = INDS(NBITM1+1,NIB)+1
                INDS(NBITM1+1,NIB) = IPOS
                !BUFOUT(IPOS+(NIB-1)*NBSIZ1) = FINI
                !INDOUT(IBOFF1+IPOS+(NIB-1)*RTOI*NBSIZ1) = (NA-1)*NVIRT+NC
                BUFS(IPOS,NIB) = FINI
                INDS(IPOS,NIB) = (NA-1)*NVIRT+NC
                if (IPOS >= NBITM1) then
                  JDISK = IDISK
                  !call dDAFILE(Lu_60,1,INDOUT(1+(NIB-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
                  call iDAFILE(Lu_60,1,INDS(1,NIB),NBITM1+2,IDISK)
                  call dDAFILE(Lu_60,1,BUFS(1,NIB),NBITM1,IDISK)
                  !INDOUT(IBDA1+(NIB-1)*RTOI*NBSIZ1) = JDISK
                  !INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1) = 0
                  INDS(NBITM1+1,NIB) = 0
                  INDS(NBITM1+2,NIB) = JDISK
                end if
                if (NA == NB) GO TO 306
                NAT = NA
                NA = NB
                NB = NAT
                NIB = (NI-1)*NVIRT+NB+1
                !IPOS = INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)+1
                !INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1) = IPOS
                IPOS = INDS(NBITM1+1,NIB)+1
                INDS(NBITM1+1,NIB) = IPOS
                !BUFOUT(IPOS+(NIB-1)*NBSIZ1) = FINI
                !INDOUT(IBOFF1+IPOS+(NIB-1)*RTOI*NBSIZ1) = (NA-1)*NVIRT+NC
                BUFS(IPOS,NIB) = FINI
                INDS(IPOS,NIB) = (NA-1)*NVIRT+NC
                if (IPOS >= NBITM1) then
                  JDISK = IDISK
                  !call dDAFILE(Lu_60,1,INDOUT(1+(NIB-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
                  call iDAFILE(Lu_60,1,INDS(1,NIB),NBITM1+2,IDISK)
                  call dDAFILE(Lu_60,1,BUFS(1,NIB),NBITM1,IDISK)
                  !INDOUT(IBDA1+(NIB-1)*RTOI*NBSIZ1) = JDISK
                  !INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1) = 0
                  INDS(NBITM1+1,NIB) = 0
                  INDS(NBITM1+2,NIB) = JDISK
                end if
                goto 306
                ! IJKL
109             IIJ = IROW(NI)+NJ
                KL = IROW(NK)+NL
                IJKL = IIJ*(IIJ-1)/2+KL
                IJ = 1
                !IPOS = INDOUT(IBBC1+(IJ-1)*RTOI*NBSIZ1)+1
                !INDOUT(IBBC1+(IJ-1)*RTOI*NBSIZ1) = IPOS
                IPOS = INDS(NBITM1+1,IJ)+1
                INDS(NBITM1+1,IJ) = IPOS
                !BUFOUT(IPOS+(IJ-1)*NBSIZ1) = FINI
                !INDOUT(IBOFF1+IPOS+(IJ-1)*RTOI*NBSIZ1) = IJKL
                BUFS(IPOS,IJ) = FINI
                INDS(IPOS,IJ) = IJKL
                if (IPOS == NBITM1) then
                  JDISK = IDISK
                  !call dDAFILE(Lu_60,1,INDOUT(1+(IJ-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
                  call iDAFILE(Lu_60,1,INDS(1,IJ),NBITM1+2,IDISK)
                  call dDAFILE(Lu_60,1,BUFS(1,IJ),NBITM1,IDISK)
                  !INDOUT(IBDA1+(IJ-1)*RTOI*NBSIZ1) = JDISK
                  !INDOUT(IBBC1+(IJ-1)*RTOI*NBSIZ1) = 0
                  INDS(NBITM1+1,IJ) = 0
                  INDS(NBITM1+2,IJ) = JDISK
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
if (NChn1 > mChain) then
  write(6,*) 'SORTA Error: NCHN1 > MCHAIN (See code).'
  call QUIT(_RC_GENERAL_ERROR_)
end if
do I=1,NCHN1
  JDISK = IDISK
  !call dDAFILE(Lu_60,1,INDOUT(1+(I-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
  call iDAFILE(Lu_60,1,INDS(1,I),NBITM1+2,IDISK)
  call dDAFILE(Lu_60,1,BUFS(1,I),NBITM1,IDISK)
  LASTAD(I) = JDISK
end do
! IJKL
IDISK = 0
IBUFIJ = 0
ISRTAD = -1
IADR = LASTAD(1)
201 continue
!call dDAFILE(Lu_60,2,INDOUT,NBSIZ1,IADR)
call iDAFILE(Lu_60,2,INDS,NBITM1+2,IADR)
call dDAFILE(Lu_60,2,BUFS,NBITM1,IADR)
!LENGTH = INDOUT(IBBC1)
!IADR = INDOUT(IBDA1)
LENGTH = INDS(NBITM1+1,1)
IADR = INDS(NBITM1+2,1)
do I=1,LENGTH
  IBUFIJ = IBUFIJ+1
  !VALSRT(IBUFIJ) = BUFOUT(I)
  !INDSRT(IBUFIJ) = INDOUT(IBOFF1+I)
  VALSRT(IBUFIJ) = BUFS(I,1)
  INDSRT(IBUFIJ) = INDS(I,1)
  if (IBUFIJ < NSRTMX) GO TO 202
  NSRTCN = NSRTMX
  JDISK = IDISK

  INDSRT(NSRTMX+1) = NSRTCN
  INDSRT(NSRTMX+2) = ISRTAD
  call dDAFILE(Lu_70,1,VALSRT,NSRTMX,IDISK)
  call iDAFILE(Lu_70,1,INDSRT,NSRTMX+2,IDISK)

  ISRTAD = JDISK
  IBUFIJ = 0
202 continue
end do
if (IADR /= -1) GO TO 201
! EMPTY LAST BUFFER
NSRTCN = IBUFIJ
JDISK = IDISK
!
INDSRT(NSRTMX+1) = NSRTCN
INDSRT(NSRTMX+2) = ISRTAD
call dDAFILE(Lu_70,1,VALSRT,NSRTMX,IDISK)
call iDAFILE(Lu_70,1,INDSRT,NSRTMX+2,IDISK)

LASTAD(1) = JDISK
! ABCI
IAD15 = IDISK
IADABCI = IAD15
INSOUT = 0
IADD10 = IAD10(4)
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
IN = 2
NSAVE = ICOP1(IN)
100 NI = NSAVE
IOUT = 0
110 IN = IN+1
if (IN > LEN) then
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  LEN = ICOP1(nCOP+1)
  if (LEN <= 0) GO TO 6
  IN = 1
end if
if (ICHK /= 0) GO TO 460
if (ICOP1(IN) == 0) then
  ICHK = 1
else
  IOUT = IOUT+1
end if
GO TO 110
460 ICHK = 0
NSAVE = ICOP1(IN)
6 continue
NIB = 1+(NI-1)*NVIRT
! LOOP OVER VIRTUAL ORBITAL INDEX B:
do NB=1,NVIRT
  NSIB = MUL(NSM(LN+NB),NSM(NI))
  INS = NVPAIR(NSIB)
  if (INS /= 0) then
    call FZERO(BIAC,INS)
    call FZERO(BICA,INS)
  end if
  NIB = NIB+1
  ! READ & PROCESS INTEGRAL BUFFERS ON UNIT 14:
  IADR = LASTAD(NIB)
203 continue
  !call dDAFILE(Lu_60,2,INDOUT,NBSIZ1,IADR)
  call iDAFILE(Lu_60,2,INDS,NBITM1+2,IADR)
  call dDAFILE(Lu_60,2,BUFS,NBITM1,IADR)
  !LENGTH = INDOUT(IBBC1)
  !IADR = INDOUT(IBDA1)
  LENGTH = INDS(NBITM1+1,1)
  IADR = INDS(NBITM1+2,1)
  do KK=1,LENGTH
    !INND = INDOUT(IBOFF1+KK)
    INND = INDS(KK,1)
    NA = (INND-1)/NVIRT+1
    NC = INND-(NA-1)*NVIRT
    IACS = ISAB(NA+(NC-1)*NVIRT)
    !BIAC(IACS) = BIAC(IACS)+BUFOUT(KK)
    BIAC(IACS) = BIAC(IACS)+BUFS(KK,1)
    !if (NA > NC) BICA(IACS) = BICA(IACS)-BUFOUT(KK)
    !if (NA < NC) BICA(IACS) = BICA(IACS)+BUFOUT(KK)
    if (NA > NC) BICA(IACS) = BICA(IACS)-BUFS(KK,1)
    if (NA < NC) BICA(IACS) = BICA(IACS)+BUFS(KK,1)
  end do
  if (IADR /= -1) GO TO 203
  do ILOOP=0,1
    INSB = INS
73  INB = KBUFF1-INSOUT
    INUMB = INSB
    if (INSB > INB) INUMB = INB
    IST = INS-INSB+1
    if (ILOOP == 0) call DCOPY_(INUMB,BIAC(IST),1,BUFBI(INSOUT+1),1)
    if (ILOOP == 1) call DCOPY_(INUMB,BICA(IST),1,BUFBI(INSOUT+1),1)
    INSOUT = INSOUT+INUMB
    if (INSOUT > KBUFF1) then
      write(6,*) 'SortA: INSOUT > KBUFF1'
      write(6,*) 'INSOUT=',INSOUT
      write(6,*) 'KBUFF1=',KBUFF1
      call ABEND()
    end if
    if (INSOUT == KBUFF1) then
      call dDAFILE(Lu_70,1,BUFBI,KBUFF1,IAD15)
      INSOUT = 0
    end if
    INSB = INSB-INUMB
    if (INSB > 0) GO TO 73
  end do
end do
if (LEN >= 0) GO TO 100
! EMPTY LAST BUFFER IF NOT EMPTY
if (INSOUT > 0) call dDAFILE(Lu_70,1,BUFBI,KBUFF1,IAD15)

return

end subroutine SORTA
