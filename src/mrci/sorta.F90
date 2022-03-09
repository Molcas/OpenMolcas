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

use mrci_global, only: IADABCI, ICH, IFIRST, INDSRT, IPRINT, IROW, ISMAX, KBUFF1, LASTAD, Lu_60, Lu_70, LUSYMB, LUTRA, LN, MCHAIN, &
                       NBITM1, NCHN1, NORB, NSM, NSRTMX, NSYM, NTIBUF, NVIRT, NVPAIR, TIBUF, VALSRT
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: BUFS(NBITM1,NCHN1), BUFBI(KBUFF1), BIAC(ISMAX), BICA(ISMAX)
integer(kind=iwp), intent(out) :: INDS(NBITM1+2,NCHN1), NINTGR
integer(kind=iwp), intent(in) :: ISAB(*)
#include "tratoc.fh"
#include "warnings.h"
integer(kind=iwp) :: I, IACS, IAD15, IAD50, IADD10, IADR, IBUFIJ, ICHK, IDISK, IIJ, IIN, IJ, IJKL, ILEN, ILOOP, INB, INND, INS, &
                     INSB, INSOUT, INUMB, IOUT, IPOS, ISRTAD, IST, JDISK, KK, KL, LENGTH, M, NA, NAT, NB, NC, NI, NIB, NJ, NK, NL, &
                     NOP, NOQ, NOR, NORB0(9), NORBP, NOS, NSAVE, NSIB, NSP, NSPQ, NSPQR, NSQ, NSR, NSRTCN, NSS, NSSM, NT, NTM, NU, &
                     NUMAX, NUMIN, NV, NX, NXM
real(kind=wp) :: FINI
logical(kind=iwp) :: Skip

call COUNT_MRCI(NINTGR,NSYM,NORB,MUL)
if (IPRINT >= 6) write(u6,1234) NINTGR

IAD50 = 0
call iDAFILE(LUTRA,2,iTraToc,nTraToc,IAD50)

IDISK = 0
ICHK = 0
INDS(NBITM1+1,:) = 0
INDS(NBITM1+2,:) = -1
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
                if (NI <= 0) cycle
                NJ = ICH(NORB0(NSQ)+NU)
                if (NJ <= 0) cycle
                NK = ICH(NORB0(NSR)+NV)
                if (NK <= 0) cycle
                NL = ICH(NORB0(NSS)+NX)
                if (NL <= 0) cycle
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
                if (NI > LN) then
                  if (NK <= LN) cycle
                  if (IFIRST /= 0) cycle
                  if (NJ > LN) then
                    if (NL > LN) cycle
                    ! ABCI
                    NA = NI-LN
                    NB = NJ-LN
                    NC = NK-LN
                    NI = NL
                  else
                    if (NL <= LN) cycle
                    ! CIAB
                    NA = NK-LN
                    NB = NL-LN
                    NC = NI-LN
                    NI = NJ
                  end if
                  NIB = (NI-1)*NVIRT+NB+1
                  IPOS = INDS(NBITM1+1,NIB)+1
                  INDS(NBITM1+1,NIB) = IPOS
                  BUFS(IPOS,NIB) = FINI
                  INDS(IPOS,NIB) = (NA-1)*NVIRT+NC
                  if (IPOS >= NBITM1) then
                    JDISK = IDISK
                    call iDAFILE(Lu_60,1,INDS(:,NIB),NBITM1+2,IDISK)
                    call dDAFILE(Lu_60,1,BUFS(:,NIB),NBITM1,IDISK)
                    INDS(NBITM1+1,NIB) = 0
                    INDS(NBITM1+2,NIB) = JDISK
                  end if
                  if (NA /= NB) then
                    NAT = NA
                    NA = NB
                    NB = NAT
                    NIB = (NI-1)*NVIRT+NB+1
                    IPOS = INDS(NBITM1+1,NIB)+1
                    INDS(NBITM1+1,NIB) = IPOS
                    BUFS(IPOS,NIB) = FINI
                    INDS(IPOS,NIB) = (NA-1)*NVIRT+NC
                    if (IPOS >= NBITM1) then
                      JDISK = IDISK
                      call iDAFILE(Lu_60,1,INDS(:,NIB),NBITM1+2,IDISK)
                      call dDAFILE(Lu_60,1,BUFS(:,NIB),NBITM1,IDISK)
                      INDS(NBITM1+1,NIB) = 0
                      INDS(NBITM1+2,NIB) = JDISK
                    end if
                  end if
                else
                  ! IJKL
                  IIJ = IROW(NI)+NJ
                  KL = IROW(NK)+NL
                  IJKL = IIJ*(IIJ-1)/2+KL
                  IJ = 1
                  IPOS = INDS(NBITM1+1,IJ)+1
                  INDS(NBITM1+1,IJ) = IPOS
                  BUFS(IPOS,IJ) = FINI
                  INDS(IPOS,IJ) = IJKL
                  if (IPOS == NBITM1) then
                    JDISK = IDISK
                    call iDAFILE(Lu_60,1,INDS(:,IJ),NBITM1+2,IDISK)
                    call dDAFILE(Lu_60,1,BUFS(:,IJ),NBITM1,IDISK)
                    INDS(NBITM1+1,IJ) = 0
                    INDS(NBITM1+2,IJ) = JDISK
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
! EMPTY LAST BUFFERS
if (NCHN1 > MCHAIN) then
  write(u6,*) 'SORTA Error: NCHN1 > MCHAIN (See code).'
  call QUIT(_RC_GENERAL_ERROR_)
end if
do I=1,NCHN1
  JDISK = IDISK
  call iDAFILE(Lu_60,1,INDS(:,I),NBITM1+2,IDISK)
  call dDAFILE(Lu_60,1,BUFS(:,I),NBITM1,IDISK)
  LASTAD(I) = JDISK
end do
! IJKL
IDISK = 0
IBUFIJ = 0
ISRTAD = -1
IADR = LASTAD(1)
do
  call iDAFILE(Lu_60,2,INDS,NBITM1+2,IADR)
  call dDAFILE(Lu_60,2,BUFS,NBITM1,IADR)
  LENGTH = INDS(NBITM1+1,1)
  IADR = INDS(NBITM1+2,1)
  do I=1,LENGTH
    IBUFIJ = IBUFIJ+1
    VALSRT(IBUFIJ) = BUFS(I,1)
    INDSRT(IBUFIJ) = INDS(I,1)
    if (IBUFIJ < NSRTMX) cycle
    NSRTCN = NSRTMX
    JDISK = IDISK

    INDSRT(NSRTMX+1) = NSRTCN
    INDSRT(NSRTMX+2) = ISRTAD
    call dDAFILE(Lu_70,1,VALSRT,NSRTMX,IDISK)
    call iDAFILE(Lu_70,1,INDSRT,NSRTMX+2,IDISK)

    ISRTAD = JDISK
    IBUFIJ = 0
  end do
  if (IADR == -1) exit
end do
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
      call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      ILEN = ICOP1(nCOP+1)
      if (ILEN <= 0) then
        Skip = .true.
        exit
      end if
      IIN = 1
    end if
    if (ICHK /= 0) exit
    if (ICOP1(IIN) == 0) then
      ICHK = 1
    else
      IOUT = IOUT+1
    end if
  end do
  if (.not. Skip) then
    ICHK = 0
    NSAVE = ICOP1(IIN)
  end if
  NIB = 1+(NI-1)*NVIRT
  ! LOOP OVER VIRTUAL ORBITAL INDEX B:
  do NB=1,NVIRT
    NSIB = MUL(NSM(LN+NB),NSM(NI))
    INS = NVPAIR(NSIB)
    if (INS /= 0) then
      BIAC(1:INS) = Zero
      BICA(1:INS) = Zero
    end if
    NIB = NIB+1
    ! READ & PROCESS INTEGRAL BUFFERS ON UNIT 14:
    IADR = LASTAD(NIB)
    do
      call iDAFILE(Lu_60,2,INDS,NBITM1+2,IADR)
      call dDAFILE(Lu_60,2,BUFS,NBITM1,IADR)
      LENGTH = INDS(NBITM1+1,1)
      IADR = INDS(NBITM1+2,1)
      do KK=1,LENGTH
        INND = INDS(KK,1)
        NA = (INND-1)/NVIRT+1
        NC = INND-(NA-1)*NVIRT
        IACS = ISAB(NA+(NC-1)*NVIRT)
        BIAC(IACS) = BIAC(IACS)+BUFS(KK,1)
        if (NA > NC) BICA(IACS) = BICA(IACS)-BUFS(KK,1)
        if (NA < NC) BICA(IACS) = BICA(IACS)+BUFS(KK,1)
      end do
      if (IADR == -1) exit
    end do
    do ILOOP=0,1
      INSB = INS
      do
        INB = KBUFF1-INSOUT
        INUMB = INSB
        if (INSB > INB) INUMB = INB
        IST = INS-INSB+1
        if (ILOOP == 0) call DCOPY_(INUMB,BIAC(IST),1,BUFBI(INSOUT+1),1)
        if (ILOOP == 1) call DCOPY_(INUMB,BICA(IST),1,BUFBI(INSOUT+1),1)
        INSOUT = INSOUT+INUMB
        if (INSOUT > KBUFF1) then
          write(u6,*) 'SortA: INSOUT > KBUFF1'
          write(u6,*) 'INSOUT=',INSOUT
          write(u6,*) 'KBUFF1=',KBUFF1
          call ABEND()
        end if
        if (INSOUT == KBUFF1) then
          call dDAFILE(Lu_70,1,BUFBI,KBUFF1,IAD15)
          INSOUT = 0
        end if
        INSB = INSB-INUMB
        if (INSB <= 0) exit
      end do
    end do
  end do
  if (ILEN < 0) exit
end do
! EMPTY LAST BUFFER IF NOT EMPTY
if (INSOUT > 0) call dDAFILE(Lu_70,1,BUFBI,KBUFF1,IAD15)

return

1234 format(//6X,'NUMBER OF TWO-ELECTRON INTEGRALS',I10)

end subroutine SORTA
