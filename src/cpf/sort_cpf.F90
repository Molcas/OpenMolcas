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

use cpf_global, only: IAD25S, ICH, IPRINT, IROW, ITOC17, LASTAD, LBUF, LN, Lu_25, Lu_TiABIJ, Lu_TraInt, Lu_TraOne, MADR, N, NORB, &
                      NORBT, NSM, NSYM, NSYS, NTIBUF, NVIR, NVIRT, POTNUC
use Symmetry_Info, only: Mul
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6, RtoI

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: BUFOUT(*), FC(*), FIJ(*), FJI(*), TIBUF(NTIBUF)
integer(kind=iwp), intent(_OUT_) :: INDOUT(*), ICAD(*), IBUFL(*)
#include "tratoc.fh"
integer(kind=iwp) :: I, IAD50, IADD17, IADD25, IBUF, ICP, ICPP, ICQ, ID, IDISK, IDIV, IEXP, II, IIJ, IIN, IJ, IJT, INAV, IND, &
                     IORBI, IOUT, IPOF(65), IREC, ISYM, IVEC(20), J, JDISK, JK, JNAV, JORBI, KORBI, LBUF0, LBUF1, LBUF2, M1, M2, &
                     M3, M4, N1, N2, N3, N4, NAV, NBV, NI, NJ, NK, NL, NOB2, NOP, NOQ, NOR, NORB0(9), NORBP, NORBTT, NOS, NOT2, &
                     NOTT, NOV, NOVST, NSA, NSB, NSIJT, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NSSM, NT, NTM, NTMP, NU, NUMAX, NUMIN, &
                     NV, NVT, NX, NXM
real(kind=wp) :: DFINI, EMY, FINI, ONEHAM

IAD50 = 0
call iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
NVT = IROW(NVIRT+1)
do I=1,20
  IVEC(I) = 0
end do
IIN = 1
do I=1,NSYM
  call IPO_CPF(IPOF(IIN),NVIR,MUL,NSYM,I,-1)
  IIN = IIN+NSYM
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
FC(1:NOB2) = Zero
IBUF = 0
KORBI = 0
do ISYM=1,NSYM
  do JORBI=KORBI+1,KORBI+NORB(ISYM)
    do IORBI=KORBI+1,JORBI
      IBUF = IBUF+1
      ONEHAM = FIJ(IBUF)
      NI = ICH(IORBI)
      NJ = ICH(JORBI)
      if ((NI == 0) .or. (NJ == 0)) cycle
      if (NI < NJ) then
        NTMP = NI
        NI = NJ
        NJ = NTMP
      end if
      if (NJ > 0) then
        IJT = IROW(NI)+NJ
        FC(IJT) = FC(IJT)+ONEHAM
      else if (NI == NJ) then
        EMY = EMY+Two*ONEHAM
      end if
    end do
  end do
  KORBI = KORBI+NORB(ISYM)
end do
FIJ(1:NOB2) = Zero
FJI(1:NOB2) = Zero
if (IPRINT >= 20) then
  call TRIPRT('FC IN SORT BEFORE TWOEL',' ',FC,NORBT)
  write(u6,'(A,F20.8)') ' EMY:',EMY
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
                if ((M1 == 0) .or. (M2 == 0) .or. (M3 == 0) .or. (M4 == 0)) cycle
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
                if ((NI > 0) .and. (NJ > 0) .and. (NK > 0) .and. (NL > 0)) then
                  DFINI = abs(FINI)
                  IEXP = int(-log10(DFINI+1.0e-20_wp)+5)
                  if (IEXP > 20) IEXP = 20
                  if (IEXP < 1) IEXP = 1
                  IVEC(IEXP) = IVEC(IEXP)+1
                  if ((NI == NJ) .and. (NK == NL)) then
                    IJ = IROW(NI)+NK
                    FIJ(IJ) = FINI
                    ! SKIP (AA/II) INTEGRALS
                  else
                    if ((NI == NK) .and. (NJ == NL)) then
                      IJ = IROW(NI)+NJ
                      FJI(IJ) = FINI
                    end if
                    if (NI >= LN) then
                      if (NJ <= LN) then
                        if (NK <= LN) then
                          ! AIJK
                          JK = NOTT+IROW(NK)+NL
                          IBUFL(JK) = IBUFL(JK)+1
                          ICQ = ICAD(JK)
                          ICP = ICQ/IDIV+IBUFL(JK)
                          BUFOUT(ICP) = FINI
                          ICPP = ICQ+LBUF0+IBUFL(JK)
                          INDOUT(ICPP) = IROW(NI)+NJ
                          if (IBUFL(JK) >= LBUF) then
                            INDOUT(ICQ+LBUF1) = LBUF
                            JDISK = IDISK
                            call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
                            INDOUT(ICQ+LBUF2) = JDISK
                            IBUFL(JK) = 0
                          end if
                        else if (NL <= LN) then
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
                          if (NL <= NJ) then
                            INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                          else
                            INAV = IPOF(NSIJT+NSB)+(NAV-1)*NVIR(NSB)+NBV
                          end if
                          ICPP = ICQ+LBUF0+IBUFL(IIJ)
                          INDOUT(ICPP) = INAV
                          if (IBUFL(IIJ) >= LBUF) then
                            INDOUT(ICQ+LBUF1) = LBUF
                            JDISK = IDISK
                            call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
                            INDOUT(ICQ+LBUF2) = JDISK
                            IBUFL(IIJ) = 0
                          end if
                          if ((NJ == NL) .and. (NI /= NK)) then
                            JNAV = IROW(NI)+NK
                            FC(JNAV) = FC(JNAV)-FINI
                          end if
                        end if
                      else if (NK <= LN) then
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
                        if (IBUFL(IIJ) >= LBUF) then
                          INDOUT(ICQ+LBUF1) = LBUF
                          JDISK = IDISK
                          call iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
                          INDOUT(ICQ+LBUF2) = JDISK
                          IBUFL(IIJ) = 0
                        end if
                        if ((NK == NL) .and. (NI /= NJ)) then
                          JNAV = IROW(NI)+NJ
                          FC(JNAV) = FC(JNAV)+Two*FINI
                        end if
                      end if
                    end if
                  end if
                else
                  ! CHECK FOR FOCK-MATRIX CONTRIBUTION
                  if (NI == NJ) then
                    II = 1
                    call IFOCK(FC,NI,NK,NL,FINI,II)
                    if (NK == NL) then
                      if ((NI <= 0) .and. (NK <= 0)) then
                        EMY = EMY+Two*FINI
                        if (NI /= NK) EMY = EMY+Two*FINI
                      end if
                      if (NI /= NK) call IFOCK(FC,NK,NI,NJ,FINI,II)
                    end if
                  else if (NK == NL) then
                    II = 1
                    call IFOCK(FC,NK,NI,NJ,FINI,II)
                  end if
                  II = 0
                  if (NI == NK) then
                    call IFOCK(FC,NI,NJ,NL,FINI,II)
                    if (NJ == NL) then
                      if ((NI <= 0) .and. (NJ <= 0)) then
                        EMY = EMY-FINI
                        if (NI /= NJ) EMY = EMY-FINI
                      end if
                      if (NI /= NJ) call IFOCK(FC,NJ,NI,NK,FINI,II)
                    end if
                  else if (NI == NL) then
                    call IFOCK(FC,NI,NJ,NK,FINI,II)
                  else if (NJ == NK) then
                    call IFOCK(FC,NJ,NI,NL,FINI,II)
                  else if (NJ == NL) then
                    call IFOCK(FC,NJ,NI,NK,FINI,II)
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
if ((NOVST+NOV) > MADR) then
  write(u6,*) 'SORT Error: NOVST+NOV>MADR (See code).'
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
write(u6,154)
write(u6,155) (IVEC(I),I=1,20)

return

154 format(//6X,'STATISTICS FOR INTEGRALS, FIRST ENTRY 10**3-10**4',/)
155 format(6X,5I10)

end subroutine SORT_CPF
