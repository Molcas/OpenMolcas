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

subroutine SORT_MRCI(BUFS,INDS,FC,FIIJJ,FIJIJ)

use mrci_global, only: IAD25S, ICH, IPRINT, IROW, ITOC17, LASTAD, LN, Lu_25, Lu_60, LUONE, LUTRA, MCHAIN, NBITM3, NBTRI, NCHN3, &
                       NELEC, NORB, NORBT, NSM, NSYM, NTIBUF, NVIR, NVIRP, NVIRT, POTNUC, TIBUF
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: BUFS(NBITM3,NCHN3), FC(NBTRI)
integer(kind=iwp), intent(out) :: INDS(NBITM3+2,NCHN3)
real(kind=wp), intent(_OUT_) :: FIIJJ(*), FIJIJ(*)
#include "tratoc.fh"
#include "warnings.h"
integer(kind=iwp) :: I, IAD50, IADD17, IADD25, IBUF, IDISK, IEXP, IIJ, IIN, IJ, IJT, IKT, INAV, IND, IORBI, IOUT, IPOF(65), IPOS, &
                     ISYM, IVEC(20), J, JDISK, JK, JORBI, KORBI, M1, M2, M3, M4, N1, N2, N3, N4, NAV, NBV, NI, NJ, NK, NL, NOP, &
                     NOQ, NOR, NORB0(9), NORBP, NORBTT, NOS, NOT2, NOTT, NOVST, NSA, NSB, NSIJT, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, &
                     NSSM, NT, NTM, NTMP, NU, NUMAX, NUMIN, NV, NVT, NX, NXM
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

INDS(NBITM3+1,:) = 0
INDS(NBITM3+2,:) = -1
NORB0(1) = 0
do I=1,NSYM
  NORB0(I+1) = NORB0(I)+NORB(I)
end do
! READ ONE-ELECTRON ORBITALS. USE FIIJJ TEMPORARILY AS READ BUFFER.
NORBTT = 0
do ISYM=1,NSYM
  NORBTT = NORBTT+(NORB(ISYM)*(NORB(ISYM)+1))/2
end do
EFROZ = POTNUC
FC(:) = Zero
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
        EFROZ = EFROZ+2*ONEHAM
      end if
    end do
  end do
  KORBI = KORBI+NORB(ISYM)
end do
if (IPRINT >= 20) then
  call TRIPRT('FC IN SORT_MRCI BEFORE TWOEL',' ',FC,NORBT)
  write(u6,'(A,F20.8)') ' EFROZ:',EFROZ
end if
FIIJJ(1:NBTRI) = Zero
FIJIJ(1:NBTRI) = Zero
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
                M1 = ICH(NORB0(NSP)+NT)
                M2 = ICH(NORB0(NSQ)+NU)
                M3 = ICH(NORB0(NSR)+NV)
                M4 = ICH(NORB0(NSS)+NX)
                if ((M1 == 0) .or. (M2 == 0)) cycle
                if ((M3 == 0) .or. (M4 == 0)) cycle
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
                if ((NI <= 0) .or. (NJ <= 0) .or. (NK <= 0) .or. (NL <= 0)) then
                  ! CHECK FOR FOCK-MATRIX, AND FROZEN ENERGY, CONTRIBUTIONS
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
                  cycle
                end if
                DFINI = abs(FINI)+1.0e-20_wp
                IEXP = int(-log10(DFINI))+5
                if (IEXP <= 20) IVEC(IEXP) = IVEC(IEXP)+1
                if ((NI == NJ) .and. (NK == NL)) then
                  IJ = IROW(NI)+NK
                  FIIJJ(IJ) = FINI
                  ! SKIP (AA/II) INTEGRALS
                else
                  if ((NI == NK) .and. (NJ == NL)) then
                    IJ = IROW(NI)+NJ
                    FIJIJ(IJ) = FINI
                  end if
                  if (NI <= LN) then
                  else if (NJ > LN) then
                    if (NK > LN) cycle
                    ! ABIJ
                    IIJ = IROW(NK)+NL
                    IPOS = INDS(NBITM3+1,IIJ)+1
                    INDS(NBITM3+1,IIJ) = IPOS
                    BUFS(IPOS,IIJ) = FINI
                    NSA = NSM(NI)
                    NAV = NI-LN-NVIRP(NSA)
                    NSB = NSM(NJ)
                    NBV = NJ-LN-NVIRP(NSB)
                    NSIJT = (MUL(NSM(NK),NSM(NL))-1)*NSYM
                    INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                    INDS(IPOS,IIJ) = INAV
                    if (IPOS >= NBITM3) then
                      JDISK = IDISK
                      call iDAFILE(Lu_60,1,INDS(1,IIJ),NBITM3+2,IDISK)
                      call dDAFILE(Lu_60,1,BUFS(1,IIJ),NBITM3,IDISK)
                      INDS(NBITM3+1,IIJ) = 0
                      INDS(NBITM3+2,IIJ) = JDISK
                    end if
                    if (NK /= NL) cycle
                    if (NI == NJ) cycle
                    IJT = IROW(NI)+NJ
                    FC(IJT) = FC(IJT)+2*FINI
                  else if (NK > LN) then
                    if (NL > LN) cycle
                    ! AIBJ
                    IIJ = NOT2+IROW(NJ)+NL
                    if (NL > NJ) IIJ = NOT2+IROW(NL)+NJ
                    IPOS = INDS(NBITM3+1,IIJ)+1
                    INDS(NBITM3+1,IIJ) = IPOS
                    BUFS(IPOS,IIJ) = FINI
                    NSA = NSM(NI)
                    NAV = NI-LN-NVIRP(NSA)
                    NSB = NSM(NK)
                    NBV = NK-LN-NVIRP(NSB)
                    NSIJT = (MUL(NSM(NJ),NSM(NL))-1)*NSYM
                    if (NL <= NJ) then
                      INAV = IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
                    else
                      INAV = IPOF(NSIJT+NSB)+(NAV-1)*NVIR(NSB)+NBV
                    end if
                    INDS(IPOS,IIJ) = INAV
                    if (IPOS >= NBITM3) then
                      JDISK = IDISK
                      call iDAFILE(Lu_60,1,INDS(1,IIJ),NBITM3+2,IDISK)
                      call dDAFILE(Lu_60,1,BUFS(1,IIJ),NBITM3,IDISK)
                      INDS(NBITM3+1,IIJ) = 0
                      INDS(NBITM3+2,IIJ) = JDISK
                    end if
                    if (NJ /= NL) cycle
                    if (NI == NK) cycle
                    IKT = IROW(NI)+NK
                    FC(IKT) = FC(IKT)-FINI
                  else
                    ! AIJK
                    JK = NOTT+IROW(NK)+NL
                    IPOS = INDS(NBITM3+1,JK)+1
                    INDS(NBITM3+1,JK) = IPOS
                    BUFS(IPOS,JK) = FINI
                    INDS(IPOS,JK) = IROW(NI)+NJ
                    if (IPOS >= NBITM3) then
                      JDISK = IDISK
                      call iDAFILE(Lu_60,1,INDS(1,JK),NBITM3+2,IDISK)
                      call dDAFILE(Lu_60,1,BUFS(1,JK),NBITM3,IDISK)
                      INDS(NBITM3+1,JK) = 0
                      INDS(NBITM3+2,JK) = JDISK
                    end if
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
if ((NOVST+NCHN3) > MCHAIN) then
  write(u6,*) 'SORT_MRCI Error: NOVST+NCHN3>MCHAIN (See code).'
  call QUIT(_RC_GENERAL_ERROR_)
end if
do I=1,NCHN3
  JDISK = IDISK
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
write(u6,155) (IVEC(I),I=1,20)
!end if

return

154 format(//6X,'STATISTICS FOR INTEGRALS, FIRST ENTRY 10**3-10**4',/)
155 format(6X,5I10)

end subroutine SORT_MRCI
