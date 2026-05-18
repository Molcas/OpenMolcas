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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKBF(DREF,NDREF,PREF,NPREF,FP)

use definitions, only: iwp, wp
use constants, only: Half, Four
use SUPERINDEX, only: MTU, MTGEU, KTU, KTGTU
use caspt2_global, only: ipea_shift
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT, IDBMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NTU, NTUES, NASHT, NTGEU, EASUM, NTGTU, NTGEUES, NTGTUES

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp) PREF(NPREF), FP(NPREF), DREF(NDREF)
real(kind=wp), allocatable :: BF(:), BFP(:), SDP(:), SP(:), BFM(:), SDM(:), SM(:)
integer(kind=iwp) ISYM, NINP, NAS, NBF, ITU, ITUABS, ITABS, IUABS, IXY, IXYABS, IXABS, IYABS, IBADR, ITX, IUY, IP1, IP2, IP, NASP, &
                  NBFP, NSP, IDSP, IDIAG, I, NASM, NBFM, NSM, IDSM, IBMADR, IBPADR, IDISK, IDT, IDU, INSM, ITGEU, ITGEUABS, ITGTU, &
                  IXGEY, IXGEYABS, IXGTY, IYX
real(kind=wp) BTUXY, BTUYX

! Set up the matrices BFP(tu,xy) and BFM(tu,xy)
! Formulae used:
!    BF(tu,xy)= 2*(Ftxuy - EASUM*Gtxuy)
!    BFP(tu,xy)=BF(tu,xy)+BF(tu,yx)
!    BFM(tu,xy)=BF(tu,xy)-BF(tu,yx)

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NINP = NINDEP(ISYM,8)
  if (NINP == 0) cycle
  NAS = NTU(ISYM)
  NBF = (NAS*(NAS+1))/2
  if (NBF > 0) call mma_allocate(BF,NBF,LABEL='BF')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      IBADR = (ITU*(ITU-1))/2+IXY
      ITX = ITABS+NASHT*(IXABS-1)
      IUY = IUABS+NASHT*(IYABS-1)
      IP1 = max(ITX,IUY)
      IP2 = min(ITX,IUY)
      IP = (IP1*(IP1-1))/2+IP2
      BF(IBADR) = Four*(FP(IP)-EASUM*PREF(IP))
    end do
  end do
  NASP = NTGEU(ISYM)
  NBFP = (NASP*(NASP+1))/2
  if (NBFP > 0) then
    call mma_allocate(BFP,NBFP,Label='BFP')
    !GG.Nov03  Load in SDP the diagonal elements of SFP matrix:
    NSP = (NASP*(NASP+1))/2
    call mma_allocate(SP,NSP,Label='SP')
    call mma_allocate(SDP,NASP,Label='SDP')
    IDSP = IDSMAT(ISYM,8)
    call DDAFILE(LUSBT,2,SP,NSP,IDSP)
    IDIAG = 0
    do I=1,NASP
      IDIAG = IDIAG+I
      SDP(I) = SP(IDIAG)
    end do
    call mma_deallocate(SP)
    !GG End
  end if
  NASM = NTGTU(ISYM)
  NBFM = (NASM*(NASM+1))/2
  if (NBFM > 0) then
    call mma_allocate(BFM,NBFM,Label='BFM')
    !GG.Nov03  Load in SDM the diagonal elements of SFM matrix:
    NSM = (NASM*(NASM+1))/2
    call mma_allocate(SM,NSM,Label='SM')
    call mma_allocate(SDM,NASM,Label='SDM')
    IDSM = IDSMAT(ISYM,9)
    call DDAFILE(LUSBT,2,SM,NSM,IDSM)
    IDIAG = 0
    do I=1,NASM
      IDIAG = IDIAG+I
      SDM(I) = SM(IDIAG)
    end do
    call mma_deallocate(SM)
    !GG End
  end if
  INSM = 1
  do ITGEU=1,NASP
    ITGEUABS = ITGEU+NTGEUES(ISYM)
    ITABS = MTGEU(1,ITGEUABS)
    IUABS = MTGEU(2,ITGEUABS)
    ITU = KTU(ITABS,IUABS)-NTUES(ISYM)
    do IXGEY=1,ITGEU
      IXGEYABS = IXGEY+NTGEUES(ISYM)
      IXABS = MTGEU(1,IXGEYABS)
      IYABS = MTGEU(2,IXGEYABS)
      IXY = KTU(IXABS,IYABS)-NTUES(ISYM)
      IYX = KTU(IYABS,IXABS)-NTUES(ISYM)
      if (ITU >= IXY) then
        IBADR = (ITU*(ITU-1))/2+IXY
      else
        IBADR = (IXY*(IXY-1))/2+ITU
      end if
      BTUXY = BF(IBADR)
      if (ITU >= IYX) then
        IBADR = (ITU*(ITU-1))/2+IYX
      else
        IBADR = (IYX*(IYX-1))/2+ITU
      end if
      BTUYX = BF(IBADR)
      IBPADR = (ITGEU*(ITGEU-1))/2+IXGEY
      BFP(IBPADR) = BTUXY+BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = (ITABS*(ITABS+1))/2
        IDU = (IUABS*(IUABS+1))/2
        BFP(IBPADR) = BFP(IBPADR)+ipea_shift*Half*(Four-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
      end if
      !GG End
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      IBMADR = (ITGTU*(ITGTU-1))/2+IXGTY
      BFM(IBMADR) = BTUXY-BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = (ITABS*(ITABS+1))/2
        IDU = (IUABS*(IUABS+1))/2
        BFM(IBMADR) = BFM(IBMADR)+ipea_shift*Half*(Four-DREF(IDT)-DREF(IDU))*SDM(INSM)
        INSM = INSM+1
      end if

    end do
  end do
  if (NBF > 0) call mma_deallocate(BF)

  ! Write to disk
  if ((NBFP > 0) .and. (NINDEP(ISYM,8) > 0)) then
    IDISK = IDBMAT(ISYM,8)
    call DDAFILE(LUSBT,1,BFP,NBFP,IDISK)
    call mma_deallocate(BFP)
    !GG.Nov03 DisAlloc SDP
    call mma_deallocate(SDP)
    !GG End
  end if
  if (NBFM > 0) then
    if (NINDEP(ISYM,9) > 0) then
      IDISK = IDBMAT(ISYM,9)
      call DDAFILE(LUSBT,1,BFM,NBFM,IDISK)
    end if
    call mma_deallocate(BFM)
    !GG.Nov03 DisAlloc SDM
    call mma_deallocate(SDM)
    !GG End
  end if
end do

end subroutine MKBF
