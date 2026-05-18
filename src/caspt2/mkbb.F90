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

subroutine MKBB(DREF,NDREF,PREF,NPREF,FD,FP)

use definitions, only: iwp, wp
use constants, only: Half, Two, Four, Eight
use SUPERINDEX, only: MTU, MTGEU, KTU, KTGTU
use caspt2_global, only: ipea_shift
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT, IDBMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NTU, NTUES, NASHT, EASUM, EPSA, NTGEU, NTGEUES, NTGTU, NTGTUES

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(in) :: FD(NDREF), FP(NPREF)
real(kind=wp), allocatable :: BB(:), BBP(:), SP(:), SDP(:), BBM(:), SM(:), SDM(:)
integer(kind=iwp) ISYM, NINP, NAS, NBB, ITU, ITUABS, ITABS, IUABS, IXY, IXYABS, IXABS, IYABS, IBADR, IXT, IYU, IP1, IP2, IP, I, &
                  IBMADR, IBPADR, ID, ID1, ID2, IDIAG, IDISK, IDSM, IDSP, IDT, IDU, INSM, ITGEU, ITGEUABS, ITGTU, IXGEY, IXGEYABS, &
                  IXGTY, IYX, NASM, NASP, NBBM, NBBP, NSM, NSP
real(kind=wp) ET, EU, EX, EY, ATUXY, ATUX, ATYU, ATUY, ATYX, BTUYX, BTUXY, value

! Set up the matrices BBP(tu,xy) and BBM(tu,xy)
! Formulae used:
!    BB(tu,xy)= 2*( Fyuxt - (A-Et-Eu-Ex-Ey)*Gyuxt )
!      + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
!      + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
!      - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
!      - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
!      + 8*dxt*dyu (Et+Ey)
!      - 4*dxu*dyt (Et+Ex)
! where A= EASUM= sum over active w of (Ew*Dww).
!    BBP(tu,xy)=BB(tu,xy)+BB(tu,yx)
!    BBM(tu,xy)=BB(tu,xy)-BB(tu,yx)

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NINP = NINDEP(ISYM,2)
  if (NINP == 0) cycle
  NAS = NTU(ISYM)
  NBB = (NAS*(NAS+1))/2
  if (NBB > 0) call mma_allocate(BB,NBB,Label='BB')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    ET = EPSA(ITABS)
    EU = EPSA(IUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      EX = EPSA(IXABS)
      EY = EPSA(IYABS)
      IBADR = (ITU*(ITU-1))/2+IXY
      IXT = IXABS+NASHT*(ITABS-1)
      IYU = IYABS+NASHT*(IUABS-1)
      IP1 = max(IXT,IYU)
      IP2 = min(IXT,IYU)
      IP = (IP1*(IP1-1))/2+IP2
      ATUXY = EASUM-ET-EU-EX-EY
      value = Four*(FP(IP)-ATUXY*PREF(IP))
      ! Add  + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
      if (IXABS == ITABS) then
        ID1 = max(IYABS,IUABS)
        ID2 = min(IYABS,IUABS)
        ID = (ID1*(ID1-1))/2+ID2
        ATYU = EASUM-ET-EY-EU
        value = value+Four*(ATYU*DREF(ID)-FD(ID))
        ! Add  + 8*dxt*dyu (Et+Ey)
        if (IYABS == IUABS) value = value+Eight*(ET+EY)
      end if
      ! Add  + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
      if (IYABS == IUABS) then
        ID1 = max(IXABS,ITABS)
        ID2 = min(IXABS,ITABS)
        ID = (ID1*(ID1-1))/2+ID2
        ATYX = EASUM-ET-EY-EX
        value = value+Four*(ATYX*DREF(ID)-FD(ID))
      end if
      ! Add  - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
      if (IYABS == ITABS) then
        ID1 = max(IXABS,IUABS)
        ID2 = min(IXABS,IUABS)
        ID = (ID1*(ID1-1))/2+ID2
        ATUX = EASUM-ET-EU-EX
        value = value-Two*(ATUX*DREF(ID)-FD(ID))
        ! Add  - 4*dxu*dyt (Et+Ex)
        if (IXABS == IUABS) value = value-Four*(ET+EX)
      end if
      ! Add  - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
      if (IXABS == IUABS) then
        ID1 = max(IYABS,ITABS)
        ID2 = min(IYABS,ITABS)
        ID = (ID1*(ID1-1))/2+ID2
        ATUY = EASUM-ET-EU-EY
        value = value-Two*(ATUY*DREF(ID)-FD(ID))
      end if
      BB(IBADR) = value
    end do
  end do
  NASP = NTGEU(ISYM)
  NBBP = (NASP*(NASP+1))/2
  if (NBBP > 0) then
    call mma_allocate(BBP,NBBP,Label='BBP')
    !GG.Nov03  Load in SDP the diagonal elements of SBP matrix:
    NSP = (NASP*(NASP+1))/2
    call mma_allocate(SP,NSP,Label='SP')
    call mma_allocate(SDP,NASP,Label='SDP')
    IDSP = IDSMAT(ISYM,2)
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
  NBBM = (NASM*(NASM+1))/2
  if (NBBM > 0) then
    call mma_allocate(BBM,NBBM,Label='BBM')
    !GG.Nov03  Load in SDM the diagonal elements of SBM matrix:
    NSM = (NASM*(NASM+1))/2
    call mma_allocate(SDM,NASM,Label='SDM')
    if (NINDEP(ISYM,3) > 0) then
      call mma_allocate(SM,NSM,Label='SM')
      IDSM = IDSMAT(ISYM,3)
      call DDAFILE(LUSBT,2,SM,NSM,IDSM)
      IDIAG = 0
      do I=1,NASM
        IDIAG = IDIAG+I
        SDM(I) = SM(IDIAG)
      end do
      call mma_deallocate(SM)
    end if
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
      BTUXY = BB(IBADR)
      if (ITU >= IYX) then
        IBADR = (ITU*(ITU-1))/2+IYX
      else
        IBADR = (IYX*(IYX-1))/2+ITU
      end if
      BTUYX = BB(IBADR)
      IBPADR = (ITGEU*(ITGEU-1))/2+IXGEY
      BBP(IBPADR) = BTUXY+BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = (ITABS*(ITABS+1))/2
        IDU = (IUABS*(IUABS+1))/2
        BBP(IBPADR) = BBP(IBPADR)+ipea_shift*half*(DREF(IDT)+DREF(IDU))*SDP(ITGEU)
      end if
      !GG End
      if (NINDEP(ISYM,3) < 1) cycle
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      IBMADR = (ITGTU*(ITGTU-1))/2+IXGTY
      BBM(IBMADR) = BTUXY-BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = (ITABS*(ITABS+1))/2
        IDU = (IUABS*(IUABS+1))/2
        BBM(IBMADR) = BBM(IBMADR)+ipea_shift*half*(DREF(IDT)+DREF(IDU))*SDM(INSM)
        INSM = INSM+1
      end if
      !GG.End
    end do
  end do
  if (NBB > 0) call mma_deallocate(BB)

  ! Write to disk, and save size and address.
  if ((NBBP > 0) .and. (NINDEP(ISYM,2) > 0)) then
    IDISK = IDBMAT(ISYM,2)
    call DDAFILE(LUSBT,1,BBP,NBBP,IDISK)
    call mma_deallocate(BBP)
    !GG.Nov03 DisAlloc SDP
    call mma_deallocate(SDP)
    !GG End
  end if
  if (NBBM > 0) then
    if (NINDEP(ISYM,3) > 0) then
      IDISK = IDBMAT(ISYM,3)
      call DDAFILE(LUSBT,1,BBM,NBBM,IDISK)
    end if
    call mma_deallocate(BBM)
    !GG.Nov03 DisAlloc SDM
    call mma_deallocate(SDM)
    !GG End
  end if
end do

end subroutine MKBB
