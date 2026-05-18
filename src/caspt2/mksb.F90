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

!***********************************************************************
! Case B (ICASE=2,3)
!***********************************************************************
subroutine MKSB(DREF,NDREF,PREF,NPREF)

use definitions, only: iwp, wp
use constants, only: Two, Four, Eight
use SUPERINDEX, only: MTU, MTGEU, KTU, KTGTU
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NASHT, NSYM, NINDEP, NTU, NTUES, NTGEU, NTGTU, NTGEUES, NTGTUES

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), allocatable :: SB(:), SBP(:), SBM(:)
integer(kind=iwp) ISYM, NINP, NAS, NSB, ITUABS, ITABS, IUABS, IXY, IXYABS, IXABS, IYABS, ISADR, IXT, IYU, IP1, IP2, IP, ID1, ID2, &
                  ID, IDISK, ISMADR, ISPADR, ITGEU, ITGEUABS, ITGTU, ITU, IXGEY, IXGEYABS, IXGTY, IYX, NASM, NASP, NSBM, NSBP
real(kind=wp) value, STUXY, STUYX

! Set up the matrices SBP(tu,xy) and SBM(tu,xy)
! Formulae used:
!    SB(tu,xy)= 4 Pxtyu -4dxt Dyu -4dyu Dxt +2dyt Dxu + 8 dxt dyu -4dxu dyt + 2dxu Dyt
!    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
!    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NINP = NINDEP(ISYM,2)
  if (NINP == 0) cycle
  NAS = NTU(ISYM)
  NSB = (NAS*(NAS+1))/2
  if (NSB > 0) call mma_allocate(SB,NSB,Label='SB')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      ISADR = (ITU*(ITU-1))/2+IXY
      IXT = IXABS+NASHT*(ITABS-1)
      IYU = IYABS+NASHT*(IUABS-1)
      IP1 = max(IXT,IYU)
      IP2 = min(IXT,IYU)
      IP = (IP1*(IP1-1))/2+IP2
      value = Four*PREF(IP)
      ! Add  -4 dxt Dyu + 8dxt dyu
      if (IXABS == ITABS) then
        ID1 = max(IYABS,IUABS)
        ID2 = min(IYABS,IUABS)
        ID = (ID1*(ID1-1))/2+ID2
        value = value-Four*DREF(ID)
        if (IYABS == IUABS) value = value+Eight
      end if
      ! Add  -4 dyu Dxt
      if (IYABS == IUABS) then
        ID1 = max(IXABS,ITABS)
        ID2 = min(IXABS,ITABS)
        ID = (ID1*(ID1-1))/2+ID2
        value = value-Four*DREF(ID)
      end if
      ! Add  +2 dyt Dxu
      if (IYABS == ITABS) then
        ID1 = max(IXABS,IUABS)
        ID2 = min(IXABS,IUABS)
        ID = (ID1*(ID1-1))/2+ID2
        value = value+Two*DREF(ID)
      end if
      ! Add  -4dxu dyt + 2dxu Dyt
      if (IXABS == IUABS) then
        ID1 = max(IYABS,ITABS)
        ID2 = min(IYABS,ITABS)
        ID = (ID1*(ID1-1))/2+ID2
        value = value+Two*DREF(ID)
        if (IYABS == ITABS) value = value-Four
      end if
      ISADR = (ITU*(ITU-1))/2+IXY
      SB(ISADR) = value
    end do
  end do
  NASP = NTGEU(ISYM)
  NSBP = (NASP*(NASP+1))/2
  if (NSBP > 0) call mma_allocate(SBP,NSBP,Label='SBP')
  NASM = NTGTU(ISYM)
  NSBM = (NASM*(NASM+1))/2
  if (NSBM > 0) call mma_allocate(SBM,NSBM,Label='SBM')
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
        ISADR = (ITU*(ITU-1))/2+IXY
      else
        ISADR = (IXY*(IXY-1))/2+ITU
      end if
      STUXY = SB(ISADR)
      if (ITU >= IYX) then
        ISADR = (ITU*(ITU-1))/2+IYX
      else
        ISADR = (IYX*(IYX-1))/2+ITU
      end if
      STUYX = SB(ISADR)
      ISPADR = (ITGEU*(ITGEU-1))/2+IXGEY
      SBP(ISPADR) = STUXY+STUYX
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      ISMADR = (ITGTU*(ITGTU-1))/2+IXGTY
      SBM(ISMADR) = STUXY-STUYX
    end do
  end do
  if (NSB > 0) call mma_deallocate(SB)

  ! Write to disk, and save size and address.
  if (NSBP > 0) then
    IDISK = IDSMAT(ISYM,2)
    call DDAFILE(LUSBT,1,SBP,NSBP,IDISK)
    call mma_deallocate(SBP)
  end if
  if (NSBM > 0) then
    if (NINDEP(ISYM,3) > 0) then
      IDISK = IDSMAT(ISYM,3)
      call DDAFILE(LUSBT,1,SBM,NSBM,IDISK)
    end if
    call mma_deallocate(SBM)
  end if
end do

end subroutine MKSB
