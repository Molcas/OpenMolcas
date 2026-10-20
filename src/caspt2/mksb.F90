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

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: KTGTU, KTU, MTGEU, MTU
use EQSOLV, only: IDSMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: NASHT, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
use stdalloc, only: mma_allocate, mma_deallocate
use constants, only: Two, Four, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
integer(kind=iwp) :: ID, IDISK, IP, ISADR, ISMADR, ISPADR, ISYM, ITABS, ITGEU, ITGEUABS, ITGTU, ITU, ITUABS, IUABS, IXABS, IXGEY, &
                     IXGEYABS, IXGTY, IXT, IXY, IXYABS, IYABS, IYU, IYX, NAS, NASM, NASP, NINP, NSB, NSBM, NSBP
real(kind=wp) :: STUXY, STUYX, Val
real(kind=wp), allocatable :: SB(:), SBM(:), SBP(:)

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
  NSB = nTri_Elem(NAS)
  if (NSB > 0) call mma_allocate(SB,NSB,Label='SB')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      ISADR = iTri(ITU,IXY)
      IXT = IXABS+NASHT*(ITABS-1)
      IYU = IYABS+NASHT*(IUABS-1)
      IP = iTri(IXT,IYU)
      Val = Four*PREF(IP)
      ! Add  -4 dxt Dyu + 8dxt dyu
      if (IXABS == ITABS) then
        ID = iTri(IYABS,IUABS)
        Val = Val-Four*DREF(ID)
        if (IYABS == IUABS) Val = Val+Eight
      end if
      ! Add  -4 dyu Dxt
      if (IYABS == IUABS) then
        ID = iTri(IXABS,ITABS)
        Val = Val-Four*DREF(ID)
      end if
      ! Add  +2 dyt Dxu
      if (IYABS == ITABS) then
        ID = iTri(IXABS,IUABS)
        Val = Val+Two*DREF(ID)
      end if
      ! Add  -4dxu dyt + 2dxu Dyt
      if (IXABS == IUABS) then
        ID = iTri(IYABS,ITABS)
        Val = Val+Two*DREF(ID)
        if (IYABS == ITABS) Val = Val-Four
      end if
      ISADR = iTri(ITU,IXY)
      SB(ISADR) = Val
    end do
  end do
  NASP = NTGEU(ISYM)
  NSBP = nTri_Elem(NASP)
  if (NSBP > 0) call mma_allocate(SBP,NSBP,Label='SBP')
  NASM = NTGTU(ISYM)
  NSBM = nTri_Elem(NASM)
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
      ISADR = iTri(ITU,IXY)
      STUXY = SB(ISADR)
      ISADR = iTri(ITU,IYX)
      STUYX = SB(ISADR)
      ISPADR = iTri(ITGEU,IXGEY)
      SBP(ISPADR) = STUXY+STUYX
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      ISMADR = iTri(ITGTU,IXGTY)
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
