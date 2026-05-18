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

subroutine MKSF(PREF,NPREF)
! Set up the matrices SFP(tu,xy) and SFM(tu,xy)
! Formulae used:
!    SF(tu,xy)= 4 Ptxuy
!    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
!    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)

use definitions, only: iwp, wp
use constants, only: Four
use SUPERINDEX, only: MTU, MTGEU, KTU, MTGEU, KTGTU
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NTU, NTUES, NASHT, NTGEU, NTGEUES, NTGTU, NTGTUES

implicit none
integer(kind=iwp), intent(in) :: NPREF
real(kind=wp), intent(in) :: PREF(NPREF)
real(kind=wp), allocatable :: SF(:), SFP(:), SFM(:)
integer(kind=iwp) ISYM, NINP, NAS, NSF, ITU, ITUABS, ITABS, IUABS, IXY, IXYABS, IXABS, IYABS, ISADR, ITX, IUY, IP1, IP2, IP, &
                  IDISK, ISMADR, ISPADR, ITGEU, ITGEUABS, ITGTU, IXGEY, IXGEYABS, IXGTY, IYX, NASM, NASP, NSFM, NSFP
real(kind=wp) value, STUXY, STUYX

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NINP = NINDEP(ISYM,8)
  if (NINP == 0) cycle
  NAS = NTU(ISYM)
  NSF = (NAS*(NAS+1))/2
  if (NSF > 0) call mma_allocate(SF,NSF,Label='SF')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      ISADR = (ITU*(ITU-1))/2+IXY
      ITX = ITABS+NASHT*(IXABS-1)
      IUY = IUABS+NASHT*(IYABS-1)
      IP1 = max(ITX,IUY)
      IP2 = min(ITX,IUY)
      IP = (IP1*(IP1-1))/2+IP2
      value = Four*PREF(IP)
      SF(ISADR) = value
    end do
  end do
  NASP = NTGEU(ISYM)
  NSFP = (NASP*(NASP+1))/2
  if (NSFP > 0) call mma_allocate(SFP,NSFP,Label='SFP')
  NASM = NTGTU(ISYM)
  NSFM = (NASM*(NASM+1))/2
  if (NSFM > 0) call mma_allocate(SFM,NSFM,Label='SFM')
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
      STUXY = SF(ISADR)
      if (ITU >= IYX) then
        ISADR = (ITU*(ITU-1))/2+IYX
      else
        ISADR = (IYX*(IYX-1))/2+ITU
      end if
      STUYX = SF(ISADR)
      ISPADR = (ITGEU*(ITGEU-1))/2+IXGEY
      SFP(ISPADR) = STUXY+STUYX
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      ISMADR = (ITGTU*(ITGTU-1))/2+IXGTY
      SFM(ISMADR) = STUXY-STUYX
    end do
  end do
  if (NSF > 0) call mma_deallocate(SF)

  ! Write to disk
  if ((NSFP > 0) .and. (NINDEP(ISYM,8) > 0)) then
    IDISK = IDSMAT(ISYM,8)
    call DDAFILE(LUSBT,1,SFP,NSFP,IDISK)
    call mma_deallocate(SFP)
  end if
  if (NSFM > 0) then
    if (NINDEP(ISYM,9) > 0) then
      IDISK = IDSMAT(ISYM,9)
      call DDAFILE(LUSBT,1,SFM,NSFM,IDISK)
    end if
    call mma_deallocate(SFM)
  end if
end do

end subroutine MKSF
