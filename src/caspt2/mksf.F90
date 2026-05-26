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

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: KTGTU, KTU, MTGEU, MTGEU, MTU
use EQSOLV, only: IDSMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: NASHT, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NPREF
real(kind=wp), intent(in) :: PREF(NPREF)
integer(kind=iwp) :: IDISK, IP, ISADR, ISMADR, ISPADR, ISYM, ITABS, ITGEU, ITGEUABS, ITGTU, ITU, ITUABS, ITX, IUABS, IUY, IXABS, &
                     IXGEY, IXGEYABS, IXGTY, IXY, IXYABS, IYABS, IYX, NAS, NASM, NASP, NINP, NSF, NSFM, NSFP
real(kind=wp) :: STUXY, STUYX
real(kind=wp), allocatable :: SF(:), SFM(:), SFP(:)

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NINP = NINDEP(ISYM,8)
  if (NINP == 0) cycle
  NAS = NTU(ISYM)
  NSF = nTri_Elem(NAS)
  if (NSF > 0) call mma_allocate(SF,NSF,Label='SF')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      ISADR = iTri(ITU,IXY)
      ITX = ITABS+NASHT*(IXABS-1)
      IUY = IUABS+NASHT*(IYABS-1)
      IP = iTri(ITX,IUY)
      SF(ISADR) = Four*PREF(IP)
    end do
  end do
  NASP = NTGEU(ISYM)
  NSFP = nTri_Elem(NASP)
  if (NSFP > 0) call mma_allocate(SFP,NSFP,Label='SFP')
  NASM = NTGTU(ISYM)
  NSFM = nTri_Elem(NASM)
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
      ISADR = iTri(ITU,IXY)
      STUXY = SF(ISADR)
      ISADR = iTri(ITU,IYX)
      STUYX = SF(ISADR)
      ISPADR = iTri(ITGEU,IXGEY)
      SFP(ISPADR) = STUXY+STUYX
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      ISMADR = iTri(ITGTU,IXGTY)
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
