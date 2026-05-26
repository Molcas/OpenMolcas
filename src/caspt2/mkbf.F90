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

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: KTGTU, KTU, MTGEU, MTU
use EQSOLV, only: IDBMAT, IDSMAT
use caspt2_global, only: ipea_shift, LUSBT
use caspt2_module, only: EASUM, NASHT, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF), FP(NPREF)
integer(kind=iwp) :: I, IBADR, IBMADR, IBPADR, IDIAG, IDISK, IDSM, IDSP, IDT, IDU, INSM, IP, ISYM, ITABS, ITGEU, ITGEUABS, ITGTU, &
                     ITU, ITUABS, ITX, IUABS, IUY, IXABS, IXGEY, IXGEYABS, IXGTY, IXY, IXYABS, IYABS, IYX, NAS, NASM, NASP, NBF, &
                     NBFM, NBFP, NINP
real(kind=wp) :: BTUXY, BTUYX
real(kind=wp), allocatable :: BF(:), BFM(:), BFP(:), SDM(:), SDP(:), SM(:), SP(:)

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
  NBF = nTri_Elem(NAS)
  if (NBF > 0) call mma_allocate(BF,NBF,LABEL='BF')
  do ITU=1,NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      IBADR = iTri(ITU,IXY)
      ITX = ITABS+NASHT*(IXABS-1)
      IUY = IUABS+NASHT*(IYABS-1)
      IP = iTri(ITX,IUY)
      BF(IBADR) = Four*(FP(IP)-EASUM*PREF(IP))
    end do
  end do
  NASP = NTGEU(ISYM)
  NBFP = nTri_Elem(NASP)
  if (NBFP > 0) then
    call mma_allocate(BFP,NBFP,Label='BFP')
    !GG.Nov03  Load in SDP the diagonal elements of SFP matrix:
    call mma_allocate(SP,NBFP,Label='SP')
    call mma_allocate(SDP,NASP,Label='SDP')
    IDSP = IDSMAT(ISYM,8)
    call DDAFILE(LUSBT,2,SP,NBFP,IDSP)
    IDIAG = 0
    do I=1,NASP
      IDIAG = IDIAG+I
      SDP(I) = SP(IDIAG)
    end do
    call mma_deallocate(SP)
    !GG End
  end if
  NASM = NTGTU(ISYM)
  NBFM = nTri_Elem(NASM)
  if (NBFM > 0) then
    call mma_allocate(BFM,NBFM,Label='BFM')
    !GG.Nov03  Load in SDM the diagonal elements of SFM matrix:
    call mma_allocate(SM,NBFM,Label='SM')
    call mma_allocate(SDM,NASM,Label='SDM')
    IDSM = IDSMAT(ISYM,9)
    call DDAFILE(LUSBT,2,SM,NBFM,IDSM)
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
      IBADR = iTri(ITU,IXY)
      BTUXY = BF(IBADR)
      IBADR = iTri(ITU,IYX)
      BTUYX = BF(IBADR)
      IBPADR = iTri(ITGEU,IXGEY)
      BFP(IBPADR) = BTUXY+BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = nTri_Elem(ITABS)
        IDU = nTri_Elem(IUABS)
        BFP(IBPADR) = BFP(IBPADR)+ipea_shift*Half*(Four-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
      end if
      !GG End
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      IBMADR = iTri(ITGTU,IXGTY)
      BFM(IBMADR) = BTUXY-BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = nTri_Elem(ITABS)
        IDU = nTri_Elem(IUABS)
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
