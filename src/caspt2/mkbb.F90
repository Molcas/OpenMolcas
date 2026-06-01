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

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: KTGTU, KTU, MTGEU, MTU
use EQSOLV, only: IDBMAT, IDSMAT
use caspt2_global, only: ipea_shift, LUSBT
use caspt2_module, only: EASUM, EPSA, NASHT, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Four, Eight, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF), FD(NDREF), FP(NPREF)
integer(kind=iwp) :: I, IBADR, IBMADR, IBPADR, ID, IDIAG, IDISK, IDSM, IDSP, IDT, IDU, INSM, IP, ISYM, ITABS, ITGEU, ITGEUABS, &
                     ITGTU, ITU, ITUABS, IUABS, IXABS, IXGEY, IXGEYABS, IXGTY, IXT, IXY, IXYABS, IYABS, IYU, IYX, NAS, NASM, NASP, &
                     NBB, NBBM, NBBP, NINP
real(kind=wp) :: ATUX, ATUXY, ATUY, ATYU, ATYX, BTUXY, BTUYX, ET, EU, EX, EY, Val
real(kind=wp), allocatable :: BB(:), BBP(:), SP(:), SDP(:), BBM(:), SM(:), SDM(:)

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
  NBB = nTri_Elem(NAS)
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
      IBADR = iTri(ITU,IXY)
      IXT = IXABS+NASHT*(ITABS-1)
      IYU = IYABS+NASHT*(IUABS-1)
      IP = iTri(IXT,IYU)
      ATUXY = EASUM-ET-EU-EX-EY
      Val = Four*(FP(IP)-ATUXY*PREF(IP))
      ! Add  + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
      if (IXABS == ITABS) then
        ID = iTri(IYABS,IUABS)
        ATYU = EASUM-ET-EY-EU
        Val = Val+Four*(ATYU*DREF(ID)-FD(ID))
        ! Add  + 8*dxt*dyu (Et+Ey)
        if (IYABS == IUABS) Val = Val+Eight*(ET+EY)
      end if
      ! Add  + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
      if (IYABS == IUABS) then
        ID = iTri(IXABS,ITABS)
        ATYX = EASUM-ET-EY-EX
        Val = Val+Four*(ATYX*DREF(ID)-FD(ID))
      end if
      ! Add  - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
      if (IYABS == ITABS) then
        ID = iTri(IXABS,IUABS)
        ATUX = EASUM-ET-EU-EX
        Val = Val-Two*(ATUX*DREF(ID)-FD(ID))
        ! Add  - 4*dxu*dyt (Et+Ex)
        if (IXABS == IUABS) Val = Val-Four*(ET+EX)
      end if
      ! Add  - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
      if (IXABS == IUABS) then
        ID = iTri(IYABS,ITABS)
        ATUY = EASUM-ET-EU-EY
        Val = Val-Two*(ATUY*DREF(ID)-FD(ID))
      end if
      BB(IBADR) = Val
    end do
  end do
  NASP = NTGEU(ISYM)
  NBBP = nTri_Elem(NASP)
  if (NBBP > 0) then
    call mma_allocate(BBP,NBBP,Label='BBP')
    !GG.Nov03  Load in SDP the diagonal elements of SBP matrix:
    call mma_allocate(SP,NBBP,Label='SP')
    call mma_allocate(SDP,NASP,Label='SDP')
    IDSP = IDSMAT(ISYM,2)
    call DDAFILE(LUSBT,2,SP,NBBP,IDSP)
    IDIAG = 0
    do I=1,NASP
      IDIAG = IDIAG+I
      SDP(I) = SP(IDIAG)
    end do
    call mma_deallocate(SP)
    !GG End
  end if

  NASM = NTGTU(ISYM)
  NBBM = nTri_Elem(NASM)
  if (NBBM > 0) then
    call mma_allocate(BBM,NBBM,Label='BBM')
    !GG.Nov03  Load in SDM the diagonal elements of SBM matrix:
    call mma_allocate(SDM,NASM,Label='SDM')
    if (NINDEP(ISYM,3) > 0) then
      call mma_allocate(SM,NBBM,Label='SM')
      IDSM = IDSMAT(ISYM,3)
      call DDAFILE(LUSBT,2,SM,NBBM,IDSM)
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
      IBADR = iTri(ITU,IXY)
      BTUXY = BB(IBADR)
      IBADR = iTri(ITU,IYX)
      BTUYX = BB(IBADR)
      IBPADR = iTri(ITGEU,IXGEY)
      BBP(IBPADR) = BTUXY+BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = nTri_Elem(ITABS)
        IDU = nTri_Elem(IUABS)
        BBP(IBPADR) = BBP(IBPADR)+ipea_shift*half*(DREF(IDT)+DREF(IDU))*SDP(ITGEU)
      end if
      !GG End
      if (NINDEP(ISYM,3) < 1) cycle
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      IBMADR = iTri(ITGTU,IXGTY)
      BBM(IBMADR) = BTUXY-BTUYX
      !GG.Nov03
      if (ITGEU == IXGEY) then
        IDT = nTri_Elem(ITABS)
        IDU = nTri_Elem(IUABS)
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
