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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine MKBNEVB(nAshT,NG3,Hbar,Gact,G1,G2,G3,idxG3)

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: LUSBT
use caspt2_module, only: IASYM, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
use EQSOLV, only: IDBMAT
use NEVPT2_mod, only: tG2, tG3
use SUPERINDEX, only: KTGTU, KTU, MTGEU, MTU
use Symmetry_Info, only: Mul
use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: nAshT, NG3
real(kind=wp), intent(in) :: Hbar(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(*)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: IBADR, IBMADR, IBPADR, ID, IDISK, iG3, IST, ISU, ISV, ISX, ISY, ISYM, ISZ, IT, ITABS, ITGEU, ITGEUABS, ITGTU, &
                     ITU, ITUABS, ITUVS, IU, IUABS, IV, IVABS, IVX, IW, IWABS, IX, IXABS, IXGEY, IXGEYABS, IXGTY, IXY, IXYABS, &
                     IXYZS, IY, IYABS, IYX, IYZ, IZ, IZABS, NAS, NASM, NASP, NBB, NBBM, NBBP
real(kind=wp) :: BTUXY, BTUYX, G3VAL, Val
real(kind=wp), allocatable :: BB(:), BBM(:), BBP(:)

do ISYM=1,NSYM
  if (NINDEP(ISYM,2) == 0) cycle
  NAS = NTU(ISYM)
  NBB = nTri_Elem(NAS)
  if (NBB > 0) call mma_Allocate(BB,NBB,LABEL='BB')
  BB(:) = Zero
  call Init_Tsk(ID,NAS)

  ! tilde{R}^{(2)} --> tG2
  ! The abcd <--> acbd symmetry is considered in the tG2 subroutine
  !do ITU=1,NAS
  !  if (.not. Rsv_Tsk(ID,ITU)) cycle
  do while (Rsv_Tsk(ID,ITU))
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    !iSymT = IASYM(ITABS)
    !iSymU = IASYM(IUABS)
    do IXY=1,ITU
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      !iSymX = IASYM(IXABS)
      !iSymY = IASYM(IYABS)
      IBADR = iTri(ITU,IXY)
      Val = Zero
      ! Hbar term
      do IV=1,NASHT
        IVABS = IV!+NAES(ISYM)
        Val = Val+Hbar(IVABS,IXABS)*tG2(ITABS,IUABS,IVABS,IYABS,G1,G2)
        Val = Val+Hbar(IVABS,IYABS)*tG2(ITABS,IUABS,IXABS,IVABS,G1,G2)
      end do
      ! RDM2 term
      do IV=1,NASHT
        IVABS = IV!+NAES(ISYM)
        do IW=1,NASHT
          IWABS = IW!+NAES(ISYM)
          do IZ=1,NASHT
            IZABS = IZ!+NAES(ISYM)
            if (IVABS == IZABS) Val = Val+Two*Gact(IVABS,IZABS,IWABS,IXABS)*tG2(ITABS,IUABS,IWABS,IYABS,G1,G2)
            if (IWABS == IZABS) Val = Val-Half*Gact(IVABS,IZABS,IWABS,IXABS)*tG2(ITABS,IUABS,IVABS,IYABS,G1,G2)
            if (IYABS == IZABS) Val = Val-Half*Gact(IVABS,IZABS,IWABS,IXABS)*tG2(ITABS,IUABS,IWABS,IVABS,G1,G2)

            if (IVABS == IZABS) Val = Val+Two*Gact(IVABS,IZABS,IWABS,IYABS)*tG2(ITABS,IUABS,IXABS,IWABS,G1,G2)
            if (IXABS == IZABS) Val = Val-Half*Gact(IVABS,IZABS,IWABS,IYABS)*tG2(ITABS,IUABS,IVABS,IWABS,G1,G2)
            if (IWABS == IZABS) Val = Val-Half*Gact(IVABS,IZABS,IWABS,IYABS)*tG2(ITABS,IUABS,IXABS,IVABS,G1,G2)
          end do
        end do
      end do
      BB(IBADR) = Val
      !write(u6,'(4i3,f20.10)') itabs,iuabs,ixabs,iyabs,Val
    end do
  end do
  call Free_Tsk(ID)

  ! K (tu,xy) = -(vz|wx)*D(tuz,wyv) - (vz|wy)*D(tuz,xwv)
  ! K1(tu,wy) = -(zv|xw)*D(tuv,xyz)
  !>K1(tx,wv) = -(zu|yw)*D(txu,yvz)
  ! K (tu,xy) = -(vz|wy)*D(tuz,xwv)
  ! K2(tu,xw) = -(zv|yw)*D(tuv,xyz)
  !>K2(tx,yw) = -(zu|vw)*D(txu,yvz)
  do IG3=1,NG3
    iT = idxG3(1,iG3)
    iU = idxG3(2,iG3)
    iV = idxG3(3,iG3)
    iX = idxG3(4,iG3)
    iY = idxG3(5,iG3)
    iZ = idxG3(6,iG3)
    iST = IASYM(iT)
    iSU = IASYM(iU)
    iSV = IASYM(iV)
    iSX = IASYM(iX)
    iSY = IASYM(iY)
    iSZ = IASYM(iZ)
    ituvs = MUL(IST,MUL(ISU,ISV))
    ixyzs = MUL(ISX,MUL(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU = iT+NASHT*(iU-1)
    iVX = iV+NASHT*(iX-1)
    iYZ = iY+NASHT*(iZ-1)

    ! G3(tuvxyz) = <Psi|t+ v+ y+ z x u|Psi>
    G3VAL = tG3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))

    ! D(tuvxyz)
    if (MUL(iST,iSV) == iSym) call Add_MKBNEVB(iT,iU,iV,iX,iY,iZ)

    if ((iTU /= iVX) .or. (iVX /= iYZ)) then
      if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
        ! D(vxtuyz)
        if (MUL(iSV,iST) == iSym) call Add_MKBNEVB(iV,iX,iT,iU,iY,iZ)
        ! D(yzvxtu)
        if (MUL(iSY,iSV) == iSym) call Add_MKBNEVB(iY,iZ,iV,iX,iT,iU)
        ! D(tuyzvx)
        if (MUL(iST,iSY) == iSym) call Add_MKBNEVB(iT,iU,iY,iZ,iV,iX)
      end if
      ! D(yztuvx)
      if (MUL(iSY,iST) == iSym) call Add_MKBNEVB(iY,iZ,iT,iU,iV,iX)
      ! D(vxyztu)
      if (MUL(iSV,iSY) == iSym) call Add_MKBNEVB(iV,iX,iY,iZ,iT,iU)
    end if

    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle

    ! D(utxvzy)
    if (MUL(iSU,iSX) == iSym) call Add_MKBNEVB(iU,iT,iX,iV,iZ,iY)

    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
      ! D(xvutzy)
      if (MUL(iSX,iSU) == iSym) call Add_MKBNEVB(iX,iV,iU,iT,iZ,iY)
      ! D(zyxvut)
      if (MUL(iSZ,iSX) == iSym) call Add_MKBNEVB(iZ,iY,iX,iV,iU,iT)
      ! D(utzyxv)
      if (MUL(iSU,iSZ) == iSym) call Add_MKBNEVB(iU,iT,iZ,iY,iX,iV)
    end if
    ! D(zyutxv)
    if (MUL(iSZ,iSU) == iSym) call Add_MKBNEVB(iZ,iY,iU,iT,iX,iV)
    ! D(xvzyut)
    if (MUL(iSX,iSZ) == iSym) call Add_MKBNEVB(iX,iV,iZ,iY,iU,iT)
  end do

  BB(:) = Two*BB(:)

  NASP = NTGEU(ISYM)
  NBBP = nTri_Elem(NASP)
  if (NBBP > 0) call mma_allocate(BBP,NBBP,Label='BBP')
  NASM = NTGTU(ISYM)
  NBBM = nTri_Elem(NASM)
  if (NBBM > 0) call mma_allocate(BBM,NBBM,Label='BBM')

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
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      IBMADR = iTri(ITGTU,IXGTY)
      BBM(IBMADR) = BTUXY-BTUYX
    end do
  end do
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    call GADGOP(BBP,NBBP,'+')
    call GADGOP(BBM,NBBM,'+')
  end if
# endif

  if (NBB > 0) call mma_deallocate(BB)

  if ((NBBP > 0) .and. (NINDEP(ISYM,2) > 0)) then
    IDISK = IDBMAT(ISYM,2)
    call DDAFILE(LUSBT,1,BBP,NBBP,IDISK)
    call mma_deallocate(BBP)
  end if
  if (NBBM > 0) then
    if (NINDEP(ISYM,3) > 0) then
      IDISK = IDBMAT(ISYM,3)
      call DDAFILE(LUSBT,1,BBM,NBBM,IDISK)
    end if
    call mma_deallocate(BBM)
  end if
end do

contains

subroutine Add_MKBNEVB(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_)

  integer(kind=iwp), intent(in) :: ITABS_, IXABS_, IUABS_, IYABS_, IVABS_, IZABS_
  integer(kind=iwp) :: IBADR, ITU_, IWABS, IXY_

  iTU_ = KTU(ITABS_,IUABS_)-NTUES(ISYM)
  do IWABS=1,nAshT
    if (MUL(IASYM(IWABS),IASYM(IYABS_)) == iSym) then
      iXY_ = KTU(IWABS,IYABS_)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTri(iTU_,iXY_)
        BB(IBADR) = BB(IBADR)-G3VAL*Gact(IZABS_,IVABS_,IXABS_,IWABS)
      end if
    end if
    if (MUL(IASYM(IXABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IXABS_,IWABS)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTri(iTU_,iXY_)
        BB(IBADR) = BB(IBADR)-G3VAL*Gact(IZABS_,IVABS_,IYABS_,IWABS)
      end if
    end if
  end do

end subroutine Add_MKBNEVB

end subroutine MKBNEVB
