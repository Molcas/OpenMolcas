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

subroutine MKBNEVF(nAshT,NG3,Hbar,Gact,G2,G3,idxG3)

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: LUSBT
use caspt2_module, only: IASYM, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
use EQSOLV, only: IDBMAT
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
real(kind=wp), intent(in) :: Hbar(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: IBADR, IBMADR, IBPADR, ID, IDISK, iG3, IST, ISU, ISV, ISX, ISY, ISYM, ISZ, IT, ITABS, ITGEU, ITGEUABS, ITGTU, &
                     ITU, ITUABS, ITUVS, IU, IUABS, IV, IVABS, IVX, IW, IWABS, IX, IXABS, IXGEY, IXGEYABS, IXGTY, IXY, IXYABS, &
                     IXYZS, IY, IYABS, IYX, IYZ, IZ, IZABS, NAS, NASM, NASP, NBF, NBFM, NBFP
real(kind=wp) :: BTUXY, BTUYX, G3VAL, Val
real(kind=wp), allocatable :: BF(:), BFM(:), BFP(:)

do ISYM=1,NSYM
  if (NINDEP(ISYM,8) == 0) cycle
  NAS = NTU(ISYM)
  NBF = nTri_Elem(NAS)
  if (NBF > 0) call mma_Allocate(BF,NBF,LABEL='BF')
  BF(:) = Zero
  call Init_Tsk(ID,NAS)

  ! K(tu,xy) = - h(xv)D(tu,vy) - h(yv)D(tu,xv)
  ! K(tu,xy) = - 0.5*(vw|xz) ( del(vz)D(tu,wy) + del(yv)D(tu,zw) )
  !            - 0.5*(vw|yz) ( del(xv)D(tu,wz) + del(vz)D(tu,xw) )
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
      !ITX = ITABS+NASHT*(IXABS-1)
      !IUY = IUABS+NASHT*(IYABS-1)
      !IP = iTri(ITX,IUY)
      Val = Zero
      ! Hbar term
      do IV=1,NASHT
        IVABS = IV!+NAES(ISYM)
        Val = Val-Hbar(IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
        Val = Val-Hbar(IYABS,IVABS)*G2(ITABS,IXABS,IUABS,IVABS)
      end do
      !do IV=1,nAsh(iSymX)
      !  IVABS = IV+NAES(iSymX)
      !  Val = Val-Hbar(IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
      !end do
      !do IV=1,nAsh(iSymY)
      !  IVABS = IV+NAES(iSymY)
      !  Val = Val-Hbar(IYABS,IVABS)*G2(ITABS,IXABS,IUABS,IVABS)
      !end do
      ! RDM3 term
      ! RDM2 term
      do IV=1,NASHT
        IVABS = IV!+NAES(ISYM)
        do IW=1,NASHT
          IWABS = IW!+NAES(ISYM)
          do IZ=1,NASHT
            IZABS = IZ!+NAES(ISYM)
            if (IVABS == IZABS) Val = Val-Half*Gact(IVABS,IWABS,IXABS,IZABS)*G2(ITABS,IWABS,IUABS,IYABS)
            if (IYABS == IVABS) Val = Val-Half*Gact(IVABS,IWABS,IXABS,IZABS)*G2(ITABS,IZABS,IUABS,IWABS)
            if (IXABS == IVABS) Val = Val-Half*Gact(IVABS,IWABS,IYABS,IZABS)*G2(ITABS,IWABS,IUABS,IZABS)
            if (IVABS == IZABS) Val = Val-Half*Gact(IVABS,IWABS,IYABS,IZABS)*G2(ITABS,IXABS,IUABS,IWABS)
          end do
        end do
      end do
      BF(IBADR) = Val
      !write(u6,'(4i3,f20.10)') itabs,iuabs,ixabs,iyabs,Val
    end do
  end do
  call Free_Tsk(ID)

  ! The relation with the "standard" RDM and "MOLCAS" RDM is:
  ! D(MOLCAS)(tuvxyz) = D(standard)(tvy,uxz)
  ! D(standard)(tuvxyz) = D(MOLCAS)(txu,yvz)
  ! K(tu,xy) = - (vw|xz) D(tuv,zyw)
  !            - (vw|yz) D(tuv,xzw)
  ! DK(db,a'b') = RDM3(a'b'c',abc)*G(d,a,c',c)
  ! DK(wy,t u ) = RDM3(tuv,xyz)*G(w,x,v,z)
  !>DK(wv,t x ) = RDM3(txu,yvz)*G(w,y,u,z)

  ! DK(ad,a'b') = RDM3(a'b'c',abc)*G(d,b,c',c)
  ! DK(xw,t u ) = RDM3(tuv,xyz)*G(w,y,v,z)
  !>DK(yw,t x ) = RDM3(txu,yvz)*G(w,v,u,z)
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

    G3VAL = G3(iG3)

    ! F(tuvxyz) -> BC(vut,xyz)
    if (MUL(iST,iSV) == iSym) call Add_MKBNEVF(iT,iU,iV,iX,iY,iZ)

    if ((iTU /= iVX) .or. (iVX /= iYZ)) then
      if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
        ! F(vxtuyz) -> BC(txv,uyz)
        if (MUL(iSV,iST) == iSym) call Add_MKBNEVF(iV,iX,iT,iU,iY,iZ)
        ! F(yzvxtu) -> BC(vzy,xtu)
        if (MUL(iSY,iSV) == iSym) call Add_MKBNEVF(iY,iZ,iV,iX,iT,iU)
        ! F(tuyzvx) -> BC(yut,zvx)
        if (MUL(iST,iSY) == iSym) call Add_MKBNEVF(iT,iU,iY,iZ,iV,iX)
      end if
      ! F(yztuvx) -> BC(tzy,uvx)
      if (MUL(iSY,iST) == iSym) call Add_MKBNEVF(iY,iZ,iT,iU,iV,iX)
      ! F(vxyztu) -> BC(yxv,ztu)
      if (MUL(iSV,iSY) == iSym) call Add_MKBNEVF(iV,iX,iY,iZ,iT,iU)
    end if

    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle

    ! F(utxvzy) -> BC(xtu,vzy)
    if (MUL(iSU,iSX) == iSym) call Add_MKBNEVF(iU,iT,iX,iV,iZ,iY)

    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
      ! F(xvutzy) -> BC(uvx,tzy)
      if (MUL(iSX,iSU) == iSym) call Add_MKBNEVF(iX,iV,iU,iT,iZ,iY)
      ! F(zyxvut) -> BC(xyz,vut)
      if (MUL(iSZ,iSX) == iSym) call Add_MKBNEVF(iZ,iY,iX,iV,iU,iT)
      ! F(utzyxv) -> BC(ztu,yxv)
      if (MUL(iSU,iSZ) == iSym) call Add_MKBNEVF(iU,iT,iZ,iY,iX,iV)
    end if
    ! F(zyutxv) -> BC(uyz,txv)
    if (MUL(iSZ,iSU) == iSym) call Add_MKBNEVF(iZ,iY,iU,iT,iX,iV)
    ! F(xvzyut) -> BC(zvx,yut)
    if (MUL(iSX,iSZ) == iSym) call Add_MKBNEVF(iX,iV,iZ,iY,iU,iT)
  end do

  BF(:) = Two*BF(:)

  NASP = NTGEU(ISYM)
  NBFP = nTri_Elem(NASP)
  if (NBFP > 0) call mma_allocate(BFP,NBFP,Label='BFP')
  NASM = NTGTU(ISYM)
  NBFM = nTri_Elem(NASM)
  if (NBFM > 0) call mma_allocate(BFM,NBFM,Label='BFM')

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
      if (ITABS == IUABS) cycle
      if (IXABS == IYABS) cycle
      ITGTU = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
      IXGTY = KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
      IBMADR = iTri(ITGTU,IXGTY)
      BFM(IBMADR) = BTUXY-BTUYX
    end do
  end do

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    call GADGOP(BFP,NBFP,'+')
    call GADGOP(BFM,NBFM,'+')
  end if
# endif

  if (NBF > 0) call mma_deallocate(BF)

  if ((NBFP > 0) .and. (NINDEP(ISYM,8) > 0)) then
    IDISK = IDBMAT(ISYM,8)
    call DDAFILE(LUSBT,1,BFP,NBFP,IDISK)
    call mma_deallocate(BFP)
  end if
  if (NBFM > 0) then
    if (NINDEP(ISYM,9) > 0) then
      IDISK = IDBMAT(ISYM,9)
      call DDAFILE(LUSBT,1,BFM,NBFM,IDISK)
    end if
    call mma_deallocate(BFM)
  end if
end do

contains

subroutine Add_MKBNEVF(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_)

  integer(kind=iwp), intent(in) :: ITABS_, IXABS_, IUABS_, IYABS_, IVABS_, IZABS_
  integer(kind=iwp) :: IBADR, ITU_, IWABS, IXY_

  iTU_ = KTU(ITABS_,IUABS_)-NTUES(ISYM)
  do IWABS=1,nAshT
    if (MUL(IASYM(IWABS),IASYM(IYABS_)) == iSym) then
      iXY_ = KTU(IWABS,IYABS_)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTri(iTU_,iXY_)
        BF(IBADR) = BF(IBADR)-G3VAL*Gact(IWABS,IXABS_,IVABS_,IZABS_)
      end if
    end if
    if (MUL(IASYM(IXABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IXABS_,IWABS)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTri(iTU_,iXY_)
        BF(IBADR) = BF(IBADR)-G3VAL*Gact(IWABS,IYABS_,IVABS_,IZABS_)
      end if
    end if
  end do

end subroutine Add_MKBNEVF

end subroutine MKBNEVF
