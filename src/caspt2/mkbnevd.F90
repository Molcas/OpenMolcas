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

subroutine MKBNEVD(nAshT,NG3,Hbar,Gact,G1,G2,G3,idxG3)

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: LUSBT
use caspt2_module, only: IASYM, NINDEP, NSYM, NTU, NTUES
use EQSOLV, only: IDBMAT
use SUPERINDEX, only: KTU, MTU
use nevpt2_mod, only: ex2, ex3
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
integer(kind=iwp) :: IB11, IB12, IB21, IB22, id, IDISK, iG3, IST, ISU, ISV, ISX, ISY, ISYM, ISZ, IT, ITABS, ITU, ITU2, ITUABS, &
                     ITUVS, IU, IUABS, IV, IVABS, IVX, IW, IWABS, IX, IXABS, IXY, IXY2, IXYABS, IXYZS, IY, IYABS, IYZ, IZ, IZABS, &
                     NAS, NBD
real(kind=wp) :: B11, B22, G3VAL, tmp
real(kind=wp), allocatable :: BD(:), BD1(:), BD2(:)

do ISYM=1,NSYM
  if (NINDEP(ISYM,5) == 0) cycle
  NAS = NTU(ISYM)
  NBD = nTri_Elem(2*NAS)
  if (NBD > 0) then
    call mma_Allocate(bD,NBD,LABEL='BD')
    call mma_Allocate(bD1,nTri_Elem(NAS),LABEL='BD1')
    call mma_Allocate(bD2,nTri_Elem(NAS),LABEL='BD2')
  end if
  BD(:) = Zero
  BD1(:) = Zero
  BD2(:) = Zero
  call Init_Tsk(ID,Nas)

  !do ITU=1,NAS
  !  if (.not. Rsv_Tsk(ID,ITU)) cycle
  do while (Rsv_Tsk(id,ITU))
    !ITU2 = ITU+NAS
    ITUABS = ITU+NTUes(ISYM)
    ITABS = MTU(1,ITuabS)
    IUABS = MTU(2,ITuabS)
    do IXY=1,ITU
      !IXY2 = IXY+NAS
      IXYABS = IXY+NtueS(ISYM)
      IXABS = MTU(1,ixyABS)
      IYABS = MTU(2,ixyABS)
      IB11 = iTri(ITU,IXY)
      !IB21 = iTri(ITU2,IXY)
      !IB12 = iTri(IXY2,ITU)
      !IB22 = iTri(ITU2,IXY2)

      !! A matrix
      B11 = Zero
      do IV=1,NASHT
        IVABS = IV
        !B11 = B11+hbar(IVABS,IXABS)*G2(IUABS,ITABS,IVABS,IYABS)
        !if (ITABS == IVABS) B11 = B11+Hbar(IVABS,IXABS)*G1(IUABS,IYABS)
        !B11 = B11-hbar(IYABS,IVABS)*G2(IUABS,ITABS,IXABS,IVABS)
        !if (ITABS == IXABS) B11 = B11-Hbar(IYABS,IVABS)*G1(IUABS,IVABS)
        B11 = B11+Hbar(IVABS,IXABS)*ex2(IUABS,ITABS,IVABS,IYABS,G1,G2)
        B11 = B11-Hbar(IYABS,IVABS)*ex2(IUABS,ITABS,IXABS,IVABS,G1,G2)
      end do

      !! D matrix (hbar terms)
      B22 = Zero
      do IV=1,NASHT
        IVABS = IV
        !! tilde{E}_{ab} = 2 \delta_{a,b} - E_{ba}
        tmp = -ex2(IxabS,ITABS,IUABS,IVABS,G1,G2)
        if (ITABS == ixABS) tmp = tmp+Two*G1(IUABS,IVABS)
        if (IUABS == itABS) tmp = tmp+G1(IVABS,IXABS)
        B22 = B22-Hbar(IYABS,IVABS)*tmp
        tmp = -ex2(IvabS,ITABS,IUABS,IYABS,G1,G2)
        if (ITABS == ivABS) tmp = tmp+Two*G1(IUABS,IYABS)
        if (IUABS == itABS) tmp = tmp+G1(IYABS,IVABS)
        B22 = B22+Hbar(IVABS,IXABS)*tmp
      end do

      !! D matrix (erI terms)
      do IV=1,NASHT
        IVABS = IV
        do IW=1,NASHt
          IWABS = IW
          do IZ=1,NAsht
            IZABS = iz
            tmp = Zero
            if (ITABs == IXABS) tmp = tmp+Two*ex2(IVABS,IZABS,IUABS,IWABS,G1,G2)
            if (ITABs == IXABS) tmp = tmp+Two*ex2(IUABS,IWABS,IVABS,IZABS,G1,G2)
            if (IUABs == ITABS) tmp = tmp+ex2(IVABS,IZABS,IXABS,IWABS,G1,G2)
            if (IUABs == ITABS) tmp = tmp+ex2(IXABS,IWABS,IVABS,IZABS,G1,G2)
            if ((ITABs == IVABS) .and. (IZABS == IXABS)) tmp = tmp+Two*G1(IUABS,IWABS)
            if (ITABs == IVABS) tmp = tmp-ex2(IXABS,IZABS,IUABS,IWABS,G1,G2)
            if ((IUABs == IZABS) .and. (IXABS == ITABS)) tmp = tmp-Two*G1(IVABS,IWABS)
            if (IUABs == IZABS) tmp = tmp+ex2(IXABS,ITABS,IVABS,IWABS,G1,G2)
            B22 = B22-half*tmp*Gact(IVABS,IZABS,IYABS,IWABS)

            tmp = Zero
            if (ITABs == IWABS) tmp = tmp+Two*ex2(IVABS,IZABS,IUABS,IYABS,G1,G2)
            if (ITABs == IWABS) tmp = tmp+Two*ex2(IUABS,IYABS,IVABS,IZABS,G1,G2)
            if (IUABs == ITABS) tmp = tmp+ex2(IVABS,IZABS,IWABS,IYABS,G1,G2)
            if (IUABs == ITABS) tmp = tmp+ex2(IWABS,IYABS,IVABS,IZABS,G1,G2)
            if ((ITABs == IVABS) .and. (IZABS == IWABS)) tmp = tmp+Two*G1(IUABS,IYABS)
            if (ITABs == IVABS) tmp = tmp-ex2(IWABS,IZABS,IUABS,IYABS,G1,G2)
            if ((IUABs == IZABS) .and. (ITABS == IWABS)) tmp = tmp-Two*G1(IVABS,IYABS)
            if (IUABs == IZABS) tmp = tmp+ex2(IWABS,ITABS,IVABS,IYABS,G1,G2)
            B22 = B22+half*tmp*Gact(IVABS,IZABS,IWABS,IXABS)
          end do
        end do
      end do

      BD1(IB11) = B11
      BD2(IB11) = B22
    end do
  end do
  call Free_Tsk(ID)

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
    ituvs = MUL(IST,mul(ISU,ISV))
    ixyzs = MUL(ISX,mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU = iT+NASHT*(iu-1)
    iVX = iV+NASHT*(ix-1)
    iYZ = iY+NASHT*(iz-1)

    ! G3(tuvxyz) = <psi|t+ v+ y+ z x u|Psi>
    !G3VAL = tG3(IT,iu,IV,IX,IY,IZ,G1,G2,G3(iG3))

    ! D(tuvxyz)
    if ((MUL(iSU,iST) == iSym) .or. (MUL(iSU,ISV) == iSym) .or. (MUL(iSX,iSY) == iSym)) then
      G3VAL = ex3(IT,iu,IV,IX,IY,IZ,G1,G2,G3(iG3))
      call Add_MKBNEvd(iT,iU,iV,iX,iY,iZ)
    end if

    if ((iTU /= iVX) .or. (iVX /= iYZ)) then
      if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
        ! D(vxtuyz)
        if ((MUL(iSX,isv) == iSym) .or. (MUL(iSX,iST) == iSym) .or. (MUL(iSU,iSY) == iSym)) then
          G3VAL = ex3(iV,IX,IT,IU,IY,IZ,G1,G2,G3(iG3))
          call Add_MkbnEVD(iV,iX,iT,iU,iY,iZ)
        end if
        ! D(yzvxtu)
        if ((MUL(iSZ,isy) == iSym) .or. (MUL(iSZ,iSV) == iSym) .or. (MUL(iSX,iST) == iSym)) then
          G3VAL = ex3(iY,IZ,IV,IX,IT,IU,G1,G2,G3(iG3))
          call Add_MkbnEVD(iY,iZ,iV,iX,iT,iU)
        end if
        ! D(tuyzvx)
        if ((MUL(iSU,ist) == iSym) .or. (MUL(iSU,iSY) == iSym) .or. (MUL(iSZ,iSV) == iSym)) then
          G3VAL = ex3(iT,IU,IY,IZ,IV,IX,G1,G2,G3(iG3))
          call Add_MkbnEVD(iT,iU,iY,iZ,iV,iX)
        end if
      end if
      ! D(yztuvx)
      if ((MUL(iSZ,iSy) == iSym) .or. (MUL(iSZ,iST) == iSym) .or. (MUL(iSU,iSV) == iSym)) then
        G3VAL = ex3(iy,IZ,IT,IU,IV,IX,G1,G2,G3(iG3))
        call Add_MKBnevD(iY,iZ,iT,iU,iV,iX)
      end if
      ! D(vxyztu)
      if ((MUL(iSX,iSv) == iSym) .or. (MUL(iSX,iSY) == iSym) .or. (MUL(iSZ,iST) == iSym)) then
        G3VAL = ex3(iv,IX,IY,IZ,IT,IU,G1,G2,G3(iG3))
        call Add_MKBnevD(iV,iX,iY,iZ,iT,iU)
      end if
    end if

    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle

    ! D(utxvzy)
    if ((MUL(iST,iSU) == iSym) .or. (MUL(iST,iSX) == iSym) .or. (MUL(iSV,iSZ) == iSym)) then
      G3VAL = ex3(IU,it,IX,IV,IZ,IY,G1,G2,G3(iG3))
      call Add_MKBNEvd(iU,iT,iX,iV,iZ,iY)
    end if

    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
      ! D(xvutzy)
      if ((MUL(iSV,iSx) == iSym) .or. (MUL(iSV,iSU) == iSym) .or. (MUL(iST,iSZ) == iSym)) then
        G3VAL = ex3(ix,IV,IU,IT,IZ,IY,G1,G2,G3(iG3))
        call Add_MKBnevD(iX,iV,iU,iT,iZ,iY)
      end if
      ! D(zyxvut)
      if ((MUL(iSY,iSz) == iSym) .or. (MUL(iSY,iSX) == iSym) .or. (MUL(iSV,iSU) == iSym)) then
        G3VAL = ex3(iz,IY,IX,IV,IU,IT,G1,G2,G3(iG3))
        call Add_MKBnevD(iZ,iY,iX,iV,iU,iT)
      end if
      ! D(utzyxv)
      if ((MUL(iST,iSu) == iSym) .or. (MUL(iST,iSZ) == iSym) .or. (MUL(iSY,iSX) == iSym)) then
        G3VAL = ex3(iu,IT,IZ,IY,IX,IV,G1,G2,G3(iG3))
        call Add_MKBnevD(iU,iT,iZ,iY,iX,iV)
      end if
    end if
    ! D(zyutxv)
    if ((MUL(iSY,iSZ) == iSym) .or. (MUL(iSY,iSU) == iSym) .or. (MUL(iST,iSX) == iSym)) then
      G3VAL = ex3(IZ,iy,IU,IT,IX,IV,G1,G2,G3(iG3))
      call Add_MKBNEvd(iZ,iY,iU,iT,iX,iV)
    end if
    ! D(xvzyut)
    if ((MUL(iSV,iSX) == iSym) .or. (MUL(iSV,iSZ) == iSym) .or. (MUL(iSY,iSU) == iSym)) then
      G3VAL = ex3(IX,iv,IZ,IY,IU,IT,G1,G2,G3(iG3))
      call Add_MKBNEvd(iX,iV,iZ,iY,iU,iT)
    end if
  end do

  do ITU=1,NAS
    ITU2 = ITU+NAS
    ITUABS = ITU+NTUes(ISYM)
    ITABS = MTU(1,ITuabS)
    IUABS = MTU(2,ITuabS)
    do IXY=1,ITU
      IXY2 = IXY+NAS
      IXYABS = IXY+NtueS(ISYM)
      IXABS = MTU(1,ixyABS)
      IYABS = MTU(2,ixyABS)
      IB11 = iTri(ITU,IXY)
      IB21 = iTri(ITU2,IXY)
      IB12 = iTri(IXY2,ITU)
      IB22 = iTri(ITU2,IXY2)

      B11 = Two*BD1(ib11)
      B22 = BD2(IB11)

      BD(IB11) = B11
      BD(IB21) = -Half*B11
      BD(IB12) = -Half*B11
      BD(IB22) = B22
    end do
  end do

# ifdef _MOLCAS_MPP_
  if (is_real_par()) call GADGOP(BD,NBD,'+')
# endif

  if ((NBD > 0) .and. (ninDEP(ISYM,5) > 0)) then
    IDISK = IDBMAT(Isym,5)
    call DDAFILE(LUSbt,1,BD,NBD,IDISK)
    call mma_deallocate(BD)
    call mma_deallocate(BD1)
    call mma_deallocate(BD2)
  end if
end do

contains

subroutine Add_MKBNEvd(ITABS_,IUABS_,IVABS_,IXABS_,IYABS_,IZABS_)

  integer(kind=iwp), intent(in) :: ITABS_, IUABS_, IVABS_, IXABS_, IYABS_, IZABS_
  integer(kind=iwp) :: IBADR, ITU_, IWABS, IXY_

  !! A
  iTU_ = KTU(IUABS_,itaBS_)-NTUES(ISYM)
  if (MUL(IASYM(IUABs_),IASYM(ITABS_)) == iSym) then
    do IWABS=1,nAshT
      !! A(a'b',ab) = (ce|da) E(b'a',ce,db) --> A(ut,wz) = (vx|yw) E(tu,vx,yz)
      if (MUL(IASYM(iwaBS),IASYM(IZABS_)) == iSym) then
        iXY_ = KTU(IwabS,IZABS_)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD1(IBADR) = BD1(IBADR)+Half*G3VAL*Gact(IVABS_,IXABS_,IYABS_,IWABS)
        end if
      end if
      !! A(a'b',ab) = (ce|da) E(b'a',db,ce) --> A(ut,wx) = (yz|vw) E(tu,vx,yz)
      if (MUL(IASYM(iwaBS),IASYM(IXABS_)) == iSym) then
        iXY_ = KTU(IwabS,IXABS_)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD1(IBADR) = BD1(IBADR)+Half*G3VAL*Gact(IYABS_,IZABS_,IVABS_,IWABS)
        end if
      end if

      !! A(a'b',ab) = (be|cd) E(b'a',ae,cd) --> A(ut,vw) = (wx|yz) E(tu,vx,yz)
      if (MUL(IASYM(ivaBS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IvabS_,IWABS)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD1(IBADR) = BD1(IBADR)-Half*G3VAL*Gact(IWABS,IXABS_,IYABS_,IZABS_)
        end if
      end if
      !! A(a'b',ab) = (be|cd) E(b'a',cd,ae) --> A(ut,yw) = (wz|vx) E(tu,vx,yz)
      if (MUL(IASYM(iyaBS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IyabS_,IWABS)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD1(IBADR) = BD1(IBADR)-Half*G3VAL*Gact(IWABS,IZABS_,IVABS_,IXABS_)
        end if
      end if
    end do
  end if

  !! D
  do IWABS=1,nAshT
    !! D(a'b',ab) = -(ce|bd) E(ce,aa',b'd) --> D(xy,vw) = -(tu|wz) E(tu,vx,yz)
    iTU_ = KTU(IXABS_,iYABS_)-NTUES(ISYM)
    if (MUL(IASYM(IXabs_),IASYM(IYABS_)) == iSym) then
      if (MUL(IASYM(ivaBS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IvabS_,IWABS)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD2(IBADR) = BD2(IBADR)+Half*G3VAL*Gact(ITABS_,IUABS_,IWABS,IZABS_)
        end if
      end if
    end if
    !! D(a'b',ab) = -(ce|bd) E(aa',b'd,ce) --> D(uv,tw) = -(yz|wx) E(tu,vx,yz)
    iTU_ = KTU(IUABS_,iVABS_)-NTUES(ISYM)
    if (MUL(IASYM(IUabs_),IASYM(IVABS_)) == iSym) then
      if (MUL(IASYM(itaBS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(ItabS_,IWABS)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD2(IBADR) = BD2(IBADR)+Half*G3VAL*Gact(IYABS_,IZABS_,IWABS,IXABS_)
        end if
      end if
    end if

    !! D(a'b',ab) = +(ce|da) E(ce,da',b'b) --> D(xy,wz) = +(tu|vw) E(tu,vx,yz)
    iTU_ = KTU(IXABS_,iYABS_)-NTUES(ISYM)
    if (MUL(IASYM(IXabs_),IASYM(IYABS_)) == iSym) then
      if (MUL(IASYM(iwaBS),IASYM(IZABS_)) == iSym) then
        iXY_ = KTU(IwabS,IZABS_)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD2(IBADR) = BD2(IBADR)-Half*G3VAL*Gact(ITABS_,IUABS_,IVABS_,IWABS)
        end if
      end if
    end if
    !! D(a'b',ab) = +(ce|da) E(da',b'b,ce) --> D(uv,wx) = +(yz|tw) E(tu,vx,yz)
    iTU_ = KTU(IUABS_,iVABS_)-NTUES(ISYM)
    if (MUL(IASYM(IUabs_),IASYM(IVABS_)) == iSym) then
      if (MUL(IASYM(iwaBS),IASYM(IXABS_)) == iSym) then
        iXY_ = KTU(IwabS,IXABS_)-NTUES(ISYM)
        if (iTU_ >= ixy_) then
          IBADR = iTri(iTu_,iXY_)
          BD2(IBADR) = BD2(IBADR)-Half*G3VAL*Gact(IYABS_,IZABS_,ITABS_,IWABS)
        end if
      end if
    end if
  end do

end subroutine Add_MkbnEVD

end subroutine MKBNEvd
