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
  use caspt2_global, only: LUSBT
  use caspt2_module, only: IASYM, NINDEP, NSYM, NTU, NTUES
  use Constants, only: Zero, Half, Two
  use definitions, only: iwp,wp,byte
  use EQSOLV, only: IDBMAT
  use stdalloc, only: mma_allocate, mma_deallocate
  USE SUPERINDEX, only : KTU, MTU
  use Symmetry_Info, only: Mul
  use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
#ifdef _MOLCAS_MPP_
  USE Para_Info, ONLY: Is_Real_Par
#endif

  use nevpt2_mod, only: ex2, ex3

  implicit none

  integer(kind=iwp), intent(in) :: nAshT, NG3
  real(kind=wp), intent(in) :: Hbar(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), &
                               G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(*)
  integer(kind=byte), intent(in) :: idxG3(6,NG3)

  integer(kind=iwp) :: IB11, IB12, IB21, IB22, IDISK, ISYM, IT, ITABS, ITU, ITUABS, ITU2, &
                       IU, IUABS, IV, IVABS, IVX, IW, IWABS, IX, IXABS, IXY, IXYABS, IXY2, IY, &
                       IYABS, IZ, IZABS, NAS, NBD
  integer(kind=iwp) :: iG3, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IXYZS, IYZ
  integer(kind=iwp) :: ID
  real(kind=wp) :: B11, B22, G3VAL, tmp

  real(kind=wp), allocatable :: BD(:), BD1(:), BD2(:)

  DO ISYM=1,NSYM
    IF(NINDEP(ISYM,5) == 0) cycle
    NAS=NTU(ISYM)
    NBD=(2*NAS*(2*NAS+1))/2
    if (NBD > 0) then
      CALL mma_Allocate(BD,NBD,LABEL='BD')
      CALL mma_Allocate(BD1,NAS*(NAS+1)/2,LABEL='BD1')
      CALL mma_Allocate(BD2,NAS*(NAS+1)/2,LABEL='BD2')
    end if
    BD(:) = Zero
    BD1(:) = Zero
    BD2(:) = Zero
    Call Init_Tsk(ID,NAS)

!   DO ITU=1,NAS
!     if (.not.Rsv_Tsk(ID,ITU)) cycle
    do while (Rsv_Tsk(ID,ITU))
!     ITU2=ITU+NAS
      ITUABS=ITU+NTUES(ISYM)
      ITABS=MTU(1,ITUABS)
      IUABS=MTU(2,ITUABS)
      DO IXY=1,ITU
!       IXY2=IXY+NAS
        IXYABS=IXY+NTUES(ISYM)
        IXABS=MTU(1,IXYABS)
        IYABS=MTU(2,IXYABS)
        IB11=(ITU*(ITU-1))/2+IXY
!       IB21=(ITU2*(ITU2-1))/2+IXY
!       IB12=(IXY2*(IXY2-1))/2+ITU
!       IB22=(ITU2*(ITU2-1))/2+IXY2

        !! A matrix
        B11 = Zero
        DO IV = 1, NASHT
          IVABS = IV
!         B11 = B11 + Hbar(IVABS,IXABS)*G2(IUABS,ITABS,IVABS,IYABS)
!         IF (ITABS == IVABS) B11 = B11 + Hbar(IVABS,IXABS)*G1(IUABS,IYABS)
!         B11 = B11 - Hbar(IYABS,IVABS)*G2(IUABS,ITABS,IXABS,IVABS)
!         IF (ITABS == IXABS) B11 = B11 - Hbar(IYABS,IVABS)*G1(IUABS,IVABS)
          B11 = B11 + Hbar(IVABS,IXABS)*ex2(IUABS,ITABS,IVABS,IYABS,G1,G2)
          B11 = B11 - Hbar(IYABS,IVABS)*ex2(IUABS,ITABS,IXABS,IVABS,G1,G2)
        END DO

        !! D matrix (hbar terms)
        B22 = Zero
        DO IV = 1, NASHT
          IVABS = IV
          !! tilde{E}_{ab} = 2 \delta_{a,b} - E_{ba}
          tmp = -ex2(IXABS,ITABS,IUABS,IVABS,G1,G2)
          if (ITABS == IXABS) tmp = tmp + Two*G1(IUABS,IVABS)
          if (IUABS == ITABS) tmp = tmp + G1(IVABS,IXABS)
          B22 = B22 - Hbar(IYABS,IVABS)*tmp
          tmp = -ex2(IVABS,ITABS,IUABS,IYABS,G1,G2)
          if (ITABS == IVABS) tmp = tmp + Two*G1(IUABS,IYABS)
          if (IUABS == ITABS) tmp = tmp + G1(IYABS,IVABS)
          B22 = B22 + Hbar(IVABS,IXABS)*tmp
        END DO

        !! D matrix (ERI terms)
        DO IV = 1, NASHT
          IVABS = IV
          DO IW = 1, NASHT
            IWABS = IW
            DO IZ = 1, NASHT
              IZABS = IZ
              tmp = Zero
              if (ITABS == IXABS) tmp = tmp + Two*ex2(IVABS,IZABS,IUABS,IWABS,G1,G2)
              if (ITABS == IXABS) tmp = tmp + Two*ex2(IUABS,IWABS,IVABS,IZABS,G1,G2)
              if (IUABS == ITABS) tmp = tmp + ex2(IVABS,IZABS,IXABS,IWABS,G1,G2)
              if (IUABS == ITABS) tmp = tmp + ex2(IXABS,IWABS,IVABS,IZABS,G1,G2)
              if (ITABS == IVABS .and. IZABS == IXABS) tmp = tmp + Two*G1(IUABS,IWABS)
              if (ITABS == IVABS) tmp = tmp - ex2(IXABS,IZABS,IUABS,IWABS,G1,G2)
              if (IUABS == IZABS .and. IXABS == ITABS) tmp = tmp - Two*G1(IVABS,IWABS)
              if (IUABS == IZABS) tmp = tmp + ex2(IXABS,ITABS,IVABS,IWABS,G1,G2)
              B22 = B22 - Half*tmp*Gact(IVABS,IZABS,IYABS,IWABS)

              tmp = Zero
              if (ITABS == IWABS) tmp = tmp + Two*ex2(IVABS,IZABS,IUABS,IYABS,G1,G2)
              if (ITABS == IWABS) tmp = tmp + Two*ex2(IUABS,IYABS,IVABS,IZABS,G1,G2)
              if (IUABS == ITABS) tmp = tmp + ex2(IVABS,IZABS,IWABS,IYABS,G1,G2)
              if (IUABS == ITABS) tmp = tmp + ex2(IWABS,IYABS,IVABS,IZABS,G1,G2)
              if (ITABS == IVABS .and. IZABS == IWABS) tmp = tmp + Two*G1(IUABS,IYABS)
              if (ITABS == IVABS) tmp = tmp - ex2(IWABS,IZABS,IUABS,IYABS,G1,G2)
              if (IUABS == IZABS .and. ITABS == IWABS) tmp = tmp - Two*G1(IVABS,IYABS)
              if (IUABS == IZABS) tmp = tmp + ex2(IWABS,ITABS,IVABS,IYABS,G1,G2)
              B22 = B22 + Half*tmp*Gact(IVABS,IZABS,IWABS,IXABS)
            END DO
          END DO
        END DO

        BD1(IB11) = B11
        BD2(IB11) = B22
      END DO
    END DO
    CALL Free_Tsk(ID)

    do IG3 = 1, NG3
      iT=idxG3(1,iG3)
      iU=idxG3(2,iG3)
      iV=idxG3(3,iG3)
      iX=idxG3(4,iG3)
      iY=idxG3(5,iG3)
      iZ=idxG3(6,iG3)
      iST=IASYM(iT)
      iSU=IASYM(iU)
      iSV=IASYM(iV)
      iSX=IASYM(iX)
      iSY=IASYM(iY)
      iSZ=IASYM(iZ)
      ituvs=MUL(IST,MUL(ISU,ISV))
      ixyzs=MUL(ISX,MUL(ISY,ISZ))
      if (ituvs /= ixyzs) cycle
      iTU=iT+NASHT*(iU-1)
      iVX=iV+NASHT*(iX-1)
      iYZ=iY+NASHT*(iZ-1)

      ! G3(tuvxyz) = <Psi|t+ v+ y+ z x u|Psi>
!     G3VAL = tG3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))

      ! D(tuvxyz)
      if (MUL(iSU,iST) == iSym .or. MUL(iSU,ISV) == iSym .or. MUL(iSX,iSY) == iSym) then
        G3VAL = ex3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))
        call Add_MKBNEVD(iT,iU,iV,iX,iY,iZ)
      end if

      if (iTU /= iVX .or. iVX /= iYZ) then
        if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
          ! D(vxtuyz)
          if (MUL(iSX,iSV) == iSym .or. MUL(iSX,iST) == iSym .or. MUL(iSU,iSY) == iSym) then
            G3VAL = ex3(IV,IX,IT,IU,IY,IZ,G1,G2,G3(iG3))
            call Add_MKBNEVD(iV,iX,iT,iU,iY,iZ)
          end if
          ! D(yzvxtu)
          if (MUL(iSZ,iSY) == iSym .or. MUL(iSZ,iSV) == iSym .or. MUL(iSX,iST) == iSym) then
            G3VAL = ex3(IY,IZ,IV,IX,IT,IU,G1,G2,G3(iG3))
            call Add_MKBNEVD(iY,iZ,iV,iX,iT,iU)
          end if
          ! D(tuyzvx)
          if (MUL(iSU,iST) == iSym .or. MUL(iSU,iSY) == iSym .or. MUL(iSZ,iSV) == iSym) then
            G3VAL = ex3(IT,IU,IY,IZ,IV,IX,G1,G2,G3(iG3))
            call Add_MKBNEVD(iT,iU,iY,iZ,iV,iX)
          end if
        end if
        ! D(yztuvx)
        if (MUL(iSZ,iSY) == iSym .or. MUL(iSZ,iST) == iSym .or. MUL(iSU,iSV) == iSym) then
          G3VAL = ex3(IY,IZ,IT,IU,IV,IX,G1,G2,G3(iG3))
          call Add_MKBNEVD(iY,iZ,iT,iU,iV,iX)
        end if
        ! D(vxyztu)
        if (MUL(iSX,iSV) == iSym .or. MUL(iSX,iSY) == iSym .or. MUL(iSZ,iST) == iSym) then
          G3VAL = ex3(IV,IX,IY,IZ,IT,IU,G1,G2,G3(iG3))
          call Add_MKBNEVD(iV,iX,iY,iZ,iT,iU)
        end if
      end if

      if (iT == iU .and. iV == iX .and. iY == iZ) cycle
      if (iT == iU .and. iV == iZ .and. iX == iY) cycle
      if (iX == iV .and. iT == iZ .and. iU == iY) cycle
      if (iZ == iY .and. iV == iU .and. iX == iT) cycle

      ! D(utxvzy)
      if (MUL(iST,iSU) == iSym .or. MUL(iST,iSX) == iSym .or. MUL(iSV,iSZ) == iSym) then
        G3VAL = ex3(IU,IT,IX,IV,IZ,IY,G1,G2,G3(iG3))
        call Add_MKBNEVD(iU,iT,iX,iV,iZ,iY)
      end if

      if (iTU == iVX .and. iVX == iYZ) cycle
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
        ! D(xvutzy)
        if (MUL(iSV,iSX) == iSym .or. MUL(iSV,iSU) == iSym .or. MUL(iST,iSZ) == iSym) then
          G3VAL = ex3(IX,IV,IU,IT,IZ,IY,G1,G2,G3(iG3))
          call Add_MKBNEVD(iX,iV,iU,iT,iZ,iY)
        end if
        ! D(zyxvut)
        if (MUL(iSY,iSZ) == iSym .or. MUL(iSY,iSX) == iSym .or. MUL(iSV,iSU) == iSym) then
          G3VAL = ex3(IZ,IY,IX,IV,IU,IT,G1,G2,G3(iG3))
          call Add_MKBNEVD(iZ,iY,iX,iV,iU,iT)
        end if
        ! D(utzyxv)
        if (MUL(iST,iSU) == iSym .or. MUL(iST,iSZ) == iSym .or. MUL(iSY,iSX) == iSym) then
          G3VAL = ex3(IU,IT,IZ,IY,IX,IV,G1,G2,G3(iG3))
          call Add_MKBNEVD(iU,iT,iZ,iY,iX,iV)
        end if
      end if
      ! D(zyutxv)
      if (MUL(iSY,iSZ) == iSym .or. MUL(iSY,iSU) == iSym .or. MUL(iST,iSX) == iSym) then
        G3VAL = ex3(IZ,IY,IU,IT,IX,IV,G1,G2,G3(iG3))
        call Add_MKBNEVD(iZ,iY,iU,iT,iX,iV)
      end if
      ! D(xvzyut)
      if (MUL(iSV,iSX) == iSym .or. MUL(iSV,iSZ) == iSym .or. MUL(iSY,iSU) == iSym) then
        G3VAL = ex3(IX,IV,IZ,IY,IU,IT,G1,G2,G3(iG3))
        call Add_MKBNEVD(iX,iV,iZ,iY,iU,iT)
      end if
    end do

    DO ITU=1,NAS
      ITU2=ITU+NAS
      ITUABS=ITU+NTUES(ISYM)
      ITABS=MTU(1,ITUABS)
      IUABS=MTU(2,ITUABS)
      DO IXY=1,ITU
        IXY2=IXY+NAS
        IXYABS=IXY+NTUES(ISYM)
        IXABS=MTU(1,IXYABS)
        IYABS=MTU(2,IXYABS)
        IB11=(ITU*(ITU-1))/2+IXY
        IB21=(ITU2*(ITU2-1))/2+IXY
        IB12=(IXY2*(IXY2-1))/2+ITU
        IB22=(ITU2*(ITU2-1))/2+IXY2

        B11 = Two*BD1(IB11)
        B22 =     BD2(IB11)

        BD(IB11)= B11
        BD(IB21)=-Half*B11
        BD(IB12)=-Half*B11
        BD(IB22)= B22
      END DO
    END DO

#ifdef _MOLCAS_MPP_
    if (is_real_par()) call GADGOP(BD,NBD,'+')
#endif

    IF(NBD > 0 .and. NINDEP(ISYM,5) > 0) THEN
      IDISK=IDBMAT(ISYM,5)
      CALL DDAFILE(LUSBT,1,BD,NBD,IDISK)
      CALL mma_deallocate(BD)
      CALL mma_deallocate(BD1)
      CALL mma_deallocate(BD2)
    END IF
  end do

contains

  subroutine Add_MKBNEVD(ITABS_,IUABS_,IVABS_,IXABS_,IYABS_,IZABS_)

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IUABS_, IVABS_, IXABS_, IYABS_, IZABS_
  integer(kind=iwp) :: IBADR, ITU_, IWABS, IXY_

  !! A
  iTU_ = KTU(IUABS_,ITABS_)-NTUES(ISYM)
  if (MUL(IASYM(IUABS_),IASYM(ITABS_)) == iSym) then
    do IWABS = 1, nAshT
      !! A(a'b',ab) = (ce|da) E(b'a',ce,db) --> A(ut,wz) = (vx|yw) E(tu,vx,yz)
      if (MUL(IASYM(IWABS),IASYM(IZABS_)) == iSym) then
        iXY_ = KTU(IWABS,IZABS_)-NTUES(ISYM)
        if (iTU_ >= iXY_) then
          IBADR = iTU_*(iTU_-1)/2 + iXY_
          BD1(IBADR) = BD1(IBADR) + Half*G3VAL*Gact(IVABS_,IXABS_,IYABS_,IWABS)
        end if
      end if
      !! A(a'b',ab) = (ce|da) E(b'a',db,ce) --> A(ut,wx) = (yz|vw) E(tu,vx,yz)
      if (MUL(IASYM(IWABS),IASYM(IXABS_)) == iSym) then
        iXY_ = KTU(IWABS,IXABS_)-NTUES(ISYM)
        if (iTU_ >= iXY_) then
          IBADR = iTU_*(iTU_-1)/2 + iXY_
          BD1(IBADR) = BD1(IBADR) + Half*G3VAL*Gact(IYABS_,IZABS_,IVABS_,IWABS)
        end if
      end if

      !! A(a'b',ab) = (be|cd) E(b'a',ae,cd) --> A(ut,vw) = (wx|yz) E(tu,vx,yz)
      if (MUL(IASYM(IVABS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IVABS_,IWABS)-NTUES(ISYM)
        if (iTU_ >= iXY_) then
          IBADR = iTU_*(iTU_-1)/2 + iXY_
          BD1(IBADR) = BD1(IBADR) - Half*G3VAL*Gact(IWABS,IXABS_,IYABS_,IZABS_)
        end if
      end if
      !! A(a'b',ab) = (be|cd) E(b'a',cd,ae) --> A(ut,yw) = (wz|vx) E(tu,vx,yz)
      if (MUL(IASYM(IYABS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IYABS_,IWABS)-NTUES(ISYM)
        if (iTU_ >= iXY_) then
          IBADR = iTU_*(iTU_-1)/2 + iXY_
          BD1(IBADR) = BD1(IBADR) - Half*G3VAL*Gact(IWABS,IZABS_,IVABS_,IXABS_)
        end if
      end if
    end do
  end if

  !! D
  do IWABS = 1, nAshT
    !! D(a'b',ab) = -(ce|bd) E(ce,aa',b'd) --> D(xy,vw) = -(tu|wz) E(tu,vx,yz)
    iTU_ = KTU(IXABS_,IYABS_)-NTUES(ISYM)
    if (MUL(IASYM(IXABS_),IASYM(IYABS_)) == iSym) then
    if (MUL(IASYM(IVABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IVABS_,IWABS)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BD2(IBADR) = BD2(IBADR) + Half*G3VAL*Gact(ITABS_,IUABS_,IWABS,IZABS_)
      end if
    end if
    end if
    !! D(a'b',ab) = -(ce|bd) E(aa',b'd,ce) --> D(uv,tw) = -(yz|wx) E(tu,vx,yz)
    iTU_ = KTU(IUABS_,IVABS_)-NTUES(ISYM)
    if (MUL(IASYM(IUABS_),IASYM(IVABS_)) == iSym) then
    if (MUL(IASYM(ITABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(ITABS_,IWABS)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BD2(IBADR) = BD2(IBADR) + Half*G3VAL*Gact(IYABS_,IZABS_,IWABS,IXABS_)
      end if
    end if
    end if

    !! D(a'b',ab) = +(ce|da) E(ce,da',b'b) --> D(xy,wz) = +(tu|vw) E(tu,vx,yz)
    iTU_ = KTU(IXABS_,IYABS_)-NTUES(ISYM)
    if (MUL(IASYM(IXABS_),IASYM(IYABS_)) == iSym) then
    if (MUL(IASYM(IWABS),IASYM(IZABS_)) == iSym) then
      iXY_ = KTU(IWABS,IZABS_)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BD2(IBADR) = BD2(IBADR) - Half*G3VAL*Gact(ITABS_,IUABS_,IVABS_,IWABS)
      end if
    end if
    end if
    !! D(a'b',ab) = +(ce|da) E(da',b'b,ce) --> D(uv,wx) = +(yz|tw) E(tu,vx,yz)
    iTU_ = KTU(IUABS_,IVABS_)-NTUES(ISYM)
    if (MUL(IASYM(IUABS_),IASYM(IVABS_)) == iSym) then
    if (MUL(IASYM(IWABS),IASYM(IXABS_)) == iSym) then
      iXY_ = KTU(IWABS,IXABS_)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BD2(IBADR) = BD2(IBADR) - Half*G3VAL*Gact(IYABS_,IZABS_,ITABS_,IWABS)
      end if
    end if
    end if
  end do

  end subroutine Add_MKBNEVD

end subroutine MKBNEVD
