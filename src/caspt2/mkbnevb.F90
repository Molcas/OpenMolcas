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
  use caspt2_global, only: LUSBT
  use caspt2_module, only: IASYM, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
  use Constants, only: Zero, Two, Half
  use definitions, only: iwp,wp,byte
  use EQSOLV, only: IDBMAT
  use NEVPT2_mod, only: tG2, tG3
  use stdalloc, only: mma_allocate, mma_deallocate
  USE SUPERINDEX, only: KTU, MTU, MTGEU, KTGTU
  use Symmetry_Info, only: Mul
  use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
#ifdef _MOLCAS_MPP_
  USE Para_Info, ONLY: Is_Real_Par
#endif

  implicit none

  integer(kind=iwp), intent(in) :: nAshT, NG3
  real(kind=wp), intent(in) :: Hbar(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), &
                               G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(*)
  integer(kind=byte), intent(in) :: idxG3(6,NG3)

  integer(kind=iwp) :: IBADR, IBMADR, IBPADR, IDISK, ISYM, ITABS, ITGEU, ITGEUABS, IT, ITGTU, ITU, ITUABS, &
                       IU, IUABS, IV, IVABS, IVX, IW, IWABS, IX, IXABS, IXGEY, IXGEYABS, IXGTY, IXY, IXYABS, IY, IYABS, IYX, &
                       IZ, IZABS, &
                       NAS, NASM, NASP, NBB, NBBM, NBBP
  integer(kind=iwp) :: iG3, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IXYZS, IYZ
! integer(kind=iwp) :: iSymT, iSymU, iSymX, iSymY
  integer(kind=iwp) :: ID
  real(kind=wp) :: BTUXY, BTUYX, G3VAL, VALUE

  real(kind=wp), allocatable :: BB(:), BBP(:), BBM(:)

  DO ISYM=1,NSYM
    IF(NINDEP(ISYM,2) == 0) cycle
    NAS=NTU(ISYM)
    NBB=(NAS*(NAS+1))/2
    if (NBB > 0) CALL mma_Allocate(BB,NBB,LABEL='BB')
    BB(:) = Zero
    Call Init_Tsk(ID,NAS)

    ! tilde{R}^{(2)} --> tG2
    ! The abcd <--> acbd symmetry is considered in the tG2 subroutine
!   DO ITU=1,NAS
!     if (.not.Rsv_Tsk(ID,ITU)) cycle
    do while (Rsv_Tsk(ID,ITU))
      ITUABS=ITU+NTUES(ISYM)
      ITABS=MTU(1,ITUABS)
      IUABS=MTU(2,ITUABS)
!     iSymT = IASYM(ITABS)
!     iSymU = IASYM(IUABS)
      DO IXY=1,ITU
        IXYABS=IXY+NTUES(ISYM)
        IXABS=MTU(1,IXYABS)
        IYABS=MTU(2,IXYABS)
!       iSymX = IASYM(IXABS)
!       iSymY = IASYM(IYABS)
        IBADR=(ITU*(ITU-1))/2+IXY
        VALUE = Zero
        ! Hbar term
        DO IV=1,NASHT
          IVABS=IV!+NAES(ISYM)
          VALUE = VALUE + Hbar(IVABS,IXABS)*tG2(ITABS,IUABS,IVABS,IYABS,G1,G2)
          VALUE = VALUE + Hbar(IVABS,IYABS)*tG2(ITABS,IUABS,IXABS,IVABS,G1,G2)
        END DO
        ! RDM2 term
        DO IV=1,NASHT
          IVABS=IV!+NAES(ISYM)
          DO IW=1,NASHT
            IWABS=IW!+NAES(ISYM)
            DO IZ=1,NASHT
              IZABS=IZ!+NAES(ISYM)
              if (IVABS == IZABS) VALUE = VALUE +  Two*Gact(IVABS,IZABS,IWABS,IXABS)*tG2(ITABS,IUABS,IWABS,IYABS,G1,G2)
              if (IWABS == IZABS) VALUE = VALUE - Half*Gact(IVABS,IZABS,IWABS,IXABS)*tG2(ITABS,IUABS,IVABS,IYABS,G1,G2)
              if (IYABS == IZABS) VALUE = VALUE - Half*Gact(IVABS,IZABS,IWABS,IXABS)*tG2(ITABS,IUABS,IWABS,IVABS,G1,G2)

              if (IVABS == IZABS) VALUE = VALUE +  Two*Gact(IVABS,IZABS,IWABS,IYABS)*tG2(ITABS,IUABS,IXABS,IWABS,G1,G2)
              if (IXABS == IZABS) VALUE = VALUE - Half*Gact(IVABS,IZABS,IWABS,IYABS)*tG2(ITABS,IUABS,IVABS,IWABS,G1,G2)
              if (IWABS == IZABS) VALUE = VALUE - Half*Gact(IVABS,IZABS,IWABS,IYABS)*tG2(ITABS,IUABS,IXABS,IVABS,G1,G2)
            END DO
          END DO
        END DO
        BB(IBADR)=VALUE
!       write (*,'(4i3,f20.10)') itabs,iuabs,ixabs,iyabs,value
      END DO
    END DO
    CALL Free_Tsk(ID)

    ! K (tu,xy) = -(vz|wx)*D(tuz,wyv) - (vz|wy)*D(tuz,xwv)
    ! K1(tu,wy) = -(zv|xw)*D(tuv,xyz)
    !>K1(tx,wv) = -(zu|yw)*D(txu,yvz)
    ! K (tu,xy) = -(vz|wy)*D(tuz,xwv)
    ! K2(tu,xw) = -(zv|yw)*D(tuv,xyz)
    !>K2(tx,yw) = -(zu|vw)*D(txu,yvz)
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
      G3VAL = tG3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))

      ! D(tuvxyz)
      if (MUL(iST,iSV) == iSym) call Add_MKBNEVB(iT,iU,iV,iX,iY,iZ)

      if (iTU /= iVX .or. iVX /= iYZ) then
        if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
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

      if (iT == iU .and. iV == iX .and. iY == iZ) cycle
      if (iT == iU .and. iV == iZ .and. iX == iY) cycle
      if (iX == iV .and. iT == iZ .and. iU == iY) cycle
      if (iZ == iY .and. iV == iU .and. iX == iT) cycle

      ! D(utxvzy)
      if (MUL(iSU,iSX) == iSym) call Add_MKBNEVB(iU,iT,iX,iV,iZ,iY)

      if (iTU == iVX .and. iVX == iYZ) cycle
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
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

    NASP=NTGEU(ISYM)
    NBBP=(NASP*(NASP+1))/2
    IF(NBBP > 0) CALL mma_allocate(BBP,NBBP,Label='BBP')
    NASM=NTGTU(ISYM)
    NBBM=(NASM*(NASM+1))/2
    IF(NBBM > 0) CALL mma_allocate(BBM,NBBM,Label='BBM')

    DO ITGEU=1,NASP
      ITGEUABS=ITGEU+NTGEUES(ISYM)
      ITABS=MTGEU(1,ITGEUABS)
      IUABS=MTGEU(2,ITGEUABS)
      ITU=KTU(ITABS,IUABS)-NTUES(ISYM)
      DO IXGEY=1,ITGEU
        IXGEYABS=IXGEY+NTGEUES(ISYM)
        IXABS=MTGEU(1,IXGEYABS)
        IYABS=MTGEU(2,IXGEYABS)
        IXY=KTU(IXABS,IYABS)-NTUES(ISYM)
        IYX=KTU(IYABS,IXABS)-NTUES(ISYM)
        IF(ITU >= IXY) THEN
          IBADR=(ITU*(ITU-1))/2+IXY
        ELSE
          IBADR=(IXY*(IXY-1))/2+ITU
        END IF
        BTUXY=BB(IBADR)
        IF(ITU >= IYX) THEN
          IBADR=(ITU*(ITU-1))/2+IYX
        ELSE
          IBADR=(IYX*(IYX-1))/2+ITU
        END IF
        BTUYX=BB(IBADR)
        IBPADR=(ITGEU*(ITGEU-1))/2+IXGEY
        BBP(IBPADR)=BTUXY+BTUYX
        IF(ITABS == IUABS) CYCLE
        IF(IXABS == IYABS) CYCLE
        ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
        IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
        IBMADR=(ITGTU*(ITGTU-1))/2+IXGTY
        BBM(IBMADR)=BTUXY-BTUYX
      END DO
    END DO
#ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      call GADGOP(BBP,NBBP,'+')
      call GADGOP(BBM,NBBM,'+')
    end if
#endif

    IF(NBB > 0) CALL mma_deallocate(BB)

    IF(NBBP > 0.and.NINDEP(ISYM,2) > 0) THEN
      IDISK=IDBMAT(ISYM,2)
      CALL DDAFILE(LUSBT,1,BBP,NBBP,IDISK)
      CALL mma_deallocate(BBP)
    END IF
    IF(NBBM > 0) THEN
     IF(NINDEP(ISYM,3) > 0) THEN
      IDISK=IDBMAT(ISYM,3)
      CALL DDAFILE(LUSBT,1,BBM,NBBM,IDISK)
     END IF
     CALL mma_deallocate(BBM)
    END IF
  end do

contains

  subroutine Add_MKBNEVB(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_)

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IVABS_, IYABS_, IUABS_, IXABS_, IZABS_
  integer(kind=iwp) :: IBADR, ITU_, IWABS, IXY_

  iTU_ = KTU(ITABS_,IUABS_)-NTUES(ISYM)
  do IWABS = 1, nAshT
    if (MUL(IASYM(IWABS),IASYM(IYABS_)) == iSym) then
      iXY_ = KTU(IWABS,IYABS_)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BB(IBADR) = BB(IBADR) - G3VAL*Gact(IZABS_,IVABS_,IXABS_,IWABS)
      end if
    end if
    if (MUL(IASYM(IXABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IXABS_,IWABS)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BB(IBADR) = BB(IBADR) - G3VAL*Gact(IZABS_,IVABS_,IYABS_,IWABS)
      end if
    end if
  end do

  end subroutine Add_MKBNEVB

end subroutine MKBNEVB
