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
  use caspt2_global, only: LUSBT
  use caspt2_module, only: IASYM, NINDEP, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES, NTU, NTUES
  use Constants, only: Zero, Two, Half
  use definitions, only: iwp,wp,byte
  use EQSOLV, only: IDBMAT
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
                               G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  integer(kind=byte), intent(in) :: idxG3(6,NG3)

  integer(kind=iwp) :: IBADR, IBMADR, IBPADR, IDISK, ISYM, ITABS, ITGEU, ITGEUABS, IT, ITGTU, ITU, ITUABS, ITX, &
                       IU, IUABS, IUY, IV, IVABS, IVX, IW, IWABS, IX, IXABS, IXGEY, IXGEYABS, IXGTY, IXY, IXYABS, IY, IYABS, IYX, &
                       IZ, IZABS, &
                       NAS, NASM, NASP, NBF, NBFM, NBFP
  integer(kind=iwp) :: iG3, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IXYZS, IYZ
! integer(kind=iwp) :: iSymT, iSymU, iSymX, iSymY
  integer(kind=iwp) :: ID
  real(kind=wp) :: BTUXY, BTUYX, G3VAL, VALUE

  real(kind=wp), allocatable :: BF(:), BFP(:), BFM(:)

  DO ISYM=1,NSYM
    IF(NINDEP(ISYM,8) == 0) cycle
    NAS=NTU(ISYM)
    NBF=(NAS*(NAS+1))/2
    if (NBF > 0) CALL mma_Allocate(BF,NBF,LABEL='BF')
    BF(:) = Zero
    Call Init_Tsk(ID,NAS)

    ! K(tu,xy) = - h(xv)D(tu,vy) - h(yv)D(tu,xv)
    ! K(tu,xy) = - 0.5*(vw|xz) ( del(vz)D(tu,wy) + del(yv)D(tu,zw) )
    !            - 0.5*(vw|yz) ( del(xv)D(tu,wz) + del(vz)D(tu,xw) )
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
        ITX=ITABS+NASHT*(IXABS-1)
        IUY=IUABS+NASHT*(IYABS-1)
!       IP1=MAX(ITX,IUY)
!       IP2=MIN(ITX,IUY)
!       IP=(IP1*(IP1-1))/2+IP2
        VALUE = Zero
        ! Hbar term
        DO IV=1,NASHT
          IVABS=IV!+NAES(ISYM)
          VALUE = VALUE - Hbar(IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
          VALUE = VALUE - Hbar(IYABS,IVABS)*G2(ITABS,IXABS,IUABS,IVABS)
        END DO
!       do IV = 1, nAsh(iSymX)
!         IVABS = IV + NAES(iSymX)
!         VALUE = VALUE - Hbar(IXABS,IVABS)*G2(ITABS,IVABS,IUABS,IYABS)
!       end do
!       do IV = 1, nAsh(iSymY)
!         IVABS = IV + NAES(iSymY)
!         VALUE = VALUE - Hbar(IYABS,IVABS)*G2(ITABS,IXABS,IUABS,IVABS)
!       end do
        ! RDM3 term
        ! RDM2 term
        DO IV=1,NASHT
          IVABS=IV!+NAES(ISYM)
          DO IW=1,NASHT
            IWABS=IW!+NAES(ISYM)
            DO IZ=1,NASHT
              IZABS=IZ!+NAES(ISYM)
              if (IVABS == IZABS) VALUE = VALUE - Half*Gact(IVABS,IWABS,IXABS,IZABS)*G2(ITABS,IWABS,IUABS,IYABS)
              if (IYABS == IVABS) VALUE = VALUE - Half*Gact(IVABS,IWABS,IXABS,IZABS)*G2(ITABS,IZABS,IUABS,IWABS)
              if (IXABS == IVABS) VALUE = VALUE - Half*Gact(IVABS,IWABS,IYABS,IZABS)*G2(ITABS,IWABS,IUABS,IZABS)
              if (IVABS == IZABS) VALUE = VALUE - Half*Gact(IVABS,IWABS,IYABS,IZABS)*G2(ITABS,IXABS,IUABS,IWABS)
            END DO
          END DO
        END DO
        BF(IBADR)=VALUE
!       write (*,'(4i3,f20.10)') itabs,iuabs,ixabs,iyabs,value
      END DO
    END DO
    CALL Free_Tsk(ID)

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

      G3VAL = G3(iG3)

      ! F(tuvxyz) -> BC(vut,xyz)
      if (MUL(iST,iSV) == iSym) call Add_MKBNEVF(iT,iU,iV,iX,iY,iZ)

      if (iTU /= iVX .or. iVX /= iYZ) then
        if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
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

      if (iT == iU .and. iV == iX .and. iY == iZ) cycle
      if (iT == iU .and. iV == iZ .and. iX == iY) cycle
      if (iX == iV .and. iT == iZ .and. iU == iY) cycle
      if (iZ == iY .and. iV == iU .and. iX == iT) cycle

      ! F(utxvzy) -> BC(xtu,vzy)
      if (MUL(iSU,iSX) == iSym) call Add_MKBNEVF(iU,iT,iX,iV,iZ,iY)

      if (iTU == iVX .and. iVX == iYZ) cycle
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
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

    NASP=NTGEU(ISYM)
    NBFP=(NASP*(NASP+1))/2
    IF(NBFP > 0) CALL mma_allocate(BFP,NBFP,Label='BFP')
    NASM=NTGTU(ISYM)
    NBFM=(NASM*(NASM+1))/2
    IF(NBFM > 0) CALL mma_allocate(BFM,NBFM,Label='BFM')

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
        BTUXY=BF(IBADR)
        IF(ITU >= IYX) THEN
          IBADR=(ITU*(ITU-1))/2+IYX
        ELSE
          IBADR=(IYX*(IYX-1))/2+ITU
        END IF
        BTUYX=BF(IBADR)
        IBPADR=(ITGEU*(ITGEU-1))/2+IXGEY
        BFP(IBPADR)=BTUXY+BTUYX
        IF(ITABS == IUABS) CYCLE
        IF(IXABS == IYABS) CYCLE
        ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
        IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
        IBMADR=(ITGTU*(ITGTU-1))/2+IXGTY
        BFM(IBMADR)=BTUXY-BTUYX
      END DO
    END DO

#ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      call GADGOP(BFP,NBFP,'+')
      call GADGOP(BFM,NBFM,'+')
    end if
#endif

    IF(NBF > 0) CALL mma_deallocate(BF)

    IF(NBFP > 0.and.NINDEP(ISYM,8) > 0) THEN
      IDISK=IDBMAT(ISYM,8)
      CALL DDAFILE(LUSBT,1,BFP,NBFP,IDISK)
      CALL mma_deallocate(BFP)
    END IF
    IF(NBFM > 0) THEN
     IF(NINDEP(ISYM,9) > 0) THEN
      IDISK=IDBMAT(ISYM,9)
      CALL DDAFILE(LUSBT,1,BFM,NBFM,IDISK)
     END IF
     CALL mma_deallocate(BFM)
    END IF
  end do

contains

  subroutine Add_MKBNEVF(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_)

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IVABS_, IYABS_, IUABS_, IXABS_, IZABS_
  integer(kind=iwp) :: IBADR, ITU_, IWABS, IXY_

  iTU_ = KTU(ITABS_,IUABS_)-NTUES(ISYM)
  do IWABS = 1, nAshT
    if (MUL(IASYM(IWABS),IASYM(IYABS_)) == iSym) then
      iXY_ = KTU(IWABS,IYABS_)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BF(IBADR) = BF(IBADR) - G3VAL*Gact(IWABS,IXABS_,IVABS_,IZABS_)
      end if
    end if
    if (MUL(IASYM(IXABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IXABS_,IWABS)-NTUES(ISYM)
      if (iTU_ >= iXY_) then
        IBADR = iTU_*(iTU_-1)/2 + iXY_
        BF(IBADR) = BF(IBADR) - G3VAL*Gact(IWABS,IYABS_,IVABS_,IZABS_)
      end if
    end if
  end do

  end subroutine Add_MKBNEVF

end subroutine MKBNEVF
