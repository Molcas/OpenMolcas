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

module BDerNEV

  use caspt2_global, only: LUSOLV
  use Constants, only: Zero, Half, Two
  use Definitions, only: iwp,wp,byte,u6
  use stdalloc, only: mma_allocate, mma_deallocate

  use NEVPT2_mod, only: nAshT

  implicit none

  real(kind=wp), allocatable :: Hact(:,:), Gact(:,:,:,:), Hbar(:,:), Htilde(:,:)
  real(kind=wp), allocatable :: Hder(:,:), Gder(:,:,:,:), Hbder(:,:), Htder(:,:)
  real(kind=wp), allocatable :: BDERA_save(:), BDERC_save(:)

contains

subroutine BDerNEV_initial()

  use caspt2_module, only: NSYM, NTUV, nAshT_ => nAshT
#ifdef _MOLCAS_MPP_
  use NEVPT2_E4, only: MAXBUF
  use caspt2_module, only: MXCI
#endif

  implicit none

  integer(kind=iwp) :: iSym, ntuvtot

  nAshT = nAshT_ ! define the number of the active orbitals used in NEVPT2_MOD
#ifdef _MOLCAS_MPP_
  MAXBUF = 2000000000/(MXCI*8)
#endif

  if (nAshT /= 0) then
    call mma_allocate(Hact,nAshT,nAshT,Label='Hact')
    call mma_allocate(Gact,nAshT,nAshT,nAshT,nAshT,Label='Gact')
    call mma_allocate(Hbar,nAshT,nAshT,Label='Hbar')
    call mma_allocate(Htilde,nAshT,nAshT,Label='Htilde')

    call nevint(nAshT,Hact,Gact,Hbar,Htilde)

    call mma_allocate(Hder,nAshT,nAshT,Label='Hder')
    call mma_allocate(Gder,nAshT,nAshT,nAshT,nAshT,Label='Gder')
    call mma_allocate(Hbder,nAshT,nAshT,Label='Hbder')
    call mma_allocate(Htder,nAshT,nAshT,Label='Htder')
    Hder(:,:) = Zero
    Gder(:,:,:,:) = Zero
    Hbder(:,:) = Zero
    Htder(:,:) = Zero

    ntuvtot = 0
    do iSym = 1, nSym
      ntuvtot = ntuvtot + NTUV(iSym)**2
    end do
    call mma_allocate(BDERA_save,ntuvtot,Label='BDERA_save')
    call mma_allocate(BDERC_save,ntuvtot,Label='BDERC_save')
    BDERA_save(:) = Zero
    BDERC_save(:) = Zero
  end if

  return

end subroutine BDerNEV_initial

subroutine BDerNEV_final1(NDPT2C,DPT2C)

  use caspt2_module, only: NACTEL
#ifdef _MOLCAS_MPP_
  use Para_Info, only: Is_Real_Par
#endif

  implicit none

  integer(kind=iwp), intent(in) :: NDPT2C
  real(kind=wp), intent(inout) :: DPT2C(NDPT2C)

  integer(kind=iwp) :: ii, jj, kk

  call mma_deallocate(Hact,safe='*')
  call mma_deallocate(Gact,safe='*')
  call mma_deallocate(Hbar,safe='*')
  call mma_deallocate(Htilde,safe='*')

  do ii = 1, NASHT
    do jj = 1, NASHT
      do kk = 1, NASHT
        Gder(kk,jj,kk,ii) = Gder(kk,jj,kk,ii) - Half*Hbder(jj,ii)
        Gder(kk,jj,kk,ii) = Gder(kk,jj,kk,ii) - Htder(jj,ii)
      end do
    end do
  end do
  Hder(1:NASHT,1:NASHT) = Hder(1:NASHT,1:NASHT) + Hbder(1:NASHT,1:NASHT) + Htder(1:NASHT,1:NASHT)
  !! Compensate the division by the number of electrons in OLagNS2 or OLagNS_RI
  Hder(1:NASHT,1:NASHT) = Hder(1:NASHT,1:NASHT)*dble(max(1,nactel))

#ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    call GADGOP(Hder,NASHT**2,'+')
    call GADGOP(Gder,NASHT**4,'+')
  end if
#endif

  call AddDEPSA(NDPT2C,nAshT,DPT2C,Hder)

  call mma_deallocate(Hder,safe='*')
  call mma_deallocate(Hbder,safe='*')
  call mma_deallocate(Htder,safe='*')

  call mma_deallocate(BDERA_save,safe='*')
  call mma_deallocate(BDERC_save,safe='*')

  return

end subroutine BDerNEV_final1

subroutine BDerNEV_final2()
  implicit none
  call mma_deallocate(Gder,safe='*')
  return
end subroutine BDerNEV_final2

subroutine BDerNEV_E4(NCONF,NLEV,CLag)

  use caspt2_module, only: NSYM, NTUV, NTUVES, NINDEP
  use Definitions, only: iwp,wp

  implicit none

  integer(kind=iwp), intent(in) :: NCONF, NLEV
  real(kind=wp), intent(inout) :: CLag(NCONF)

  integer(kind=iwp) :: ICASE1, ICASE4, NASA, NASC, NBA, NBC, NINA, NINC, ISYM

  ICASE1=1
  ICASE4=4
  do ISYM=1,NSYM
    NINA=NINDEP(ISYM,ICASE1)
    NINC=NINDEP(ISYM,ICASE4)
    if (NINA == 0 .and. NINC == 0) cycle
    NASA=NTUV(ISYM)
    NASC=NTUV(ISYM)
    NBA=(NASA*(NASA+1))/2
    NBC=(NASC*(NASC+1))/2
    if (NBA <= 0 .and. NBC <= 0) cycle
    call DERE4x(NLEV,ISYM,NASA,NASC,NCONF,BDERA_save(1+NTUVES(iSym)),BDERC_save(1+NTUVES(iSym)),CLag)
  end do

end subroutine BDerNEV_E4

subroutine BDNA(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)

  use caspt2_module, only: NTUVES, IASYM
  use NEVPT2_mod, only: derex2, derex3, ex2, ex3

  implicit none

  integer(kind=iwp), intent(in) :: iSym, NAS, NG3
  real(kind=wp), intent(in) :: BDER(NAS,NAS), SDER(NAS,NAS), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DG3(NG3)

  integer(kind=iwp) :: it, itu, iu
  integer(kind=iwp) :: iLUID, iG3, IV, IVX, IX, IY, IYZ, IZ
! integer(kind=iwp) :: IST, ISU, ISV, ISX, ISY, ISZ
  real(kind=wp) :: DG3VAL, G3VAL, DG2VAL

  integer(kind=byte), allocatable :: idxG3(:,:)

  !! E4 contributions are evaluate later
  BDERA_save(1+NTUVES(iSym):NAS*NAS+NTUVES(iSym)) = reshape(BDER(1:NAS,1:NAS), (/NAS*NAS/))

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID=0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  do IG3 = 1, NG3
    iT=idxG3(1,iG3)
    iU=idxG3(2,iG3)
    iV=idxG3(3,iG3)
    iX=idxG3(4,iG3)
    iY=idxG3(5,iG3)
    iZ=idxG3(6,iG3)
!   iST=IASYM(iT)
!   iSU=IASYM(iU)
!   iSV=IASYM(iV)
!   iSX=IASYM(iX)
!   iSY=IASYM(iY)
!   iSZ=IASYM(iZ)
!   if (ituvs /= ixyzs) cycle
    iTU=iT+NASHT*(iU-1)
    iVX=iV+NASHT*(iX-1)
    iYZ=iY+NASHT*(iZ-1)

    ! D(tuvxyz)
    call Add_MKBNEVA(iT,iU,iV,iX,iY,iZ,iG3)

    if (iTU /= iVX .or. iVX /= iYZ) then
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
        ! D(vxtuyz)
        call Add_MKBNEVA(iV,iX,iT,iU,iY,iZ,iG3)
        ! D(yzvxtu)
        call Add_MKBNEVA(iY,iZ,iV,iX,iT,iU,iG3)
        ! D(tuyzvx)
        call Add_MKBNEVA(iT,iU,iY,iZ,iV,iX,iG3)
      end if
      ! D(yztuvx)
      call Add_MKBNEVA(iY,iZ,iT,iU,iV,iX,iG3)
      ! D(vxyztu)
      call Add_MKBNEVA(iV,iX,iY,iZ,iT,iU,iG3)
    end if

    if (iT == iU .and. iV == iX .and. iY == iZ) cycle
    if (iT == iU .and. iV == iZ .and. iX == iY) cycle
    if (iX == iV .and. iT == iZ .and. iU == iY) cycle
    if (iZ == iY .and. iV == iU .and. iX == iT) cycle

    ! D(utxvzy)
    call Add_MKBNEVA(iU,iT,iX,iV,iZ,iY,iG3)

    if (iTU == iVX .and. iVX == iYZ) cycle
    if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
      ! D(xvutzy)
      call Add_MKBNEVA(iX,iV,iU,iT,iZ,iY,iG3)
      ! D(zyxvut)
      call Add_MKBNEVA(iZ,iY,iX,iV,iU,iT,iG3)
      ! D(utzyxv)
      call Add_MKBNEVA(iU,iT,iZ,iY,iX,iV,iG3)
    end if
    ! D(zyutxv)
    call Add_MKBNEVA(iZ,iY,iU,iT,iX,iV,iG3)
    ! D(xvzyut)
    call Add_MKBNEVA(iX,iV,iZ,iY,iU,iT,iG3)
  end do

  call mma_deallocate(idxG3)

contains

  subroutine Add_MKBNEVA(ITABS_,IUABS_,IVABS_,IXABS_,IYABS_,IZABS_,iG3_)

  use SUPERINDEX, only: KTUV
  use Symmetry_Info, only: Mul

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IUABS_, IVABS_, IXABS_, IYABS_, IZABS_, iG3_
  integer(kind=iwp) :: iaabs, ibabs, isab, ituv_, ixyz_, iwabs

  G3VAL = ex3(itabs_,iuabs_,ivabs_,ixabs_,iyabs_,izabs_,G1,G2,G3(iG3_))
  DG3VAL = Zero
  DG2VAL = Zero

  ituv_ = KTUV(ixabs_,iuabs_,itabs_) - NTUVES(iSym)
  if (Mul(IASYM(iXABS_),Mul(IASYM(IUABS_),IASYM(ITABS_))) == iSym) then
    !! S derivative
    ixyz_ = KTUV(ivabs_,iyabs_,izabs_) - NTUVES(iSym)
    DG3VAL = DG3VAL - SDER(ituv_,ixyz_)
    if (ivabs_ == ixabs_) DG2VAL = DG2VAL + SDER(ituv_,ixyz_)
  end if

  do IAABS = 1, nAshT
    do IBABS = 1, nAshT
      ! (1.2.11) ERI(2)
      ! A: A_{aut,byz} <-- +2(vx|ab) * Etu Evx Eyz
      if (Mul(IASYM(IAABS),Mul(IASYM(IUABS_),IASYM(ITABS_))) == iSym .and. &
          Mul(IASYM(IBABS),Mul(IASYM(IYABS_),IASYM(IZABS_))) == iSym) then
        iTUV_ = KTUV(IAABS,IUABS_,ITABS_)-NTUVES(ISYM)
        iXYZ_ = KTUV(IBABS,IYABS_,IZABS_)-NTUVES(ISYM)
        Gder(iaabs,ibabs,ivabs_,ixabs_) = Gder(iaabs,ibabs,ivabs_,ixabs_) + Two*BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL + Two*BDER(ituv_,ixyz_)*Gact(iaabs,ibabs,ivabs_,ixabs_)
      end if
    end do
  end do

  iTUV_ = KTUV(IXABS_,IUABS_,ITABS_)-NTUVES(ISYM)
  if (Mul(IASYM(IXABS_),Mul(IASYM(IUABS_),IASYM(ITABS_))) == iSym) then
    do IAABS = 1, nAshT
      do IBABS = 1, nAshT
        ISAB = Mul(IASYM(IAABS),IASYM(IBABS))
        ! (1.2.11) ERI(3)
        ! A: A_{xut,avb} <-- +(ab|yz) * Etu Evx Eyz
        if (Mul(IASYM(IVABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IVABS_,IBABS)-NTUVES(ISYM)
          Gder(iaabs,ibabs,iyabs_,izabs_) = Gder(iaabs,ibabs,iyabs_,izabs_) + BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL + BDER(ituv_,ixyz_)*Gact(iaabs,ibabs,iyabs_,izabs_)
        end if
        ! (1.2.11) ERI(4)
        ! A: A_{xut,abz} <-- -(va|yb) * Etu Evx Eyz
        if (Mul(IASYM(IZABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IBABS,IZABS_)-NTUVES(ISYM)
          Gder(iaabs,ivabs_,iyabs_,ibabs) = Gder(iaabs,ivabs_,iyabs_,ibabs) - BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Gact(iaabs,ivabs_,iyabs_,ibabs)
        end if
        ! (1.2.11) ERI(5)
        ! A: A_{xut,ayb} <-- +(va|bz) * Etu Evx Eyz
        if (Mul(IASYM(IYABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IYABS_,IBABS)-NTUVES(ISYM)
          Gder(iaabs,ivabs_,izabs_,ibabs) = Gder(iaabs,ivabs_,izabs_,ibabs) + BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL + BDER(ituv_,ixyz_)*Gact(iaabs,ivabs_,izabs_,ibabs)
        end if
      end do
    end do

    G3VAL = -G3VAL
    if (IVABS_ == IXABS_) G3VAL = Two*ex2(ITABS_,IUABS_,IYABS_,IZABS_,G1,G2) + G3VAL

    dg2val = dg2val + dg3val
    dg3val = -dg3val

    do IWABS = 1, nAshT
      ! A: A(xut,vwz) <-- +h_{yw} EtuEvxEyz
      if (Mul(IASYM(IWABS),Mul(IASYM(IVABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IVABS_,IWABS,IZABS_)-NTUVES(ISYM)
        Hder(iyabs_,iwabs) = Hder(iyabs_,iwabs) + BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL + BDER(ituv_,ixyz_)*Hact(iyabs_,iwabs)
      end if
      ! A: A(xut,vyw) <-- -h_{wz} EtuEvxEyz
      if (Mul(IASYM(IWABS),Mul(IASYM(IVABS_),IASYM(IYABS_))) == iSym) then
        iXYZ_ = KTUV(IVABS_,IYABS_,IWABS)-NTUVES(ISYM)
        Htder(iwabs,izabs_) = Htder(iwabs,izabs_) - BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Htilde(iwabs,izabs_)
      end if
      ! A: A(xut,wyz) <-- +h_{wv} EtuEvxEyz
      if (Mul(IASYM(IWABS),Mul(IASYM(IYABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IWABS,IYABS_,IZABS_)-NTUVES(ISYM)
        Hder(iwabs,ivabs_) = Hder(iwabs,ivabs_) + BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL + BDER(ituv_,ixyz_)*Hact(iwabs,ivabs_)
      end if
    end do

    do IAABS = 1, nAshT
      do IBABS = 1, nAshT
        ISAB = Mul(IASYM(IAABS),IASYM(IBABS))
        ! (1.2.11) ERI(1)
        ! A: A_{xut,vab} <-- +(ya|zb) * Etu Evx Eyz
        if (Mul(IASYM(IVABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IVABS_,IAABS,IBABS)-NTUVES(ISYM)
          Gder(iaabs,iyabs_,ibabs,izabs_) = Gder(iaabs,iyabs_,ibabs,izabs_) - BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Gact(iaabs,iyabs_,ibabs,izabs_)
        end if
      end do
    end do
  end if

  dg3val = -dg3val

  call derex3(itabs_,iuabs_,ivabs_,ixabs_,iyabs_,izabs_,DG1,DG2,DG3(iG3_),DG3VAL)
  if (ivabs_ == ixabs_) then
    dg2val = -Two*(dg3val - dg2val)
    call derex2(itabs_,iuabs_,iyabs_,izabs_,DG1,DG2,DG2VAL)
  end if

  end subroutine Add_MKBNEVA

end subroutine BDNA

subroutine BDNB(iSym,iCase,NAS,NG3,BDER0,SDER0,G1,G2,G3,DG1,DG2,DG3)

  use caspt2_module, only: NTU, NTUES, NTGEUES, NTGTUES, IASYM
  use NEVPT2_mod, only: tG2, tG3, dertG2, dertG3
  use SUPERINDEX, only: KTU, MTGEU, MTGTU, MTU
  use Symmetry_Info, only: Mul
  use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk

  implicit none

  integer(kind=iwp), intent(in) :: iSym, iCase, NAS, NG3
  real(kind=wp), intent(in) :: BDER0(NAS,NAS), SDER0(NAS,NAS), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DG3(NG3)

  integer(kind=iwp) :: it, itabs, itu, ituabs, itu0, iu, iuabs, ivabs, iwabs, ixabs, ixy, ixyabs, ixy0, iyabs, iyx, izabs
  integer(kind=iwp) :: iLUID, iG3, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IV, IVX, IX, IXYZS, IY, IYZ, IZ
  integer(kind=iwp) :: iTgeUabs, iTgtUabs, iXgeYabs, iXgtYabs
  real(kind=wp) :: BDERval, DG3VAL, G3VAL, SDERval

  real(kind=wp), allocatable :: BDER(:,:), SDER(:,:)
  integer(kind=byte), allocatable :: idxG3(:,:)
  integer(kind=iwp) :: ID, nTask

  call mma_allocate(BDER,NTU(iSym),NTU(iSym),Label='BDER_B')
  call mma_allocate(SDER,NTU(iSym),NTU(iSym),Label='SDER_B')

  BDER(:,:) = Zero
  SDER(:,:) = Zero

  do itu0 = 1, NAS
    if (iCase == 2) then
      itgeuabs = itu0 + NTGEUES(iSym)
      itabs    = MTGEU(1,itgeuabs)
      iuabs    = MTGEU(2,itgeuabs)
    else if (iCase == 3) then
      itgtuabs = itu0 + NTGTUES(iSym)
      itabs    = MTGTU(1,itgtuabs)
      iuabs    = MTGTU(2,itgtuabs)
    end if
    itu = KTU(itabs,iuabs) - NTUES(iSym)
    do ixy0 = 1, NAS
      if (iCase == 2) then
        ixgeyabs = ixy0 + NTGEUES(iSym)
        ixabs    = MTGEU(1,ixgeyabs)
        iyabs    = MTGEU(2,ixgeyabs)
      else if (iCase == 3) then
        ixgtyabs = ixy0 + NTGTUES(iSym)
        ixabs    = MTGTU(1,ixgtyabs)
        iyabs    = MTGTU(2,ixgtyabs)
      end if
      ixy = KTU(ixabs,iyabs) - NTUES(iSym)
      iyx = KTU(iyabs,ixabs) - NTUES(iSym)
      if (iCase == 2) then
        BDER(itu,ixy) = BDER(itu,ixy) + BDER0(itu0,ixy0)
        BDER(itu,iyx) = BDER(itu,iyx) + BDER0(itu0,ixy0)
        SDER(itu,ixy) = SDER(itu,ixy) + SDER0(itu0,ixy0)
        SDER(itu,iyx) = SDER(itu,iyx) + SDER0(itu0,ixy0)
      end if
      if (iCase == 3) then
        if (itabs == iuabs) cycle
        if (ixabs == iyabs) cycle
        BDER(itu,ixy) = BDER(itu,ixy) + BDER0(itu0,ixy0)
        BDER(itu,iyx) = BDER(itu,iyx) - BDER0(itu0,ixy0)
        SDER(itu,ixy) = SDER(itu,ixy) + SDER0(itu0,ixy0)
        SDER(itu,iyx) = SDER(itu,iyx) - SDER0(itu0,ixy0)
      end if
    end do
  end do

  BDER(:,:) = Two*BDER(:,:)
  SDER(:,:) = Two*SDER(:,:)

  nTask = NTU(iSym)
  call Init_Tsk(ID,nTask)

  do while (Rsv_Tsk(ID,itu))
    ituabs = itu + NTUES(iSym)
    itabs  = MTU(1,ituabs)
    iuabs  = MTU(2,ituabs)
    do ixy = 1, NTU(iSym)
      ixyabs = ixy + NTUES(iSym)
      ixabs  = MTU(1,ixyabs)
      iyabs  = MTU(2,ixyabs)

      BDERval = BDER(ITU,IXY)

      !! B derivative: Hbar contributions
      do ivabs = 1, NASHT
        Hbder(ivabs,ixabs) = Hbder(ivabs,ixabs) + BDERval*tG2(itabs,iuabs,ivabs,iyabs,G1,G2)
        Hbder(ivabs,iyabs) = Hbder(ivabs,iyabs) + BDERval*tG2(itabs,iuabs,ixabs,ivabs,G1,G2)
        call dertG2(itabs,iuabs,ivabs,iyabs,DG1,DG2,+BDERval*Hbar(ixabs,ivabs))
        call dertG2(itabs,iuabs,ixabs,ivabs,DG1,DG2,+BDERval*Hbar(iyabs,ivabs))
      end do

      !! B derivative: ERI contributions
      do ivabs = 1, NASHT
        do iwabs = 1, NASHT
          do izabs = 1, NASHT
            if (ivabs == izabs) then
              Gder(ivabs,izabs,iwabs,ixabs) = Gder(ivabs,izabs,iwabs,ixabs) + Two*BDERval*tG2(itabs,iuabs,iwabs,iyabs,G1,G2)
              call dertG2(itabs,iuabs,iwabs,iyabs,DG1,DG2,Two*BDERval*Gact(ivabs,izabs,iwabs,ixabs))
            end if
            if (iwabs == izabs) then
              Gder(ivabs,izabs,iwabs,ixabs) = Gder(ivabs,izabs,iwabs,ixabs) - Half*BDERval*tG2(itabs,iuabs,ivabs,iyabs,G1,G2)
              call dertG2(itabs,iuabs,ivabs,iyabs,DG1,DG2,-Half*BDERval*Gact(ivabs,izabs,iwabs,ixabs))
            end if
            if (iyabs == izabs) then
              Gder(ivabs,izabs,iwabs,ixabs) = Gder(ivabs,izabs,iwabs,ixabs) - Half*BDERval*tG2(itabs,iuabs,iwabs,ivabs,G1,G2)
              call dertG2(itabs,iuabs,iwabs,ivabs,DG1,DG2,-Half*BDERval*Gact(ivabs,izabs,iwabs,ixabs))
            end if

            if (ivabs == izabs) then
              Gder(ivabs,izabs,iwabs,iyabs) = Gder(ivabs,izabs,iwabs,iyabs) + Two*BDERval*tG2(itabs,iuabs,ixabs,iwabs,G1,G2)
              call dertG2(itabs,iuabs,ixabs,iwabs,DG1,DG2,Two*BDERval*Gact(ivabs,izabs,iwabs,iyabs))
            end if
            if (ixabs == izabs) then
              Gder(ivabs,izabs,iwabs,iyabs) = Gder(ivabs,izabs,iwabs,iyabs) - Half*BDERval*tG2(itabs,iuabs,ivabs,iwabs,G1,G2)
              call dertG2(itabs,iuabs,ivabs,iwabs,DG1,DG2,-Half*BDERval*Gact(ivabs,izabs,iwabs,iyabs))
            end if
            if (iwabs == izabs) then
              Gder(ivabs,izabs,iwabs,iyabs) = Gder(ivabs,izabs,iwabs,iyabs) - Half*BDERval*tG2(itabs,iuabs,ixabs,ivabs,G1,G2)
              call dertG2(itabs,iuabs,ixabs,ivabs,DG1,DG2,-Half*BDERval*Gact(ivabs,izabs,iwabs,iyabs))
            end if
          end do
        end do
      end do

      SDERval = SDER(ITU,IXY)

      !! S derivative
      call dertG2(itabs,iuabs,ixabs,iyabs,DG1,DG2,SDERval)
    end do
  end do

  call Free_Tsk(ID)

  call mma_deallocate(SDER)

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID=0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  !! The remaining D3 contributions
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
    ituvs=Mul(IST,Mul(ISU,ISV))
    ixyzs=Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU=iT+NASHT*(iU-1)
    iVX=iV+NASHT*(iX-1)
    iYZ=iY+NASHT*(iZ-1)

    G3VAL = tG3(IT,IU,IV,IX,IY,IZ,G1,G2,G3(iG3))
    DG3VAL = Zero

    ! F(tuvxyz) -> BC(vut,xyz)
    if (Mul(iST,iSV) == iSym) call Add_MKBNEVB(iT,iU,iV,iX,iY,iZ,G3VAL,DG3VAL)

    if (iTU /= iVX .or. iVX /= iYZ) then
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
        ! F(vxtuyz) -> BC(txv,uyz)
        if (Mul(iSV,iST) == iSym) call Add_MKBNEVB(iV,iX,iT,iU,iY,iZ,G3VAL,DG3VAL)
        ! F(yzvxtu) -> BC(vzy,xtu)
        if (Mul(iSY,iSV) == iSym) call Add_MKBNEVB(iY,iZ,iV,iX,iT,iU,G3VAL,DG3VAL)
        ! F(tuyzvx) -> BC(yut,zvx)
        if (Mul(iST,iSY) == iSym) call Add_MKBNEVB(iT,iU,iY,iZ,iV,iX,G3VAL,DG3VAL)
      end if
      ! F(yztuvx) -> BC(tzy,uvx)
      if (Mul(iSY,iST) == iSym) call Add_MKBNEVB(iY,iZ,iT,iU,iV,iX,G3VAL,DG3VAL)
      ! F(vxyztu) -> BC(yxv,ztu)
      if (Mul(iSV,iSY) == iSym) call Add_MKBNEVB(iV,iX,iY,iZ,iT,iU,G3VAL,DG3VAL)
    end if

    if ((iT /= iU .or. iV /= iX .or. iY /= iZ) .and. (iT /= iU .or. iV /= iZ .or. iX /= iY) .and. &
        (iX /= iV .or. iT /= iZ .or. iU /= iY) .and. (iZ /= iY .or. iV /= iU .or. iX /= iT)) then
      ! F(utxvzy) -> BC(xtu,vzy)
      if (Mul(iSU,iSX) == iSym) call Add_MKBNEVB(iU,iT,iX,iV,iZ,iY,G3VAL,DG3VAL)

      if (iTU /= iVX .or. iVX /= iYZ) then
        if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
          ! F(xvutzy) -> BC(uvx,tzy)
          if (Mul(iSX,iSU) == iSym) call Add_MKBNEVB(iX,iV,iU,iT,iZ,iY,G3VAL,DG3VAL)
          ! F(zyxvut) -> BC(xyz,vut)
          if (Mul(iSZ,iSX) == iSym) call Add_MKBNEVB(iZ,iY,iX,iV,iU,iT,G3VAL,DG3VAL)
          ! F(utzyxv) -> BC(ztu,yxv)
          if (Mul(iSU,iSZ) == iSym) call Add_MKBNEVB(iU,iT,iZ,iY,iX,iV,G3VAL,DG3VAL)
        end if
        ! F(zyutxv) -> BC(uyz,txv)
        if (Mul(iSZ,iSU) == iSym) call Add_MKBNEVB(iZ,iY,iU,iT,iX,iV,G3VAL,DG3VAL)
        ! F(xvzyut) -> BC(zvx,yut)
        if (Mul(iSX,iSZ) == iSym) call Add_MKBNEVB(iX,iV,iZ,iY,iU,iT,G3VAL,DG3VAL)
      end if
    end if

    call dertG3(it,iu,iv,ix,iy,iz,DG1,DG2,DG3(iG3),DG3VAL)
  end do

  call mma_deallocate(idxG3)

  call mma_deallocate(BDER)

  return

contains

  subroutine Add_MKBNEVB(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_,g3val_,dg3val_)

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IVABS_, IYABS_, IUABS_, IXABS_, IZABS_
  real(kind=wp), intent(in) :: g3val_
  real(kind=wp), intent(inout) :: dg3val_
  integer(kind=iwp) :: ITU_, IWABS, IXY_

  iTU_ = KTU(ITABS_,IUABS_)-NTUES(ISYM)
  do IWABS = 1, nAshT
    if (MUL(IASYM(IWABS),IASYM(IYABS_)) == iSym) then
      iXY_ = KTU(IWABS,IYABS_)-NTUES(ISYM)
      dg3val_ = dg3val_ - BDER(iTU_,iXY_)*Gact(izabs_,ivabs_,ixabs_,iwabs)
      Gder(izabs_,ivabs_,ixabs_,iwabs) = Gder(izabs_,ivabs_,ixabs_,iwabs) - BDER(iTU_,iXY_)*g3val_
    end if
    if (MUL(IASYM(IXABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IXABS_,IWABS)-NTUES(ISYM)
      dg3val_ = dg3val_ - BDER(iTU_,iXY_)*Gact(izabs_,ivabs_,iyabs_,iwabs)
      Gder(izabs_,ivabs_,iyabs_,iwabs) = Gder(izabs_,ivabs_,iyabs_,iwabs) - BDER(iTU_,iXY_)*g3val_
    end if
  end do

  end subroutine Add_MKBNEVB

end subroutine BDNB

subroutine BDNC(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)

  use caspt2_module, only: NTUVES, IASYM
  use NEVPT2_mod, only: derex3, ex3

  implicit none

  integer(kind=iwp), intent(in) :: iSym, NAS, NG3
  real(kind=wp), intent(in) :: BDER(NAS,NAS), SDER(NAS,NAS), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DG3(NG3)

  integer(kind=iwp) :: it, itu, iu
  integer(kind=iwp) :: iLUID, iG3, IV, IVX, IX, IY, IYZ, IZ
! integer(kind=iwp) :: IST, ISU, ISV, ISX, ISY, ISZ
  real(kind=wp) :: DG3VAL, G3VAL

  integer(kind=byte), allocatable :: idxG3(:,:)

  !! E4 contributions are evaluate later
  BDERC_save(1+NTUVES(iSym):NAS*NAS+NTUVES(iSym)) = reshape(BDER(1:NAS,1:NAS), (/NAS*NAS/))

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID=0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  do IG3 = 1, NG3
    iT=idxG3(1,iG3)
    iU=idxG3(2,iG3)
    iV=idxG3(3,iG3)
    iX=idxG3(4,iG3)
    iY=idxG3(5,iG3)
    iZ=idxG3(6,iG3)
!   iST=IASYM(iT)
!   iSU=IASYM(iU)
!   iSV=IASYM(iV)
!   iSX=IASYM(iX)
!   iSY=IASYM(iY)
!   iSZ=IASYM(iZ)
!   if (ituvs /= ixyzs) cycle
    iTU=iT+NASHT*(iU-1)
    iVX=iV+NASHT*(iX-1)
    iYZ=iY+NASHT*(iZ-1)

    ! D(tuvxyz)
    call Add_MKBNEVC(iT,iU,iV,iX,iY,iZ,iG3)

    if (iTU /= iVX .or. iVX /= iYZ) then
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
        ! D(vxtuyz)
        call Add_MKBNEVC(iV,iX,iT,iU,iY,iZ,iG3)
        ! D(yzvxtu)
        call Add_MKBNEVC(iY,iZ,iV,iX,iT,iU,iG3)
        ! D(tuyzvx)
        call Add_MKBNEVC(iT,iU,iY,iZ,iV,iX,iG3)
      end if
      ! D(yztuvx)
      call Add_MKBNEVC(iY,iZ,iT,iU,iV,iX,iG3)
      ! D(vxyztu)
      call Add_MKBNEVC(iV,iX,iY,iZ,iT,iU,iG3)
    end if

    if (iT == iU .and. iV == iX .and. iY == iZ) cycle
    if (iT == iU .and. iV == iZ .and. iX == iY) cycle
    if (iX == iV .and. iT == iZ .and. iU == iY) cycle
    if (iZ == iY .and. iV == iU .and. iX == iT) cycle

    ! D(utxvzy)
    call Add_MKBNEVC(iU,iT,iX,iV,iZ,iY,iG3)

    if (iTU == iVX .and. iVX == iYZ) cycle
    if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
      ! D(xvutzy)
      call Add_MKBNEVC(iX,iV,iU,iT,iZ,iY,iG3)
      ! D(zyxvut)
      call Add_MKBNEVC(iZ,iY,iX,iV,iU,iT,iG3)
      ! D(utzyxv)
      call Add_MKBNEVC(iU,iT,iZ,iY,iX,iV,iG3)
    end if
    ! D(zyutxv)
    call Add_MKBNEVC(iZ,iY,iU,iT,iX,iV,iG3)
    ! D(xvzyut)
    call Add_MKBNEVC(iX,iV,iZ,iY,iU,iT,iG3)
  end do

  call mma_deallocate(idxG3)

contains

  subroutine Add_MKBNEVC(ITABS_,IUABS_,IVABS_,IXABS_,IYABS_,IZABS_,iG3_)

  use SUPERINDEX, only: KTUV
  use Symmetry_Info, only: Mul

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IUABS_, IVABS_, IXABS_, IYABS_, IZABS_, iG3_
  integer(kind=iwp) :: iaabs, ibabs, isab, ituv_, ixyz_, iwabs

  G3VAL = ex3(itabs_,iuabs_,ivabs_,ixabs_,iyabs_,izabs_,G1,G2,G3(iG3_))
  DG3VAL = Zero

  ituv_ = KTUV(ivabs_,iuabs_,itabs_) - NTUVES(iSym)
  if (Mul(IASYM(iVABS_),Mul(IASYM(IUABS_),IASYM(ITABS_))) == iSym) then
    !! S derivative
    ixyz_ = KTUV(ixabs_,iyabs_,izabs_) - NTUVES(iSym)
    DG3VAL = DG3VAL + SDER(ituv_,ixyz_)

    !! B derivative
    do IWABS = 1, nAshT
      ! C: A(vut,xwz) <-- +h_{yw} Etu*Evx*Eyz
      if (Mul(IASYM(IWABS),Mul(IASYM(IXABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IXABS_,IWABS,IZABS_)-NTUVES(ISYM)
        Hder(iyabs_,iwabs) = Hder(iyabs_,iwabs) + BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL + BDER(ituv_,ixyz_)*Hact(iyabs_,iwabs)
      end if
      ! C: A(vut,xyw) <-- -htilde_{wz} Etu*Evx*Eyz
      if (Mul(IASYM(IWABS),Mul(IASYM(IXABS_),IASYM(IYABS_))) == iSym) then
        iXYZ_ = KTUV(IXABS_,IYABS_,IWABS)-NTUVES(ISYM)
        Htder(iwabs,izabs_) = Htder(iwabs,izabs_) - BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Htilde(iwabs,izabs_)
      end if
      ! C: A(vut,wyz) <-- -htilde_{wx} Etu*Evx*Eyz
      if (Mul(IASYM(IWABS),Mul(IASYM(IYABS_),IASYM(IZABS_))) == iSym) then
        iXYZ_ = KTUV(IWABS,IYABS_,IZABS_)-NTUVES(ISYM)
        Htder(iwabs,ixabs_) = Htder(iwabs,ixabs_) - BDER(ituv_,ixyz_)*G3VAL
        DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Htilde(iwabs,ixabs_)
      end if
    end do

    do IAABS = 1, nAshT
      do IBABS = 1, nAshT
        ISAB = Mul(IASYM(IAABS),IASYM(IBABS))
        ! C: A_{vut,abx} <-- -(ab|yz) * Etu Evx Eyz
        if (Mul(IASYM(IXABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IBABS,IXABS_)-NTUVES(ISYM)
          Gder(iaabs,ibabs,iyabs_,izabs_) = Gder(iaabs,ibabs,iyabs_,izabs_) - BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Gact(iaabs,ibabs,iyabs_,izabs_)
        end if
        ! C: A_{vut,xab} <-- -(ay|bz) * Etu Evx Eyz
        if (Mul(IASYM(IXABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IXABS_,IAABS,IBABS)-NTUVES(ISYM)
          Gder(iaabs,iyabs_,ibabs,izabs_) = Gder(iaabs,iyabs_,ibabs,izabs_) - BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Gact(iaabs,iyabs_,ibabs,izabs_)
        end if
        ! C: A_{vut,abz} <-- -(ax|by) * Etu Evx Eyz
        if (Mul(IASYM(IZABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IBABS,IZABS_)-NTUVES(ISYM)
          Gder(iaabs,ixabs_,ibabs,iyabs_) = Gder(iaabs,ixabs_,ibabs,iyabs_) - BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL - BDER(ituv_,ixyz_)*Gact(iaabs,ixabs_,ibabs,iyabs_)
        end if
        ! C: A_{vut,ayb} <-- +(ax|bz) * Etu Evx Eyz
        if (Mul(IASYM(IYABS_),ISAB) == iSym) then
          iXYZ_ = KTUV(IAABS,IYABS_,IBABS)-NTUVES(ISYM)
          Gder(iaabs,ixabs_,ibabs,izabs_) = Gder(iaabs,ixabs_,ibabs,izabs_) + BDER(ituv_,ixyz_)*G3VAL
          DG3VAL = DG3VAL + BDER(ituv_,ixyz_)*Gact(iaabs,ixabs_,ibabs,izabs_)
        end if
      end do
    end do
  end if

  call derex3(itabs_,iuabs_,ivabs_,ixabs_,iyabs_,izabs_,DG1,DG2,DG3(iG3_),DG3VAL)

  end subroutine Add_MKBNEVC

end subroutine BDNC

subroutine BDND(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)

  use caspt2_module, only: NTUES, IASYM
  use nevpt2_mod, only: derex2, derex3, ex2, ex3
  use SUPERINDEX, only: MTU
  use Symmetry_Info, only: Mul
  use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk

  implicit none

  integer(kind=iwp), intent(in) :: iSym, NAS, NG3
  real(kind=wp), intent(in) :: BDER(NAS,NAS), SDER(NAS,NAS), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DG3(NG3)

  integer(kind=iwp) :: it, itabs, itu, ituabs, iu, iuabs, ivabs, iwabs, ixabs, ixy, ixyabs, iyabs, izabs
  integer(kind=iwp) :: iLUID, iG3, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IV, IVX, IX, IXYZS, IY, IYZ, IZ
  real(kind=wp) :: DG3VAL, G3VAL, tmp

  real(kind=wp), allocatable :: BDER1(:,:), BDER2(:,:), SDER1(:,:), SDER2(:,:)
  integer(kind=byte), allocatable :: idxG3(:,:)
  integer(kind=iwp) :: ID, nTask

  call mma_allocate(BDER1,NAS/2,NAS/2,Label='BDER1_D')
  call mma_allocate(BDER2,NAS/2,NAS/2,Label='BDER2_D')
  call mma_allocate(SDER1,NAS/2,NAS/2,Label='SDER1_D')
  call mma_allocate(SDER2,NAS/2,NAS/2,Label='SDER2_D')

  BDER1(1:NAS/2,1:NAS/2) = Two*BDER(1:NAS/2,1:NAS/2) - BDER(NAS/2+1:NAS,1:NAS/2) - BDER(1:NAS/2,NAS/2+1:NAS)
  BDER2(1:NAS/2,1:NAS/2) = BDER(NAS/2+1:NAS,NAS/2+1:NAS)
  SDER1(1:NAS/2,1:NAS/2) = Two*SDER(1:NAS/2,1:NAS/2) - SDER(NAS/2+1:NAS,1:NAS/2) - SDER(1:NAS/2,NAS/2+1:NAS)
  SDER2(1:NAS/2,1:NAS/2) = SDER(NAS/2+1:NAS,NAS/2+1:NAS)

  nTask = NAS/2
  call Init_Tsk(ID,nTask)

  do while (Rsv_Tsk(ID,itu))
    ituabs = itu + NTUES(iSym)
    itabs  = MTU(1,ituabs)
    iuabs  = MTU(2,ituabs)
    do ixy = 1, NAS/2
      ixyabs = ixy + NTUES(iSym)
      ixabs  = MTU(1,ixyabs)
      iyabs  = MTU(2,ixyabs)

      !! A matrix
      do ivabs = 1, NASHT
        Hbder(ivabs,ixabs) = Hbder(ivabs,ixabs) + BDER1(itu,ixy)*ex2(iuabs,itabs,ivabs,iyabs,G1,G2)
        Hbder(iyabs,ivabs) = Hbder(iyabs,ivabs) - BDER1(itu,ixy)*ex2(iuabs,itabs,ixabs,ivabs,G1,G2)
        call derex2(iuabs,itabs,ivabs,iyabs,DG1,DG2,+BDER1(itu,ixy)*Hbar(ivabs,ixabs))
        call derex2(iuabs,itabs,ixabs,ivabs,DG1,DG2,-BDER1(itu,ixy)*Hbar(iyabs,ivabs))
      end do

      !! D matrix (hbar terms)
      do ivabs = 1, NASHT
        tmp = -ex2(IXABS,ITABS,IUABS,IVABS,G1,G2)
        if (ITABS == IXABS) tmp = tmp + Two*G1(IUABS,IVABS)
        if (IUABS == ITABS) tmp = tmp + G1(IVABS,IXABS)
        Hbder(iyabs,ivabs) = Hbder(iyabs,ivabs) - BDER2(itu,ixy)*tmp
        tmp = -BDER2(itu,ixy)*Hbar(iyabs,ivabs)
        call derex2(ixabs,itabs,iuabs,ivabs,DG1,DG2,-tmp)
        if (itabs == ixabs) DG1(iuabs,ivabs) = DG1(iuabs,ivabs) + Two*tmp
        if (iuabs == itabs) DG1(ivabs,ixabs) = DG1(ivabs,ixabs) + tmp

        tmp = -ex2(IVABS,ITABS,IUABS,IYABS,G1,G2)
        if (ITABS == IVABS) tmp = tmp + Two*G1(IUABS,IYABS)
        if (IUABS == ITABS) tmp = tmp + G1(IYABS,IVABS)
        Hbder(ivabs,ixabs) = Hbder(ivabs,ixabs) + BDER2(itu,ixy)*tmp
        tmp = +BDER2(itu,ixy)*Hbar(ivabs,ixabs)
        call derex2(ivabs,itabs,iuabs,iyabs,DG1,DG2,-tmp)
        if (itabs == ivabs) DG1(iuabs,iyabs) = DG1(iuabs,iyabs) + Two*tmp
        if (iuabs == itabs) DG1(iyabs,ivabs) = DG1(iyabs,ivabs) + tmp
      end do

      !! D matrix (ERI terms)
      do ivabs = 1, NASHT
        do iwabs = 1, NASHT
          do izabs = 1, NASHT
            tmp = Zero
            if (itabs == ixabs) tmp = tmp + Two*ex2(ivabs,izabs,iuabs,iwabs,G1,G2)
            if (itabs == ixabs) tmp = tmp + Two*ex2(iuabs,iwabs,ivabs,izabs,G1,G2)
            if (iuabs == itabs) tmp = tmp + ex2(ivabs,izabs,ixabs,iwabs,G1,G2)
            if (iuabs == itabs) tmp = tmp + ex2(ixabs,iwabs,ivabs,izabs,G1,G2)
            if (ITABS == IVABS .and. IZABS == IXABS) tmp = tmp + Two*G1(IUABS,IWABS)
            if (ITABS == IVABS) tmp = tmp - ex2(IXABS,IZABS,IUABS,IWABS,G1,G2)
            if (IUABS == IZABS .and. IXABS == ITABS) tmp = tmp - Two*G1(IVABS,IWABS)
            if (IUABS == IZABS) tmp = tmp + ex2(IXABS,ITABS,IVABS,IWABS,G1,G2)
            Gder(ivabs,izabs,iyabs,iwabs) = Gder(ivabs,izabs,iyabs,iwabs) - Half*BDER2(itu,ixy)*tmp
            tmp = -Half*BDER2(itu,ixy)*Gact(ivabs,izabs,iyabs,iwabs)
            if (itabs == ixabs) call derex2(ivabs,izabs,iuabs,iwabs,DG1,DG2,Two*tmp)
            if (itabs == ixabs) call derex2(iuabs,iwabs,ivabs,izabs,DG1,DG2,Two*tmp)
            if (iuabs == itabs) call derex2(ivabs,izabs,ixabs,iwabs,DG1,DG2,tmp)
            if (iuabs == itabs) call derex2(ixabs,iwabs,ivabs,izabs,DG1,DG2,tmp)
            if (itabs == ivabs .and. izabs == ixabs) DG1(iuabs,iwabs) = DG1(iuabs,iwabs) + Two*tmp
            if (itabs == ivabs) call derex2(ixabs,izabs,iuabs,iwabs,DG1,DG2,-tmp)
            if (iuabs == izabs .and. ixabs == itabs) DG1(ivabs,iwabs) = DG1(ivabs,iwabs) - Two*tmp
            if (iuabs == izabs) call derex2(ixabs,itabs,ivabs,iwabs,DG1,DG2,+tmp)

            tmp = Zero
            if (itabs == iwabs) tmp = tmp + Two*ex2(ivabs,izabs,iuabs,iyabs,G1,G2)
            if (itabs == iwabs) tmp = tmp + Two*ex2(iuabs,iyabs,ivabs,izabs,G1,G2)
            if (iuabs == itabs) tmp = tmp + ex2(ivabs,izabs,iwabs,iyabs,G1,G2)
            if (iuabs == itabs) tmp = tmp + ex2(iwabs,iyabs,ivabs,izabs,G1,G2)
            if (itabs == ivabs .and. izabs == iwabs) tmp = tmp + Two*G1(iuabs,iyabs)
            if (itabs == ivabs) tmp = tmp - ex2(iwabs,izabs,iuabs,iyabs,G1,G2)
            if (iuabs == izabs .and. itabs == iwabs) tmp = tmp - Two*G1(ivabs,iyabs)
            if (iuabs == izabs) tmp = tmp + ex2(iwabs,itabs,ivabs,iyabs,G1,G2)
            Gder(ivabs,izabs,iwabs,ixabs) = Gder(ivabs,izabs,iwabs,ixabs) + Half*BDER2(itu,ixy)*tmp
            tmp = Half*BDER2(itu,ixy)*Gact(ivabs,izabs,iwabs,ixabs)
            if (itabs == iwabs) call derex2(ivabs,izabs,iuabs,iyabs,DG1,DG2,Two*tmp)
            if (itabs == iwabs) call derex2(iuabs,iyabs,ivabs,izabs,DG1,DG2,Two*tmp)
            if (iuabs == itabs) call derex2(ivabs,izabs,iwabs,iyabs,DG1,DG2,tmp)
            if (iuabs == itabs) call derex2(iwabs,iyabs,ivabs,izabs,DG1,DG2,tmp)
            if (itabs == ivabs .and. izabs == iwabs) DG1(iuabs,iyabs) = DG1(iuabs,iyabs) + Two*tmp
            if (itabs == ivabs) call derex2(iwabs,izabs,iuabs,iyabs,DG1,DG2,-tmp)
            if (iuabs == izabs .and. itabs == iwabs) DG1(ivabs,iyabs) = DG1(ivabs,iyabs) - Two*tmp
            if (iuabs == izabs) call derex2(iwabs,itabs,ivabs,iyabs,DG1,DG2,+tmp)
          end do
        end do
      end do

      !! Derivative of S11
      DG2(iuabs,itabs,ixabs,iyabs) = DG2(iuabs,itabs,ixabs,iyabs) + SDER1(itu,ixy)
      if (ixabs == itabs) DG1(iuabs,iyabs) = DG1(iuabs,iyabs) + SDER1(itu,ixy)
      !! Derivative of S22
      DG2(ixabs,itabs,iuabs,iyabs) = DG2(ixabs,itabs,iuabs,iyabs) - SDER2(itu,ixy)
      if (ixabs == itabs) DG1(iuabs,iyabs) = DG1(iuabs,iyabs) + Two*SDER2(itu,ixy)
    end do
  end do

  call Free_Tsk(ID)

  call mma_deallocate(SDER1)
  call mma_deallocate(SDER2)

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID=0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  !! The remaining D3 contributions
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
    ituvs=Mul(IST,Mul(ISU,ISV))
    ixyzs=Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU=iT+NASHT*(iU-1)
    iVX=iV+NASHT*(iX-1)
    iYZ=iY+NASHT*(iZ-1)

    ! D(tuvxyz)
    if (Mul(iSU,iST) == iSym .or. Mul(iSU,ISV) == iSym .or. Mul(iSX,iSY) == iSym) then
      call Add_MKBNEVD(iT,iU,iV,iX,iY,iZ,iG3)
    end if

    if (iTU /= iVX .or. iVX /= iYZ) then
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
        ! D(vxtuyz)
        if (Mul(iSX,iSV) == iSym .or. Mul(iSX,iST) == iSym .or. Mul(iSU,iSY) == iSym) then
          call Add_MKBNEVD(iV,iX,iT,iU,iY,iZ,iG3)
        end if
        ! D(yzvxtu)
        if (Mul(iSZ,iSY) == iSym .or. Mul(iSZ,iSV) == iSym .or. Mul(iSX,iST) == iSym) then
          call Add_MKBNEVD(iY,iZ,iV,iX,iT,iU,iG3)
        end if
        ! D(tuyzvx)
        if (Mul(iSU,iST) == iSym .or. Mul(iSU,iSY) == iSym .or. Mul(iSZ,iSV) == iSym) then
          call Add_MKBNEVD(iT,iU,iY,iZ,iV,iX,iG3)
        end if
      end if
      ! D(yztuvx)
      if (Mul(iSZ,iSY) == iSym .or. Mul(iSZ,iST) == iSym .or. Mul(iSU,iSV) == iSym) then
        call Add_MKBNEVD(iY,iZ,iT,iU,iV,iX,iG3)
      end if
      ! D(vxyztu)
      if (Mul(iSX,iSV) == iSym .or. Mul(iSX,iSY) == iSym .or. Mul(iSZ,iST) == iSym) then
        call Add_MKBNEVD(iV,iX,iY,iZ,iT,iU,iG3)
      end if
    end if

    if (iT == iU .and. iV == iX .and. iY == iZ) cycle
    if (iT == iU .and. iV == iZ .and. iX == iY) cycle
    if (iX == iV .and. iT == iZ .and. iU == iY) cycle
    if (iZ == iY .and. iV == iU .and. iX == iT) cycle

    ! D(utxvzy)
    if (Mul(iST,iSU) == iSym .or. Mul(iST,iSX) == iSym .or. Mul(iSV,iSZ) == iSym) then
      call Add_MKBNEVD(iU,iT,iX,iV,iZ,iY,iG3)
    end if

    if (iTU == iVX .and. iVX == iYZ) cycle
    if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
      ! D(xvutzy)
      if (Mul(iSV,iSX) == iSym .or. Mul(iSV,iSU) == iSym .or. Mul(iST,iSZ) == iSym) then
        call Add_MKBNEVD(iX,iV,iU,iT,iZ,iY,iG3)
      end if
      ! D(zyxvut)
      if (Mul(iSY,iSZ) == iSym .or. Mul(iSY,iSX) == iSym .or. Mul(iSV,iSU) == iSym) then
        call Add_MKBNEVD(iZ,iY,iX,iV,iU,iT,iG3)
      end if
      ! D(utzyxv)
      if (Mul(iST,iSU) == iSym .or. Mul(iST,iSZ) == iSym .or. Mul(iSY,iSX) == iSym) then
        call Add_MKBNEVD(iU,iT,iZ,iY,iX,iV,iG3)
      end if
    end if
    ! D(zyutxv)
    if (Mul(iSY,iSZ) == iSym .or. Mul(iSY,iSU) == iSym .or. Mul(iST,iSX) == iSym) then
      call Add_MKBNEVD(iZ,iY,iU,iT,iX,iV,iG3)
    end if
    ! D(xvzyut)
    if (Mul(iSV,iSX) == iSym .or. Mul(iSV,iSZ) == iSym .or. Mul(iSY,iSU) == iSym) then
      call Add_MKBNEVD(iX,iV,iZ,iY,iU,iT,iG3)
    end if
  end do

  call mma_deallocate(idxG3)

  call mma_deallocate(BDER1)
  call mma_deallocate(BDER2)

  return

contains

  subroutine Add_MKBNEVD(ITABS_,IUABS_,IVABS_,IXABS_,IYABS_,IZABS_,iG3_)

  use caspt2_module, only: NTUES, IASYM
  use SUPERINDEX, only: KTU
  use Symmetry_Info, only: Mul

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IUABS_, IVABS_, IXABS_, IYABS_, IZABS_, iG3_
  integer(kind=iwp) :: ITU_, IWABS, IXY_

  G3VAL = ex3(itabs_,iuabs_,ivabs_,ixabs_,iyabs_,izabs_,G1,G2,G3(iG3_))
  DG3VAL = Zero

  !! A
  iTU_ = KTU(IUABS_,ITABS_)-NTUES(ISYM)
  if (Mul(IASYM(IUABS_),IASYM(ITABS_)) == iSym) then
    do IWABS = 1, nAshT
      !! A(a'b',ab) = (ce|da) E(b'a',ce,db) --> A(ut,wz) = (vx|yw) E(tu,vx,yz)
      if (Mul(IASYM(IWABS),IASYM(IZABS_)) == iSym) then
        iXY_ = KTU(IWABS,IZABS_)-NTUES(ISYM)
        Gder(ivabs_,ixabs_,iyabs_,iwabs) = Gder(ivabs_,ixabs_,iyabs_,iwabs) + Half*BDER1(itu_,ixy_)*G3VAL
        DG3VAL = DG3VAL + Half*BDER1(itu_,ixy_)*Gact(ivabs_,ixabs_,iyabs_,iwabs)
      end if
      !! A(a'b',ab) = (ce|da) E(b'a',db,ce) --> A(ut,wx) = (yz|vw) E(tu,vx,yz)
      if (Mul(IASYM(IWABS),IASYM(IXABS_)) == iSym) then
        iXY_ = KTU(IWABS,IXABS_)-NTUES(ISYM)
        Gder(iyabs_,izabs_,ivabs_,iwabs) = Gder(iyabs_,izabs_,ivabs_,iwabs) + Half*BDER1(itu_,ixy_)*G3VAL
        DG3VAL = DG3VAL + Half*BDER1(itu_,ixy_)*Gact(iyabs_,izabs_,ivabs_,iwabs)
      end if

      !! A(a'b',ab) = (be|cd) E(b'a',ae,cd) --> A(ut,vw) = (wx|yz) E(tu,vx,yz)
      if (Mul(IASYM(IVABS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IVABS_,IWABS)-NTUES(ISYM)
        Gder(iwabs,ixabs_,iyabs_,izabs_) = Gder(iwabs,ixabs_,iyabs_,izabs_) - Half*BDER1(itu_,ixy_)*G3VAL
        DG3VAL = DG3VAL - Half*BDER1(itu_,ixy_)*Gact(iwabs,ixabs_,iyabs_,izabs_)
      end if
      !! A(a'b',ab) = (be|cd) E(b'a',cd,ae) --> A(ut,yw) = (wz|vx) E(tu,vx,yz)
      if (Mul(IASYM(IYABS_),IASYM(IWABS)) == iSym) then
        iXY_ = KTU(IYABS_,IWABS)-NTUES(ISYM)
        Gder(iwabs,izabs_,ivabs_,ixabs_) = Gder(iwabs,izabs_,ivabs_,ixabs_) - Half*BDER1(itu_,ixy_)*G3VAL
        DG3VAL = DG3VAL - Half*BDER1(itu_,ixy_)*Gact(iwabs,izabs_,ivabs_,ixabs_)
      end if
    end do
  end if

  !! D
  do IWABS = 1, nAshT
    !! D(a'b',ab) = -(ce|bd) E(ce,aa',b'd) --> D(xy,vw) = -(tu|wz) E(tu,vx,yz)
    iTU_ = KTU(IXABS_,IYABS_)-NTUES(ISYM)
    if (Mul(IASYM(IXABS_),IASYM(IYABS_)) == iSym) then
    if (Mul(IASYM(IVABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IVABS_,IWABS)-NTUES(ISYM)
      Gder(itabs_,iuabs_,iwabs,izabs_) = Gder(itabs_,iuabs_,iwabs,izabs_) + Half*BDER2(itu_,ixy_)*G3VAL
      DG3VAL = DG3VAL + Half*BDER2(itu_,ixy_)*Gact(itabs_,iuabs_,iwabs,izabs_)
    end if
    end if
    !! D(a'b',ab) = -(ce|bd) E(aa',b'd,ce) --> D(uv,tw) = -(yz|wx) E(tu,vx,yz)
    iTU_ = KTU(IUABS_,IVABS_)-NTUES(ISYM)
    if (Mul(IASYM(IUABS_),IASYM(IVABS_)) == iSym) then
    if (Mul(IASYM(ITABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(ITABS_,IWABS)-NTUES(ISYM)
      Gder(iyabs_,izabs_,iwabs,ixabs_) = Gder(iyabs_,izabs_,iwabs,ixabs_) + Half*BDER2(itu_,ixy_)*G3VAL
      DG3VAL = DG3VAL + Half*BDER2(itu_,ixy_)*Gact(iyabs_,izabs_,iwabs,ixabs_)
    end if
    end if

    !! D(a'b',ab) = +(ce|da) E(ce,da',b'b) --> D(xy,wz) = +(tu|vw) E(tu,vx,yz)
    iTU_ = KTU(IXABS_,IYABS_)-NTUES(ISYM)
    if (Mul(IASYM(IXABS_),IASYM(IYABS_)) == iSym) then
    if (Mul(IASYM(IWABS),IASYM(IZABS_)) == iSym) then
      iXY_ = KTU(IWABS,IZABS_)-NTUES(ISYM)
      Gder(itabs_,iuabs_,ivabs_,iwabs) = Gder(itabs_,iuabs_,ivabs_,iwabs) - Half*BDER2(itu_,ixy_)*G3VAL
      DG3VAL = DG3VAL - Half*BDER2(itu_,ixy_)*Gact(itabs_,iuabs_,ivabs_,iwabs)
    end if
    end if
    !! D(a'b',ab) = +(ce|da) E(da',b'b,ce) --> D(uv,wx) = +(yz|tw) E(tu,vx,yz)
    iTU_ = KTU(IUABS_,IVABS_)-NTUES(ISYM)
    if (Mul(IASYM(IUABS_),IASYM(IVABS_)) == iSym) then
    if (Mul(IASYM(IWABS),IASYM(IXABS_)) == iSym) then
      iXY_ = KTU(IWABS,IXABS_)-NTUES(ISYM)
      Gder(iyabs_,izabs_,itabs_,iwabs) = Gder(iyabs_,izabs_,itabs_,iwabs) - Half*BDER2(itu_,ixy_)*G3VAL
      DG3VAL = DG3VAL - Half*BDER2(itu_,ixy_)*Gact(iyabs_,izabs_,itabs_,iwabs)
    end if
    end if
  end do

  call derex3(itabs_,iuabs_,ivabs_,ixabs_,iyabs_,izabs_,DG1,DG2,DG3(iG3_),DG3VAL)

  end subroutine Add_MKBNEVD

end subroutine BDND

subroutine BDNE(iSym,NAS,BDER,SDER,G1,G2,DG1,DG2)

  use caspt2_module, only: NAES
  use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk

  implicit none

  integer(kind=iwp), intent(in) :: iSym, NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS), SDER(NAS,NAS), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)
  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT)

  integer(kind=iwp) :: it, itabs, iu, iuabs, ivabs, ixabs, iyabs
  integer(kind=iwp) :: ID, nTask

  nTask = NAS
  call Init_Tsk(ID,nTask)

  do while (Rsv_Tsk(ID,it))
    itabs = it + NAES(iSym)
    do iu = 1, NAS
      iuabs = iu + NAES(iSym)
      ! K(t,u) = h(v,u)*(2del(v,t)-R(t,v))
      do ivabs = 1, NASHT
        Hder(ivabs,iuabs) = Hder(ivabs,iuabs) - BDER(it,iu)*G1(itabs,ivabs)
        if (itabs == ivabs) Hder(ivabs,iuabs) = Hder(ivabs,iuabs) + Two*BDER(it,iu)
        DG1(itabs,ivabs) = DG1(itabs,ivabs) - BDER(it,iu)*Hact(ivabs,iuabs)
      end do
      ! K(t,u) = (xv|yu)*(2*del(t,y)*D(x,v)-D(y,t,x,v)-del(x,t)*D(y,v))
      do ivabs = 1, NASHT
        do ixabs = 1, NASHT
          do iyabs = 1, NASHT
            Gder(ixabs,ivabs,iyabs,iuabs) = Gder(ixabs,ivabs,iyabs,iuabs) - BDER(it,iu)*G2(iyabs,itabs,ixabs,ivabs)
            DG2(iyabs,itabs,ixabs,ivabs) = DG2(iyabs,itabs,ixabs,ivabs) - BDER(it,iu)*Gact(ixabs,ivabs,iyabs,iuabs)
            if (itabs == iyabs) then
              Gder(ixabs,ivabs,iyabs,iuabs) = Gder(ixabs,ivabs,iyabs,iuabs) + Two*BDER(it,iu)*G1(ixabs,ivabs)
              DG1(ixabs,ivabs) = DG1(ixabs,ivabs) + Two*BDER(it,iu)*Gact(ixabs,ivabs,iyabs,iuabs)
            end if
            if (ixabs == itabs) then
              Gder(ixabs,ivabs,iyabs,iuabs) = Gder(ixabs,ivabs,iyabs,iuabs) - BDER(it,iu)*G1(iyabs,ivabs)
              DG1(iyabs,ivabs) = DG1(iyabs,ivabs) - BDER(it,iu)*Gact(ixabs,ivabs,iyabs,iuabs)
            end if
          end do
        end do
      end do
      ! S(t,u) = 2del(t,u)-G1(t,u)
      DG1(itabs,iuabs) = DG1(itabs,iuabs) - SDER(it,iu)
    end do
  end do

  call Free_Tsk(ID)

  return

end subroutine BDNE

subroutine BDNF(iSym,iCase,NAS,NG3,BDER0,SDER0,G2,G3,DG2,DG3)

  use caspt2_module, only: NTU, NTUES, NTGEUES, NTGTUES, IASYM
  use SUPERINDEX, only: KTU, MTGEU, MTGTU, MTU
  use Symmetry_Info, only: Mul
  use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk

  implicit none

  integer(kind=iwp), intent(in) :: iSym, iCase, NAS, NG3
  real(kind=wp), intent(in) :: BDER0(NAS,NAS), SDER0(NAS,NAS), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  real(kind=wp), intent(inout) :: DG2(nAshT,nAshT,nAshT,nAshT), DG3(NG3)

  integer(kind=iwp) :: it, itabs, itu, ituabs, itu0, iu, iuabs, ivabs, iwabs, ixabs, ixy, ixyabs, ixy0, iyabs, iyx, izabs
  integer(kind=iwp) :: iLUID, iG3, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IV, IVX, IX, IXYZS, IY, IYZ, IZ
  integer(kind=iwp) :: iTgeUabs, iTgtUabs, iXgeYabs, iXgtYabs
  real(kind=wp) :: BDERval, DG3VAL, G3VAL, SDERval

  real(kind=wp), allocatable :: BDER(:,:), SDER(:,:)
  integer(kind=byte), allocatable :: idxG3(:,:)
  integer(kind=iwp) :: ID, nTask

  call mma_allocate(BDER,NTU(iSym),NTU(iSym),Label='BDER_F')
  call mma_allocate(SDER,NTU(iSym),NTU(iSym),Label='SDER_F')

  BDER(:,:) = Zero
  SDER(:,:) = Zero

  itabs = 0
  iuabs = 0
  ixabs = 0
  iyabs = 0
  do itu0 = 1, NAS
    if (iCase == 8) then
      itgeuabs = itu0 + NTGEUES(iSym)
      itabs    = MTGEU(1,itgeuabs)
      iuabs    = MTGEU(2,itgeuabs)
    else if (iCase == 9) then
      itgtuabs = itu0 + NTGTUES(iSym)
      itabs    = MTGTU(1,itgtuabs)
      iuabs    = MTGTU(2,itgtuabs)
    end if
    itu = KTU(itabs,iuabs) - NTUES(iSym)
    do ixy0 = 1, NAS
      if (iCase == 8) then
        ixgeyabs = ixy0 + NTGEUES(iSym)
        ixabs    = MTGEU(1,ixgeyabs)
        iyabs    = MTGEU(2,ixgeyabs)
      else if (iCase == 9) then
        ixgtyabs = ixy0 + NTGTUES(iSym)
        ixabs    = MTGTU(1,ixgtyabs)
        iyabs    = MTGTU(2,ixgtyabs)
      end if
      ixy = KTU(ixabs,iyabs) - NTUES(iSym)
      iyx = KTU(iyabs,ixabs) - NTUES(iSym)
      if (iCase == 8) then
        BDER(itu,ixy) = BDER(itu,ixy) + BDER0(itu0,ixy0)
        BDER(itu,iyx) = BDER(itu,iyx) + BDER0(itu0,ixy0)
        SDER(itu,ixy) = SDER(itu,ixy) + SDER0(itu0,ixy0)
        SDER(itu,iyx) = SDER(itu,iyx) + SDER0(itu0,ixy0)
      end if
      if (iCase == 9) then
        if (itabs == iuabs) cycle
        if (ixabs == iyabs) cycle
        BDER(itu,ixy) = BDER(itu,ixy) + BDER0(itu0,ixy0)
        BDER(itu,iyx) = BDER(itu,iyx) - BDER0(itu0,ixy0)
        SDER(itu,ixy) = SDER(itu,ixy) + SDER0(itu0,ixy0)
        SDER(itu,iyx) = SDER(itu,iyx) - SDER0(itu0,ixy0)
      end if
    end do
  end do

  BDER(:,:) = Two*BDER(:,:)
  SDER(:,:) = Two*SDER(:,:)

  nTask = NTU(iSym)
  call Init_Tsk(ID,nTask)

  do while (Rsv_Tsk(ID,itu))
    ituabs = itu + NTUES(iSym)
    itabs  = MTU(1,ituabs)
    iuabs  = MTU(2,ituabs)
    do ixy = 1, NTU(iSym)
      ixyabs = ixy + NTUES(iSym)
      ixabs  = MTU(1,ixyabs)
      iyabs  = MTU(2,ixyabs)

      BDERval = BDER(ITU,IXY)

      !! B derivative: Hbar contributions
      do ivabs = 1, NASHT
        Hbder(ixabs,ivabs) = Hbder(ixabs,ivabs) - BDERval*G2(itabs,ivabs,iuabs,iyabs)
        Hbder(iyabs,ivabs) = Hbder(iyabs,ivabs) - BDERval*G2(itabs,ixabs,iuabs,ivabs)
        DG2(itabs,ivabs,iuabs,iyabs) = DG2(itabs,ivabs,iuabs,iyabs) - BDERval*Hbar(ixabs,ivabs)
        DG2(itabs,ixabs,iuabs,ivabs) = DG2(itabs,ixabs,iuabs,ivabs) - BDERval*Hbar(iyabs,ivabs)
      end do

      !! B derivative: ERI contributions
      do ivabs = 1, NASHT
        do iwabs = 1, NASHT
          do izabs = 1, NASHT
            if (ivabs == izabs) then
              Gder(ivabs,iwabs,ixabs,izabs) = Gder(ivabs,iwabs,ixabs,izabs) - Half*BDERval*G2(itabs,iwabs,iuabs,iyabs)
              DG2(itabs,iwabs,iuabs,iyabs) = DG2(itabs,iwabs,iuabs,iyabs) - Half*BDERval*Gact(ivabs,iwabs,ixabs,izabs)
            end if
            if (iyabs == ivabs) then
              Gder(ivabs,iwabs,ixabs,izabs) = Gder(ivabs,iwabs,ixabs,izabs) - Half*BDERval*G2(itabs,izabs,iuabs,iwabs)
              DG2(itabs,izabs,iuabs,iwabs) = DG2(itabs,izabs,iuabs,iwabs) - Half*BDERval*Gact(ivabs,iwabs,ixabs,izabs)
            end if
            if (ixabs == ivabs) then
              Gder(ivabs,iwabs,iyabs,izabs) = Gder(ivabs,iwabs,iyabs,izabs) - Half*BDERval*G2(itabs,iwabs,iuabs,izabs)
              DG2(itabs,iwabs,iuabs,izabs) = DG2(itabs,iwabs,iuabs,izabs) - Half*BDERval*Gact(ivabs,iwabs,iyabs,izabs)
            end if
            if (ivabs == izabs) then
              Gder(ivabs,iwabs,iyabs,izabs) = Gder(ivabs,iwabs,iyabs,izabs) - Half*BDERval*G2(itabs,ixabs,iuabs,iwabs)
              DG2(itabs,ixabs,iuabs,iwabs) = DG2(itabs,ixabs,iuabs,iwabs) - Half*BDERval*Gact(ivabs,iwabs,iyabs,izabs)
            end if
          end do
        end do
      end do

      SDERval = SDER(ITU,IXY)

      !! S derivative
      DG2(itabs,ixabs,iuabs,iyabs) = DG2(itabs,ixabs,iuabs,iyabs) + SDERval
    end do
  end do

  call Free_Tsk(ID)

  call mma_deallocate(SDER)

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID=0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  !! The remaining D3 contributions
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
    ituvs=Mul(IST,Mul(ISU,ISV))
    ixyzs=Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU=iT+NASHT*(iU-1)
    iVX=iV+NASHT*(iX-1)
    iYZ=iY+NASHT*(iZ-1)

    G3VAL = G3(iG3)
    DG3VAL = Zero

    ! F(tuvxyz) -> BC(vut,xyz)
    if (Mul(iST,iSV) == iSym) call Add_MKBNEVF(iT,iU,iV,iX,iY,iZ,G3VAL,DG3VAL)

    if (iTU /= iVX .or. iVX /= iYZ) then
      if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
        ! F(vxtuyz) -> BC(txv,uyz)
        if (Mul(iSV,iST) == iSym) call Add_MKBNEVF(iV,iX,iT,iU,iY,iZ,G3VAL,DG3VAL)
        ! F(yzvxtu) -> BC(vzy,xtu)
        if (Mul(iSY,iSV) == iSym) call Add_MKBNEVF(iY,iZ,iV,iX,iT,iU,G3VAL,DG3VAL)
        ! F(tuyzvx) -> BC(yut,zvx)
        if (Mul(iST,iSY) == iSym) call Add_MKBNEVF(iT,iU,iY,iZ,iV,iX,G3VAL,DG3VAL)
      end if
      ! F(yztuvx) -> BC(tzy,uvx)
      if (Mul(iSY,iST) == iSym) call Add_MKBNEVF(iY,iZ,iT,iU,iV,iX,G3VAL,DG3VAL)
      ! F(vxyztu) -> BC(yxv,ztu)
      if (Mul(iSV,iSY) == iSym) call Add_MKBNEVF(iV,iX,iY,iZ,iT,iU,G3VAL,DG3VAL)
    end if

    if ((iT /= iU .or. iV /= iX .or. iY /= iZ) .and. (iT /= iU .or. iV /= iZ .or. iX /= iY) .and. &
        (iX /= iV .or. iT /= iZ .or. iU /= iY) .and. (iZ /= iY .or. iV /= iU .or. iX /= iT)) then
      ! F(utxvzy) -> BC(xtu,vzy)
      if (Mul(iSU,iSX) == iSym) call Add_MKBNEVF(iU,iT,iX,iV,iZ,iY,G3VAL,DG3VAL)

      if (iTU /= iVX .or. iVX /= iYZ) then
        if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
          ! F(xvutzy) -> BC(uvx,tzy)
          if (Mul(iSX,iSU) == iSym) call Add_MKBNEVF(iX,iV,iU,iT,iZ,iY,G3VAL,DG3VAL)
          ! F(zyxvut) -> BC(xyz,vut)
          if (Mul(iSZ,iSX) == iSym) call Add_MKBNEVF(iZ,iY,iX,iV,iU,iT,G3VAL,DG3VAL)
          ! F(utzyxv) -> BC(ztu,yxv)
          if (Mul(iSU,iSZ) == iSym) call Add_MKBNEVF(iU,iT,iZ,iY,iX,iV,G3VAL,DG3VAL)
        end if
        ! F(zyutxv) -> BC(uyz,txv)
        if (Mul(iSZ,iSU) == iSym) call Add_MKBNEVF(iZ,iY,iU,iT,iX,iV,G3VAL,DG3VAL)
        ! F(xvzyut) -> BC(zvx,yut)
        if (Mul(iSX,iSZ) == iSym) call Add_MKBNEVF(iX,iV,iZ,iY,iU,iT,G3VAL,DG3VAL)
      end if
    end if

    DG3(iG3) = DG3(iG3) + DG3VAL
  end do

  call mma_deallocate(idxG3)

  call mma_deallocate(BDER)

  return

contains

  subroutine Add_MKBNEVF(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_,g3val_,dg3val_)

  use caspt2_module, only: NTUES, IASYM
  use SUPERINDEX, only: KTU
  use Symmetry_Info, only: Mul

  implicit none

  integer(kind=iwp), intent(in) :: ITABS_, IVABS_, IYABS_, IUABS_, IXABS_, IZABS_
  real(kind=wp), intent(in) :: g3val_
  real(kind=wp), intent(inout) :: dg3val_
  integer(kind=iwp) :: ITU_, IWABS, IXY_

  iTU_ = KTU(ITABS_,IUABS_)-NTUES(ISYM)
  do IWABS = 1, nAshT
    if (Mul(IASYM(IWABS),IASYM(IYABS_)) == iSym) then
      iXY_ = KTU(IWABS,IYABS_)-NTUES(ISYM)
      dg3val_ = dg3val_ - BDER(iTU_,iXY_)*Gact(iwabs,ixabs_,ivabs_,izabs_)
      Gder(iwabs,ixabs_,ivabs_,izabs_) = Gder(iwabs,ixabs_,ivabs_,izabs_) - BDER(iTU_,iXY_)*g3val_
    end if
    if (Mul(IASYM(IXABS_),IASYM(IWABS)) == iSym) then
      iXY_ = KTU(IXABS_,IWABS)-NTUES(ISYM)
      dg3val_ = dg3val_ - BDER(iTU_,iXY_)*Gact(iwabs,iyabs_,ivabs_,izabs_)
      Gder(iwabs,iyabs_,ivabs_,izabs_) = Gder(iwabs,iyabs_,ivabs_,izabs_) - BDER(iTU_,iXY_)*g3val_
    end if
  end do

  end subroutine Add_MKBNEVF

end subroutine BDNF

subroutine BDNG(iSym,NAS,BDER,SDER,G1,G2,DG1,DG2)

  use caspt2_module, only: NAES
  use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk

  implicit none

  integer(kind=iwp), intent(in) :: iSym, NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS), SDER(NAS,NAS), G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)
  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT)

  integer(kind=iwp) :: it, itabs, iu, iuabs, ivabs, ixabs, iyabs
  integer(kind=iwp) :: ID, nTask

  nTask = NAS
  call Init_Tsk(ID,nTask)

  do while (Rsv_Tsk(ID,it))
    itabs = it + NAES(iSym)
    do iu = 1, NAS
      iuabs = iu + NAES(iSym)
      ! K(t,u) = -h(u,v)*G1(t,v)
      do ivabs = 1, NASHT
        Hder(iuabs,ivabs) = Hder(iuabs,ivabs) - BDER(it,iu)*G1(itabs,ivabs)
        DG1(itabs,ivabs) = DG1(itabs,ivabs) - BDER(it,iu)*Hact(iuabs,ivabs)
      end do
      ! K(t,u) = -(xy|uv)*G2(tv,xy)
      do ivabs = 1, NASHT
        do ixabs = 1, NASHT
          do iyabs = 1, NASHT
            Gder(ixabs,iyabs,iuabs,ivabs) = Gder(ixabs,iyabs,iuabs,ivabs) - BDER(it,iu)*G2(itabs,ivabs,ixabs,iyabs)
            DG2(itabs,ivabs,ixabs,iyabs) = DG2(itabs,ivabs,ixabs,iyabs) - BDER(it,iu)*Gact(ixabs,iyabs,iuabs,ivabs)
          end do
        end do
      end do
      ! S(t,u) = G1(t,u)
      DG1(itabs,iuabs) = DG1(itabs,iuabs) + SDER(it,iu)
    end do
  end do

  call Free_Tsk(ID)

  return

end subroutine BDNG

subroutine BDN_G3(DG1,DG2,DG3)

  use caspt2_module, only: IASYM, NG3
  use Symmetry_Info, only: Mul

  implicit none

  real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DG3(*)

  integer(kind=iwp) :: iLUID, iG3, IST, ISU, ISV, ISX, ISY, ISZ, IT, ITUVS, IU, IV, IX, IXYZS, IY, IZ
  real(kind=wp) :: G3VAL

  integer(kind=byte), allocatable :: idxG3(:,:)

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID=0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

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
    ituvs=Mul(IST,Mul(ISU,ISV))
    ixyzs=Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle

    G3VAL = DG3(iG3)

    !! remaining G3 transformation in mkfg3.f
    if (iY == iX) DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - G3VAL
    if (iV == iU) DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - G3VAL
    if (iY == iU) DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - G3VAL
    if (iY == iX .and. iV == iU) DG1(iT,iZ) = DG1(iT,iZ) - G3VAL
  end do

  call mma_deallocate(idxG3)

  return

end subroutine BDN_G3

end module BDerNEV
