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
! Copyright (C) 1991, Roland Lindh                                     *
!               1996, Per Ake Malmqvist                                *
!***********************************************************************

subroutine Drv1El()
!***********************************************************************
!                                                                      *
! Object: driver for computation of one-electron matrices.             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1991                                             *
!***********************************************************************

use AMFI_Info, only: No_AMFI
use Basis_Info, only: dbsc, nCnttp, PAMexp
use GeoList, only: Centr, Chrg
use MpmC, only: Coor_MpM
use PrpPnt, only: Den, Occ, Vec
use External_Centers, only: Dxyz, AMP_Center, EF_Centers, DMS_Centers, nDMS, nEF, nOrdEF, nOrd_XF, nWel, OAM_Center, OMQ_Center, &
                            Wel_Info, XF
use Symmetry_Info, only: iChBas
use Gateway_global, only: Primitive_Pass, PrPrt, Short, SW_FileOrb
use PAM2, only: iPAMcount, iPAMPrim, kCnttpPAM
use DKH_Info, only: BSS, DKroll
use Sizes_of_Seward, only: S
use Gateway_Info, only: Do_FckInt, DoFMM, EMFR, GIAO, kVector, lAMFI, lMXTC, lRel, NEMO, PotNuc, Vlct
use Integral_interfaces, only: int_kernel, int_mem
use Property_Label, only: PLabel
#ifdef _FDE_
use Embedding_Global, only: embInt, embPot, embPotInBasis, embPotPath
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
#include "wldata.fh"
#include "oneswi.fh"
#include "warnings.h"
integer(kind=iwp) :: i, i2, i3, iAddr, iAtom_Number, iB, iC, iChO, iChO1, iChO2, iChOx, iChOxx, iChOxy, iChOxz, iChOy, iChOyx, &
                     iChOyy, iChOyz, iChOz, iChOzx, iChOzy, iChOzz, iCmp, iCnt, iCnttp, iComp, iD, iDisk, iDMS, idum(1), iEF, &
                     iLow, iMltpl, iOpt, iPAMBas, iPAMf, iPAMltpl, iPrint, iRC, iRout, iSym, iSymBx, iSymBy, iSymBz, iSymC, &
                     iSymCX, iSymCXY, iSymCy, iSymCz, iSymD, iSymLx, iSymLy, iSymLz, iSymR(0:3), iSymRx, iSymRy, iSymRz, iSymX, &
                     iSymxLx, iSymxLy, iSymxLz, iSymXY, iSymXZ, iSymY, iSymyLx, iSymyLy, iSymyLz, iSymYZ, iSymZ, iSymzLx, iSymzLy, &
                     iSymzLz, iSyXYZ, iTemp, iTol, iWel, ix, ixyz, iy, iz, jx, jxyz, jy, jz, kCnttpPAM_, lOper, LuTmp, mCnt, &
                     mComp, mDMS, mMltpl, mOrdOp, nB, nComp, nOrdOp, nPAMltpl
real(kind=wp) :: Ccoor(3), dum(1), Fact, rHrmt
logical(kind=iwp) :: lECPnp, lECP, lPAM2np, lPAM2, lPP, lFAIEMP
character(len=8) :: Label
character(len=512) :: FName
integer(kind=iwp), external :: IrrFnc, MltLbl, n2Tri
! ipList: list of pointers to the integrals of each component
!         of the operator
! OperI: list which irreps a particular component of the operator
!        belongs to
! OperC: list the character of each component of the operator
! CoorO: list of origins of the operator, one for each component
integer(kind=iwp), allocatable :: ipList(:), OperI(:), OperC(:), iAtmNr2(:)
real(kind=wp), allocatable :: CoorO(:), Nuc(:), KnE_Int(:), NA_Int(:), FragP(:), OneHam(:), PtEl(:), PtNuc(:), SumEl(:), &
                              SumNuc(:), Charge2(:)
procedure(int_kernel) :: AMPInt, CntInt, D1Int, DMSInt, dTdmu_Int, EFInt, EMFInt, FragPint, KneInt, KneInt_GIAO, M1Int, M2Int, &
                         MltInt, MltInt_GIAO, MVeInt, NAInt, NAInt_GIAO, OAMInt, OMQInt, P_Int, PAM2Int, PPInt, PrjInt, PXInt, &
                         PXPInt, QpVInt, SROInt, VeInt, VPInt, WelInt, XFdInt
procedure(int_mem) :: AMPMem, CntMem, D1Mem, DMSMem, dTdmu_Mem, EFMem, EMFMem, FragPMem, KneMem, KneMem_GIAO, M1Mem, M2Mem, &
                      MltMem, MltMem_GIAO, MVeMem, NAMem, NAMem_GIAO, OAMMem, OMQMem, P_Mem, PAM2Mem, PPMem, PrjMem, PXMem, &
                      PXPMem, QpVMem, SROMem, VeMem, VPMem, WelMem, XFdMem
#ifdef _FDE_
! Embedding
integer(kind=iwp) :: iEMb, iunit
real(kind=wp), allocatable :: Emb_Int(:)
integer(kind=iwp), external :: isFreeUnit
procedure(int_kernel) :: embPotKernel
procedure(int_mem) :: embPotMem
#endif
#ifdef _GEN1INT_
integer(kind=iwp) :: nAtoms, jCnt
! These won't actually be called, but need to be passed around
procedure(int_kernel) :: DumInt
procedure(int_mem) :: DumMem
#endif

iRout = 131
iPrint = nPrint(iRout)

call StatusLine(' Seward:',' Computing 1-electron integrals')

call Set_Basis_Mode('Valence')
call Setup_iSD()
#ifdef _FDE_
if (embPot) call EmbPotRdRun()
#endif

lPAM2 = .false.
lECP = .false.
lPP = .false.
lFAIEMP = .false.
do i=1,nCnttp
  lPam2 = lPam2 .or. dbsc(i)%lPam2
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
  lFAIEMP = lFAIEMP .or. dbsc(i)%Frag
end do

! set center selector in OneSwi to all centers (default)

NDDO = .false.
if (Prprt .and. DKroll) then
  call WarningMessage(2,'Prprt and DKroll options cannot be combined!')
  call Quit_OnUserError()
end if

! We will always compute the following one-electron integrals per default.
! 1) Multipole moments up to quadrupole moments
! 2) Kinetic energy
! 3) Nuclear Attraction
! 4) ECP contributions
! 5) One-Electron Hamiltonian
! 6) Mass-Velocity
! 7) Darwin 1-electron contact term

lECPnp = lECP
lPAM2np = lPAM2
if (DKroll .and. Primitive_Pass) then
  lECPnp = .false.
end if
if (Prprt) then
  FName = SW_FileOrb
  call GetDens(trim(FName),short,iPrint)
  call CollapseOutput(1,'   Molecular properties:')
  write(u6,'(3X,A)') '   ---------------------'
  write(u6,*)
end if
!***********************************************************************
!***********************************************************************
!1)                                                                    *
!                                                                      *
!     Multipole Moments starting with the overlap. If SEWARD is run in *
!     the property mode we will skip the overlap integrals.            *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = One
iLow = 0
mMltpl = S%nMltpl
if (Prprt) iLow = 1

! If not Douglas-Kroll and primitive pass do no property integrals

if (Primitive_Pass) then
  iLow = 0
  if (.not. DKroll) mMltpl = -1
end if

do iMltpl=iLow,S%nMltpl
  write(Label,'(A,I2)') 'Mltpl ',iMltpl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  call DCopy_(3,Coor_MPM(1,iMltpl+1),1,Ccoor,1)
  call Allocate_Auxiliary()
  iComp = 0
  do ix=iMltpl,0,-1
    if (mod(ix,2) == 0) then
      iSymX = 1
    else
      ixyz = 1
      iSymX = 2**IrrFnc(ixyz)
      if (Ccoor(1) /= Zero) iSymX = ibset(iSymX,0)
    end if
    do iy=iMltpl-ix,0,-1
      if (mod(iy,2) == 0) then
        iSymY = 1
      else
        ixyz = 2
        iSymY = 2**IrrFnc(ixyz)
        if (Ccoor(2) /= Zero) iSymY = ibset(iSymY,0)
      end if
      iz = iMltpl-ix-iy
      if (mod(iz,2) == 0) then
        iSymZ = 1
      else
        ixyz = 4
        iSymZ = 2**IrrFnc(ixyz)
        if (Ccoor(3) /= Zero) iSymZ = ibset(iSymZ,0)
      end if
      iChO = mod(ix,2)*iChBas(2)+mod(iy,2)*iChBas(3)+mod(iz,2)*iChBas(4)

      OperI(1+iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ))
      OperC(1+iComp) = iChO
      call DCopy_(3,Coor_MPM(1,iMltpl+1),1,CoorO(1+iComp*3),1)
      iComp = iComp+1
    end do
  end do

  call MltNuc(CoorO,Chrg,Centr,S%kCentr,Nuc,iMltpl,nComp)
  !--- pow hack
  if (iMltpl == 0) then
    call Put_dScalar('Total Nuclear Charge',Nuc(1))
  end if
  !--- pow hack
  nOrdOp = iMltpl
  call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Write FMM multipole moments to disk

  if ((.not. Prprt) .and. DoFMM) then
    write(Label,'(A,I2)') 'FMMInt',iMltpl
    call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

    ! FMM overlap distribution centres:
    ! Pretend they are 1-e integrals with three (x,y,z)
    ! components and write to disk in canonical order

    if (iMltpl == 0) then
      write(Label,'(A)') 'FMMCnX'
      call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,CoorO,nOrdOp+1,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      write(Label,'(A)') 'FMMCnY'
      call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,CoorO,nOrdOp+1,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      write(Label,'(A)') 'FMMCnZ'
      call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,CoorO,nOrdOp+1,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! For picture-change corrected integrals.

  if (DKroll .and. Primitive_Pass) then
    write(Label,'(A,I2)') 'pMp   ',iMltpl
    PLabel = 'MltInt'
    call FZero(Nuc,nComp)
    call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp+2,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iMltpl == 0) then
    ! these are overlap integrals...
    ! set center selector in OneSwi to single center...
    NDDO = .true.
    write(Label,'(A,I2)') 'MltplS',iMltpl
    call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    ! reset center selector in OneSwi to all centers...
    NDDO = .false.
  end if

  call Deallocate_Auxiliary()

end do
!***********************************************************************
!***********************************************************************
!1a)                                                                   *
!                                                                      *
!     PAM integrals                                                    *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (lPAM2np .and. (.not. Primitive_Pass)) then
  call molcas_open(28,'R_vect')
  PLabel = ' '
  rHrmt = One
  iPAMcount = 1
  do kCnttpPAM_=1,nCnttp

    kCnttpPAM = kCnttpPAM_
    nPAMltpl = dbsc(kCnttpPAM)%nPAM2

    if (nPAMltpl < 0) cycle

    iAddr = 1
    do iPAMltpl=0,nPAMltpl
      nOrdOp = iPAMltpl
      nComp = (iPAMltpl+1)*(iPAMltpl+2)/2
      iPAMPrim = int(dbsc(kCnttpPAM)%PAM2(iAddr))
      iPAMBas = int(dbsc(kCnttpPAM)%PAM2(iAddr+1))

      if ((iPAMBas /= 0) .and. (iPAMPrim /= 0)) then
        Ccoor(:) = Zero
        call Allocate_Auxiliary()
        do iComp=0,nComp-1
          call dcopy_(3,dbsc(kCnttpPAM)%Coor,1,CoorO(1+3*iComp),1)
        end do

        !**** Define symmetry properties of the operator:

        iComp = 0
        do ix=iPAMltpl,0,-1
          if (mod(ix,2) == 0) then
            iSymX = 1
          else
            ixyz = 1
            iSymX = 2**IrrFnc(ixyz)
            if (Ccoor(1) /= Zero) iSymX = ibset(iSymX,0)
          end if
          do iy=iPAMltpl-ix,0,-1
            if (mod(iy,2) == 0) then
              iSymY = 1
            else
              ixyz = 2
              iSymY = 2**IrrFnc(ixyz)
              if (Ccoor(2) /= Zero) iSymY = ibset(iSymY,0)
            end if
            iz = iPAMltpl-ix-iy
            if (mod(iz,2) == 0) then
              iSymZ = 1
            else
              ixyz = 4
              iSymZ = 2**IrrFnc(ixyz)
              if (Ccoor(3) /= Zero) iSymZ = ibset(iSymZ,0)
            end if
            iChO = mod(ix,2)*iChBas(2)+mod(iy,2)*iChBas(3)+mod(iz,2)*iChBas(4)

            OperC(1+iComp) = iChO
            OperI(1+iComp) = 1

            iComp = iComp+1
          end do
        end do

        !**** Loop over basis functions

        call mma_allocate(PAMexp,iPAMPrim,2,label='PAMexp')
        call dcopy_(iPAMPrim,dbsc(kCnttpPAM)%PAM2(iAddr+2),1,PAMexp(1,1),1)
        do iPAMf=1,iPAMBas
          call dcopy_(iPAMPrim,dbsc(kCnttpPAM)%PAM2(iAddr+2+iPAMPrim*iPAMf),1,PAMexp(1,2),1)
          write(Label,'(A,I2.2,I1.1,I2.2)') 'PAM',kCnttpPAM,iPAMltpl,iPAMf

          call dcopy_(nComp,[Zero],0,Nuc,1)

          call OneEl(PAM2Int,PAM2Mem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

          !iPAMcount = iPAMcount+1

        end do
        call mma_deallocate(PAMexp)
        call Deallocate_Auxiliary()
      end if
      iAddr = iAddr+2+iPAMPrim*(iPAMBas+1)
    end do
  end do
  close(28)
end if
!***********************************************************************
!***********************************************************************
!2)                                                                    *
!                                                                      *
!     Kinetic energy, nuclear attraction and ECP/PP integrals          *
!                                                                      *
!     Mass-velocity and One-electron Darwin contact term integrals.    *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = One
nComp = 1

if (.not. Prprt) then
  call Allocate_Auxiliary()
  call dcopy_(3,[Zero],0,CoorO,1)
  OperI(1) = 1
  OperC(1) = iChBas(1)

  Label = 'Kinetic '
  nOrdOp = 2
  call OneEl(KnEInt,KnEMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  nOrdOp = 0

  Label = 'Attract '
  call OneEl(NAInt,NAMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[PotNuc],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

# ifdef _FDE_
  ! Embedding
  if (embPot) then
    Label = 'Embpot '
    call OneEl(EmbPotKernel,EmbPotMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  end if
# endif

  ! set center selector in OneSwi to two center NA Int...
  NDDO = .true.
  Label = 'AttractS'
  call OneEl(NAInt,NAMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[PotNuc],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  ! reset center selector in OneSwi to all centers...
  NDDO = .false.
  if (.not. Primitive_Pass) then
    if (lECPnp) then
      Label = 'PrjInt  '
      call OneEl(PrjInt,PrjMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      Label = 'M1Int   '
      call OneEl(M1Int,M1Mem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      Label = 'M2Int   '
      call OneEl(M2Int,M2Mem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      Label = 'SROInt  '
      call OneEl(SROInt,SROMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    end if
    if (lPP) then
      Label = 'PPInt   '
      call OneEl(PPInt,PPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    end if
    if (allocated(XF)) then
      mOrdOp = nOrd_XF
      Label = 'XFdInt  '
      call OneEl(XFdInt,XFdMem,Label,ipList,OperI,nComp,CoorO,mOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    end if
    if (lRel) then
      Label = 'MassVel '
      nOrdOp = 4
      call OneEl(MVeInt,MVeMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      Label = 'Darwin  '
      nOrdOp = 0
      call OneEl(D1Int,D1Mem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    end if
  end if

  call Deallocate_Auxiliary()
end if
!***********************************************************************
!***********************************************************************
!8a)                                                                   *
!                                                                      *
!     Velocity integrals.                                              *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = -One
if (Vlct .and. (.not. Primitive_Pass)) then
  nOrdOp = 1
  Label = 'Velocity'
  nComp = 3
  call Allocate_Auxiliary()
  call dcopy_(3*nComp,[Zero],0,CoorO,1)
  ixyz = 1
  OperI(1) = 2**IrrFnc(ixyz)
  OperC(1) = iChBas(2)
  ixyz = 2
  OperI(1+1) = 2**IrrFnc(ixyz)
  OperC(1+1) = iChBas(3)
  ixyz = 4
  OperI(1+2) = 2**IrrFnc(ixyz)
  OperC(1+2) = iChBas(4)

  call dcopy_(3,[Zero],0,Nuc,1)
  call OneEl(VeInt,VeMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if    ! Vlct
!***********************************************************************
!***********************************************************************
!8b)                                                                   *
!                                                                      *
!     Electromagnetic field radiation integrals.                       *
!                                                                      *
!     Note that the integral is not symmetric or antisymmetric!        *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = -One ! Note used
if (EMFR .and. (.not. Primitive_Pass)) then

  ! The second term in eq. 14 of Bernadotte et al

  nOrdOp = 0
  Label = 'EMFR0'
  nComp = 2
  call Allocate_Auxiliary()
  ! Here we put in the k-vector
  call FZero(CoorO,3*nComp)
  call dcopy_(3,KVector,1,CoorO,1)

  ! The electromagnetic field operator contributes to all
  ! irreducible irreps, hence OperI=255. Since the operator
  ! itself is not symmetry adapted OperC is set to a dummy value.

  OperI(1) = 255
  OperI(1+1) = 255
  OperC(1) = 0 ! Dummy
  OperC(1+1) = 0 ! Dummy

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()

  ! The first (dominating) term in eq. 14 of Bernadotte et al

  nOrdOp = 1
  Label = 'EMFR'
  nComp = 12
  call Allocate_Auxiliary()
  ! Here we put in the k-vector
  call FZero(CoorO,3*nComp)
  call dcopy_(3,KVector,1,CoorO,1)

  ! The electromagnetic field operator contributes to all
  ! irreducible irreps, hence OperI=255. Since the operator
  ! itself is not symmetry adapted OperC is set to a dummy value.

  OperI(1) = 255
  OperI(1+1) = 255
  OperI(1+2) = 255
  OperI(1+3) = 255
  OperI(1+4) = 255
  OperI(1+5) = 255
  OperI(1+6) = 255
  OperI(1+7) = 255
  OperI(1+8) = 255
  OperI(1+9) = 255
  OperI(1+10) = 255
  OperI(1+11) = 255
  OperC(1) = 0 ! Dummy
  OperC(1+1) = 0 ! Dummy
  OperC(1+2) = 0 ! Dummy
  OperC(1+3) = 0 ! Dummy
  OperC(1+4) = 0 ! Dummy
  OperC(1+5) = 0 ! Dummy
  OperC(1+6) = 0 ! Dummy
  OperC(1+7) = 0 ! Dummy
  OperC(1+8) = 0 ! Dummy
  OperC(1+9) = 0 ! Dummy
  OperC(1+10) = 0 ! Dummy
  OperC(1+11) = 0 ! Dummy

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if    ! EMFR
!***********************************************************************
!***********************************************************************
!9)                                                                    *
!10)                                                                   *
!                                                                      *
!     Electric field integrals.                                        *
!     Electric field gradient integrals.                               *
!                                                                      *
!***********************************************************************
!***********************************************************************
ixyz = 1
iSymX = 2**IrrFnc(ixyz)
ixyz = 2
iSymY = 2**IrrFnc(ixyz)
ixyz = 4
iSymZ = 2**IrrFnc(ixyz)
ixyz = 3
iSymXY = 2**IrrFnc(ixyz)
ixyz = 5
iSymXZ = 2**IrrFnc(ixyz)
ixyz = 6
iSymYZ = 2**IrrFnc(ixyz)
ixyz = 7
iSyXYZ = 2**IrrFnc(ixyz)

PLabel = ' '
rHrmt = One
do nOrdOp=0,nOrdEF

  nComp = (nOrdOp+1)*(nOrdOp+2)/2

  call Allocate_Auxiliary()

  do iEF=1,nEF

    ! Note that this parsing is a bit different here!

    write(Label,'(A,I1,I5)') 'EF',nOrdOp,iEF
    Ccoor(:) = EF_Centers(:,iEF)

    iSymR(0) = 1
    if (Ccoor(1) /= Zero) iSymR(0) = ior(iSymR(0),iSymX)
    if (Ccoor(2) /= Zero) iSymR(0) = ior(iSymR(0),iSymY)
    if (Ccoor(3) /= Zero) iSymR(0) = ior(iSymR(0),iSymZ)
    if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero)) iSymR(0) = ior(iSymR(0),iSymXY)
    if ((Ccoor(1) /= Zero) .and. (Ccoor(3) /= Zero)) iSymR(0) = ior(iSymR(0),iSymXZ)
    if ((Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymR(0) = ior(iSymR(0),iSymYZ)
    if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymR(0) = ior(iSymR(0),iSyXYZ)

    ixyz = 1
    iSym = 2**IrrFnc(ixyz)
    if (Ccoor(1) /= Zero) iSym = ibset(iSym,0)
    iSymR(1) = iSym

    ixyz = 2
    iSym = 2**IrrFnc(ixyz)
    if (Ccoor(2) /= Zero) iSym = ibset(iSym,0)
    iSymR(2) = iSym

    ixyz = 4
    iSym = 2**IrrFnc(ixyz)
    if (Ccoor(3) /= Zero) iSym = ibset(iSym,0)
    iSymR(3) = iSym

    iComp = 0
    do ix=nOrdOp,0,-1
      do iy=nOrdOp-ix,0,-1
        iz = nOrdOp-ix-iy
        iComp = iComp+1

        iSymX = 1
        if (mod(ix,2) /= 0) iSymX = iSymR(1)
        iSymCX = MltLbl(iSymR(0),iSymX)
        iSymY = 1
        if (mod(iy,2) /= 0) iSymY = iSymR(2)
        iSymCXY = MltLbl(iSymCX,iSymY)
        iSymZ = 1
        if (mod(iz,2) /= 0) iSymZ = iSymR(3)

        OperI(1+(iComp-1)) = MltLbl(iSymCXY,iSymZ)
        OperC(1+(iComp-1)) = mod(ix,2)*iChBas(2)+mod(iy,2)*iChBas(3)+mod(iz,2)*iChBas(4)

        call dcopy_(3,Ccoor,1,CoorO(1+(iComp-1)*3),1)
      end do
    end do

    call EFNuc(CoorO,Chrg,Centr,S%kCentr,Nuc,nOrdOp)
    call OneEl(EFInt,EFMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! For picture-change corrected property integrals.

    if (DKroll .and. Primitive_Pass) then
      write(Label,'(A,I1,I5)') 'PP',nOrdOp,iEF
      PLabel = 'EFInt '
      call FZero(Nuc,nComp)
      call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp+2,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do
  call Deallocate_Auxiliary()

  ! For a properties calculation, read the EF values saved in a temp file
  ! and write the sum through Add_Info

  if (PrPrt .and. (nEF > 0)) then
    call mma_allocate(PtEl,nComp,label='PtEl')
    call mma_allocate(PtNuc,nComp,label='PtNuc')
    call mma_allocate(SumEl,nComp,label='SumEl')
    call mma_allocate(SumNuc,nComp,label='SumNuc')
    call FZero(SumEl,nComp)
    call FZero(SumNuc,nComp)
    ! Read and sum the values
    LuTmp = 10
    call DaName(LuTmp,'TMPPRP')
    iDisk = 0
    do iEf=1,nEF
      call dDaFile(LuTmp,2,PtEl,nComp,iDisk)
      call dDaFile(LuTmp,2,PtNuc,nComp,iDisk)
      call DaXpY_(nComp,One,PtEl,1,SumEl,1)
      call DaXpY_(nComp,One,PtNuc,1,SumNuc,1)
    end do
    call DaClos(LuTmp)
    ! set the tolerance according to the total number of centers
    ! (assuming error scales with sqrt(nEF))
    iTol = 5
    iTol = iTol-nint(Half*log10(real(nEF,kind=wp)))
    write(label,'(a,i1,a)') 'EF',nOrdOp,'   el'
    call Add_Info(label,SumEl,nComp,iTol)
    write(label,'(a,i1,a)') 'EF',nOrdOp,'  nuc'
    call Add_Info(label,SumNuc,nComp,iTol)
    call mma_deallocate(PtEl)
    call mma_deallocate(PtNuc)
    call mma_deallocate(SumEl)
    call mma_deallocate(SumNuc)
  end if

end do
!***********************************************************************
!***********************************************************************
!12)                                                                   *
!                                                                      *
!     Orbital angular momentum integrals.                              *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = -One
if (allocated(OAM_Center) .and. (.not. Primitive_Pass)) then
  Label = 'AngMom  '
  nComp = 3
  nOrdOp = 2
  call Allocate_Auxiliary()
  call dcopy_(3,OAM_Center,1,CoorO(1),1)
  call dcopy_(3,OAM_Center,1,CoorO(1+3),1)
  call dcopy_(3,OAM_Center,1,CoorO(1+6),1)
  call dcopy_(3,OAM_Center,1,Ccoor,1)
  ixyz = 1
  iSymX = 2**IrrFnc(ixyz)
  ixyz = 2
  iSymY = 2**IrrFnc(ixyz)
  ixyz = 4
  iSymZ = 2**IrrFnc(ixyz)
  iSymCx = iSymX
  if (Ccoor(1) /= Zero) iSymCx = ibset(iSymCx,0)
  iSymCy = iSymY
  if (Ccoor(2) /= Zero) iSymCy = ibset(iSymCy,0)
  iSymCz = iSymZ
  if (Ccoor(3) /= Zero) iSymCz = ibset(iSymCz,0)

  iSymLx = ior(MltLbl(iSymCy,iSymZ),MltLbl(iSymCz,iSymY))
  iChOx = iChBas(3)+iChBas(4)
  OperI(1) = iSymLx
  OperC(1) = iChOx
  iSymLy = ior(MltLbl(iSymCz,iSymX),MltLbl(iSymCx,iSymZ))
  iChOy = iChBas(4)+iChBas(2)
  OperI(1+1) = iSymLy
  OperC(1+1) = iChOy
  iSymLz = ior(MltLbl(iSymCx,iSymY),MltLbl(iSymCy,iSymX))
  iChOz = iChBas(2)+iChBas(3)
  OperI(1+2) = iSymLz
  OperC(1+2) = iChOz

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(OAMInt,OAMMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if   ! OAM_Center
!***********************************************************************
!***********************************************************************
!12b)                                                                  *
!                                                                      *
!     Velocity quadrupole.                                             *
!     (the companion symmetric combination to the angular momentum)    *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = -One
if (Vlct .and. (S%nMltpl >= 2) .and. (.not. Primitive_Pass)) then
  Label = 'MLTPV  2'
  nComp = 6
  nOrdOp = 2
  call Allocate_Auxiliary()

  ! Use origin for quadrupole moment
  call DCopy_(nComp,Coor_MPM(1,3),0,CoorO(1),3)
  call DCopy_(nComp,Coor_MPM(2,3),0,CoorO(1+1),3)
  call DCopy_(nComp,Coor_MPM(3,3),0,CoorO(1+2),3)
  call dCopy_(3,Coor_MPM(1,3),1,Ccoor,1)

  ixyz = 1
  iSymX = 2**IrrFnc(ixyz)
  ixyz = 2
  iSymY = 2**IrrFnc(ixyz)
  ixyz = 4
  iSymZ = 2**IrrFnc(ixyz)
  iSymCx = iSymX
  if (Ccoor(1) /= Zero) iSymCx = ibset(iSymCx,0)
  iSymCy = iSymY
  if (Ccoor(2) /= Zero) iSymCy = ibset(iSymCy,0)
  iSymCz = iSymZ
  if (Ccoor(3) /= Zero) iSymCz = ibset(iSymCz,0)

  ! Calculates QpV_ij = r_i*p_j+p_i*r_j (or rater p -> nabla)

  ! QpVxx
  OperI(1) = MltLbl(iSymCx,iSymX)
  OperC(1) = iChBas(2)
  ! QpVxy
  OperI(1+1) = MltLbl(iSymCx,iSymY)
  OperC(1+1) = iChBas(2)+iChBas(3)
  ! QpVxz
  OperI(1+2) = MltLbl(iSymCx,iSymZ)
  OperC(1+2) = iChBas(2)+iChBas(4)
  ! QpVyy
  OperI(1+3) = MltLbl(iSymCy,iSymY)
  OperC(1+3) = iChBas(3)
  ! QpVyz
  OperI(1+4) = MltLbl(iSymCy,iSymZ)
  OperC(1+4) = iChBas(3)+iChBas(4)
  ! QpVzz
  OperI(1+5) = MltLbl(iSymCz,iSymZ)
  OperC(1+5) = iChBas(4)

  call DCopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(QpVInt,QpVMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if   ! Vlct .and. (S%nMltpl > 2)
!***********************************************************************
!***********************************************************************
!12c)                                                                  *
!                                                                      *
!     Orbital Magnetic Quadrupole integrals.                           *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = -One
if (allocated(OMQ_Center) .and. (.not. Primitive_Pass)) then
  Label = 'OMQ     '
  nComp = 9
  nOrdOp = 3
  call Allocate_Auxiliary()

  call dcopy_(nComp,[OMQ_Center(1)],0,CoorO(1),3)
  call dcopy_(nComp,[OMQ_Center(2)],0,CoorO(1+1),3)
  call dcopy_(nComp,[OMQ_Center(3)],0,CoorO(1+2),3)
  Ccoor(:) = OMQ_Center(:)

  ixyz = 1
  iSymX = 2**IrrFnc(ixyz)
  ixyz = 2
  iSymY = 2**IrrFnc(ixyz)
  ixyz = 4
  iSymZ = 2**IrrFnc(ixyz)
  iSymCx = iSymX
  if (Ccoor(1) /= Zero) iSymCx = ibset(iSymCx,0)
  iSymCy = iSymY
  if (Ccoor(2) /= Zero) iSymCy = ibset(iSymCy,0)
  iSymCz = iSymZ
  if (Ccoor(3) /= Zero) iSymCz = ibset(iSymCz,0)

  iSymLx = ior(MltLbl(iSymCy,iSymZ),MltLbl(iSymCz,iSymY))
  iSymLy = ior(MltLbl(iSymCz,iSymX),MltLbl(iSymCx,iSymZ))
  iSymLz = ior(MltLbl(iSymCx,iSymY),MltLbl(iSymCy,iSymX))

  ! Calculates M_ij = r_j*L_i + L_i*r_j = 2*r_j*L_i + i hbar E_ijk r_k
  ! Since the i hbar is included outside we could do
  ! M_ij = r_j*L_i + L_i*r_j = 2*r_j*L_i + E_ijk r_k
  ! We use all nine components even if we only need 6 of these later.

  ! Mxx
  iSymxLx = MltLbl(iSymCx,iSymLx)
  iChOxx = iChBas(15)
  OperI(1) = iSymxLx
  OperC(1) = iChOxx
  ! Mxy
  iSymxLy = iSymCz
  iChOxy = iChBas(4)
  OperI(1+1) = iSymxLy
  OperC(1+1) = iChOxy
  ! Mxz
  iSymxLz = iSymCy
  iChOxz = iChBas(3)
  OperI(1+2) = iSymxLz
  OperC(1+2) = iChOxz
  ! Myx
  iSymyLx = iSymCz
  iChOyx = iChBas(4)
  OperI(1+3) = iSymyLx
  OperC(1+3) = iChOyx
  ! Myy
  iSymyLy = MltLbl(iSymCy,iSymLy)
  iChOyy = iChBas(15)
  OperI(1+4) = iSymyLy
  OperC(1+4) = iChOyy
  ! Myz
  iSymyLz = iSymCx
  iChOyz = iChBas(4)
  OperI(1+5) = iSymyLz
  OperC(1+5) = iChOyz
  ! Mzx
  iSymzLx = iSymCy
  iChOzx = iChBas(3)
  OperI(1+6) = iSymzLx
  OperC(1+6) = iChOzx
  ! Mzy
  iSymzLy = iSymCx
  iChOzy = iChBas(4)
  OperI(1+7) = iSymzLy
  OperC(1+7) = iChOzy
  ! Mzz
  iSymzLz = MltLbl(iSymCz,iSymLz)
  iChOzz = iChBas(15)
  OperI(1+8) = iSymzLz
  OperC(1+8) = iChOzz

  call DCopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(OMQInt,OMQMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if   ! OMQ_Center
!***********************************************************************
!***********************************************************************
!13)                                                                   *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (DKroll .and. Primitive_Pass) then
  rHrmt = One
  Label = 'pVp     '
  PLabel = 'NAInt '
  nOrdOp = 2
  nComp = 1
  call Allocate_Auxiliary()
  call dcopy_(3,[Zero],0,CoorO,1)
  OperI(1) = 1
  OperC(1) = iChBas(1)

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()

  if (BSS) then
    rHrmt = -One
    nOrdOp = 1
    nComp = 3
    call Allocate_Auxiliary()

    call dcopy_(3*nComp,[Zero],0,CoorO,1)
    call dcopy_(3,[Zero],0,Nuc,1)

    ixyz = 1
    OperI(1) = 2**IrrFnc(ixyz)
    OperC(1) = iChBas(2)
    ixyz = 2
    OperI(1+1) = 2**IrrFnc(ixyz)
    OperC(1+1) = iChBas(3)
    ixyz = 4
    OperI(1+2) = 2**IrrFnc(ixyz)
    OperC(1+2) = iChBas(4)

    Label = 'pV      '
    PLabel = 'NAInt '
    call OneEl(PXInt,PXMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

    Label = 'Vp      '
    call OneEl(VPInt,VPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

    call Deallocate_Auxiliary()
  end if    ! BSS
end if
!***********************************************************************
!***********************************************************************
!14)                                                                   *
!                                                                      *
!     Diamagnetic shielding integrals.                                 *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = One

mDMS = nDMS
if (Primitive_Pass) mDMS = 0

do iDMS=1,mDMS
  write(Label,'(A,I2,I2)') 'DMS ',1,iDMS
  nComp = 9
  nOrdOp = 2
  CCoor(:) = DMS_Centers(:,iDMS)
  call Allocate_Auxiliary()
  iSymC = 1
  if (Ccoor(1) /= Zero) iSymC = ior(iSymC,iSymX)
  if (Ccoor(2) /= Zero) iSymC = ior(iSymC,iSymY)
  if (Ccoor(3) /= Zero) iSymC = ior(iSymC,iSymZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero)) iSymC = ior(iSymC,iSymXY)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSymXZ)
  if ((Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSymYZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSyXYZ)

  iComp = 0
  iC = 0
  do ix=1,0,-1
    do iy=1-ix,0,-1
      iz = 1-ix-iy
      iC = iC+1
      iChO1 = iChBas(iC+1)
      ixyz = 0
      if (mod(ix,2) /= 0) ixyz = ibset(ixyz,0)
      if (mod(iy,2) /= 0) ixyz = ibset(ixyz,1)
      if (mod(iz,2) /= 0) ixyz = ibset(ixyz,2)
      iSym = 2**IrrFnc(ixyz)
      if (Ccoor(iC) /= Zero) iSym = ibset(iSym,0)
      iD = 0
      do jx=1,0,-1
        do jy=1-jx,0,-1
          jz = 1-jx-jy
          iD = iD+1
          iChO2 = iChBas(iD+1)
          jxyz = 0
          if (mod(jx,2) /= 0) jxyz = ibset(jxyz,0)
          if (mod(jy,2) /= 0) jxyz = ibset(jxyz,1)
          if (mod(jz,2) /= 0) jxyz = ibset(jxyz,2)
          iSymD = 2**IrrFnc(jxyz)
          if (Dxyz(iD) /= Zero) iSymD = ibset(iSymD,0)
          if (iC == iD) then
            i2 = iD+1
            if (i2 > 3) i2 = i2-3
            i3 = iD+2
            if (i3 > 3) i3 = i3-3
            iChO = iand(iChBas(i2+1),iChBas(i3+1))
          else
            iChO = ior(iChO1,iChO2)
          end if
          OperI(1+iComp) = MltLbl(iSymD,MltLbl(iSym,iSymC))
          OperC(1+iComp) = iChO
          call dcopy_(3,Ccoor,1,CoorO(1+iComp*3),1)
          iComp = iComp+1
        end do
      end do
    end do
  end do
  call dcopy_(3,Dxyz,1,CoorO(1+3),1)

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(DMSInt,DMSMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end do
!***********************************************************************
!***********************************************************************
!15)                                                                   *
!                                                                      *
!     Nuclear attraction integrals for finite centers.                 *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = One
Label = 'Center  '
!***********************************************************************
!***********************************************************************
!16)                                                                   *
!                                                                      *
!     Spherical well integrals.                                        *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = One
if ((.not. Prprt) .and. (.not. Primitive_Pass)) then
  nComp = 1
  iWel = 0
  call Allocate_Auxiliary()
  call dcopy_(3,[Zero],0,CoorO,1)
  OperI(1) = 1
  OperC(1) = iChBas(1)
  do iWel=1,nWel
    r0 = Wel_Info(1,iWel)
    ExpB = Wel_Info(2,iWel)
    write(Label,'(A,I4)') 'Well',iWel
    call OneEl(WelInt,WelMem,Label,ipList,OperI,nComp,CoorO,iWel,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  end do
  call Deallocate_Auxiliary()
end if  ! .not. Prprt
!***********************************************************************
!***********************************************************************
!5)                                                                    *
!                                                                      *
!     One-electron Hamiltonian integrals.                              *
!                                                                      *
!***********************************************************************
!***********************************************************************
if ((.not. Prprt) .and. (.not. Primitive_Pass)) then
  call mma_allocate(KnE_Int,n2Tri(1)+4,label='KnE_Int')
  call mma_allocate(NA_Int,n2Tri(1)+4,label='NA_Int')
# ifdef _FDE_
  ! Embedding
  if (embpot) call mma_allocate(Emb_Int,n2Tri(1)+4,label='Emb_Int')
# endif
  iOpt = 0
  lOper = 0
  iRC = -1
  iCmp = 1
  Label = 'Kinetic '
  call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_READ_)
  end if
  Label = 'Attract '
  lOper = 0
  iRC = -1
  call RdOne(iRC,iOpt,Label,iCmp,NA_Int,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_READ_)
  end if
  call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
# ifdef _FDE_
  ! Embedding
  if (embpot) then
    if (embPotInBasis) then
      !write(u6,*) 'ENTER'
      ! If the potential is given in basis set representation it
      ! has not been calculated with a OneEl call and is just read
      ! from file here.
      iunit = isFreeUnit(1)
      call molcas_open(iunit,embPotPath)
      do iEmb=1,n2Tri(1)
        read(iunit,*) Emb_Int(iEmb)
        !write(u6,*) iEmb-1,': ',Emb_Int(iEmb)
      end do
      close(iunit)
    else
      Label = 'embpot  '
      lOper = 0
      iRC = -1
      call RdOne(iRC,iOpt,Label,iCmp,Emb_Int,lOper)
      if (iRC /= 0) then
        call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
        call Quit(_RC_IO_ERROR_READ_)
      end if
    end if
    call DaXpY_(n2Tri(1)+4,One,Emb_Int,1,NA_Int,1)
  end if
# endif

  !--------Add contribution from ECP

  if (lECPnp) then
    Label = 'PrjInt  '
    lOper = 0
    iRC = -1
    call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
    Label = 'M1Int   '
    lOper = 0
    iRC = -1
    call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
    Label = 'M2Int   '
    lOper = 0
    iRC = -1
    call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
    Label = 'SROInt  '
    lOper = 0
    iRC = -1
    call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
  end if   ! lECPnp

  !--------Add contributions from the Pseudo Potential

  if (lPP) then
    Label = 'PPInt   '
    lOper = 0
    iRC = -1
    call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
  end if

  !--------Add contributions from the external field

  if (allocated(XF)) then
    Label = 'XFdInt  '
    lOper = 0
    iRC = -1
    call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
  end if ! XF

  !--------Add contributions from Spherical wells

  if (nWel /= 0) then
    do iWel=1,nWel
      Fact = Wel_Info(3,iWel)
      write(Label,'(A,I4)') 'Well',iWel
      lOper = 0
      iRC = -1
      call RdOne(iRC,iOpt,Label,iCmp,KnE_Int,lOper)
      if (iRC /= 0) then
        call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
        call Quit(_RC_IO_ERROR_READ_)
      end if
      call DaXpY_(n2Tri(1)+4,Fact,KnE_Int,1,NA_Int,1)
    end do
  end if  ! nWel /= 0

  Label = 'OneHam  '
  if (iPrint >= 10) call PrMtrx(Label,[lOper],1,[1],NA_Int)
  iRC = -1
  call WrOne(iRC,iOpt,Label,1,NA_Int,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error writing ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_WRITE_)
  end if

# ifdef _FDE_
  ! Embedding
  if (embpot) then
    Label = 'embpot  '
    iRC = -1
    call WrOne(iRC,iOpt,Label,1,embInt,lOper)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error writing ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_WRITE_)
    end if
  end if
# endif

  Label = 'OneHam 0'
  iRC = -1
  call WrOne(iRC,iOpt,Label,1,NA_Int,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error writing ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_WRITE_)
  end if
  call mma_deallocate(NA_Int)
  call mma_deallocate(KnE_Int)

# ifdef _FDE_
  ! Embedding
  if (embPot) call mma_deallocate(Emb_Int)
# endif
end if
!***********************************************************************
!***********************************************************************
!17)                                                                   *
!                                                                      *
!     Angular momentum products (PAM)                                  *
!                                                                      *
!***********************************************************************
!***********************************************************************
! Hermitized products of angular momentum integrals
! Component(1) is Lx*Lx
! Component(2) is (Lx*Ly+Ly*Lx)/2, etc.
! Coded P-A Malmqvist, Garching, Nov 1996
PLabel = ' '
rHrmt = -One
if (allocated(AMP_Center) .and. (.not. Primitive_Pass)) then
  Label = 'AMProd  '
  nComp = 6
  nOrdOp = 2
  call Allocate_Auxiliary()
  call dcopy_(nComp,[AMP_Center(1)],0,CoorO(1),3)
  call dcopy_(nComp,[AMP_Center(2)],0,CoorO(1+1),3)
  call dcopy_(nComp,[AMP_Center(3)],0,CoorO(1+2),3)
  CCoor(:) = AMP_Center(:)
  ! Symmetry labels iSymX  for operator d/dx, etc.
  ! Symmetry labels iSymLx for operator Lx, etc.
  ! Characters iChOx for operator Lx, etc.
  ixyz = 1
  iSymX = 2**IrrFnc(ixyz)
  ixyz = 2
  iSymY = 2**IrrFnc(ixyz)
  ixyz = 4
  iSymZ = 2**IrrFnc(ixyz)
  iSymCx = iSymX
  if (Ccoor(1) /= Zero) iSymCx = ibset(iSymCx,0)
  iSymCy = iSymY
  if (Ccoor(2) /= Zero) iSymCy = ibset(iSymCy,0)
  iSymCz = iSymZ
  if (Ccoor(3) /= Zero) iSymCz = ibset(iSymCz,0)
  iSymLx = ibset(MltLbl(iSymCy,iSymZ),MltLbl(iSymCz,iSymY))
  iChOx = iChBas(3)+iChBas(4)
  iSymLy = ibset(MltLbl(iSymCz,iSymX),MltLbl(iSymCx,iSymZ))
  iChOy = iChBas(4)+iChBas(2)
  iSymLz = ibset(MltLbl(iSymCx,iSymY),MltLbl(iSymCy,iSymX))
  iChOz = iChBas(2)+iChBas(3)

  ! Symmetry labels and characters of products. Let G be the full
  ! molecular point group, and Gsub=subgroup of G=stabilizer of
  ! gauge origin. The totally symmetric irrep of Gsub can be
  ! decomposed into irreps of G.
  ! Then symmetry label=packed array of bits, one for each irrep
  ! of G. The bit is set, if that irrep is included in the
  ! decomposition of the totally symmetric irrep of Gsub.
  OperI(1) = MltLbl(iSymLx,iSymLx)
  OperI(1+1) = MltLbl(iSymLx,iSymLy)
  OperI(1+2) = MltLbl(iSymLx,iSymLz)
  OperI(1+3) = MltLbl(iSymLy,iSymLy)
  OperI(1+4) = MltLbl(iSymLy,iSymLz)
  OperI(1+5) = MltLbl(iSymLz,iSymLz)
  OperC(1) = 0
  OperC(1+1) = ieor(iChOx,iChOy)
  OperC(1+2) = ieor(iChOx,iChOz)
  OperC(1+3) = 0
  OperC(1+4) = ieor(iChOy,iChOz)
  OperC(1+5) = 0

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(AMPInt,AMPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if
!***********************************************************************
!***********************************************************************
!17)                                                                   *
!                                                                      *
!     Contact term integrals                                           *
!                                                                      *
!***********************************************************************
!***********************************************************************
ixyz = 1
iSymX = 2**IrrFnc(ixyz)
ixyz = 2
iSymY = 2**IrrFnc(ixyz)
ixyz = 4
iSymZ = 2**IrrFnc(ixyz)
ixyz = 3
iSymXY = 2**IrrFnc(ixyz)
ixyz = 5
iSymXZ = 2**IrrFnc(ixyz)
ixyz = 6
iSymYZ = 2**IrrFnc(ixyz)
ixyz = 7
iSyXYZ = 2**IrrFnc(ixyz)
PLabel = ' '
rHrmt = One

mCnt = 0
if (nOrdEF == 2) mCnt = nEF

nComp = 1
nOrdOp = 1
do iCnt=1,mCnt
  write(Label,'(A,I5)') 'Cnt',iCnt
  CCoor(:) = EF_Centers(:,iCnt)
  call Allocate_Auxiliary()

  iSymR(0) = 1
  if (Ccoor(1) /= Zero) iSymR(0) = ior(iSymR(0),iSymX)
  if (Ccoor(2) /= Zero) iSymR(0) = ior(iSymR(0),iSymY)
  if (Ccoor(3) /= Zero) iSymR(0) = ior(iSymR(0),iSymZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero)) iSymR(0) = ior(iSymR(0),iSymXY)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(3) /= Zero)) iSymR(0) = ior(iSymR(0),iSymXZ)
  if ((Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymR(0) = ior(iSymR(0),iSymYZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymR(0) = ior(iSymR(0),iSyXYZ)

  OperI(1) = iSymR(0)
  OperC(1) = 0

  call dcopy_(nComp,[Zero],0,Nuc,1)
  call OneEl(CntInt,CntMem,Label,ipList,OperI,nComp,Ccoor,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! For picture-change corrected integrals.

  if (DKroll .and. Primitive_Pass) then
    write(Label,'(A,I2)') 'pCp   ',iCnt
    PLabel = 'CntInt'
    call FZero(Nuc,nComp)
    call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,CCoor,nOrdOp+2,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  end if

  call Deallocate_Auxiliary()
end do

! For a properties calculation, read the CNT values saved in a temp file
! and write the sum through Add_Info

if (PrPrt .and. (mCnt > 0)) then
  call mma_allocate(PtEl,1,label='PtEl')
  call mma_allocate(PtNuc,1,label='PtNuc')
  call mma_allocate(SumEl,1,label='SumEl')
  call mma_allocate(SumNuc,1,label='SumNuc')
  SumEl(1) = Zero
  SumNuc(1) = Zero
  !     Read and sum the values
  LuTmp = 10
  call DaName(LuTmp,'TMPPRP')
  iDisk = 0
  do iCnt=1,mCnt
    call dDaFile(LuTmp,2,PtEl,1,iDisk)
    call dDaFile(LuTmp,2,PtNuc,1,iDisk)
    SumEl(1) = SumEl(1)+PtEl(1)
    SumNuc(1) = SumNuc(1)+PtNuc(1)
  end do
  call DaClos(LuTmp)
  ! set the tolerance according to the total number of centers
  ! (assuming error scales with sqrt(mCnt))
  iTol = 5
  iTol = iTol-nint(Half*log10(real(mCnt,kind=wp)))
  write(label,'(a,a)') 'CNT','   el'
  call Add_Info(label,SumEl,1,iTol)
  write(label,'(a,a)') 'CNT','  nuc'
  call Add_Info(label,SumNuc,1,iTol)
  call mma_deallocate(PtEl)
  call mma_deallocate(PtNuc)
  call mma_deallocate(SumEl)
  call mma_deallocate(SumNuc)
end if
!***********************************************************************
!***********************************************************************
!18)                                                                   *
!                                                                      *
!     Gradient of overlap integrals with respect to the magnetic field *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (GIAO .and. (.not. Primitive_Pass)) then
  PLabel = ' '
  rHrmt = -One
  nB = 3
  iLow = 0
  mMltpl = 0 ! Do only overlap.
  !mMltpl = -1 ! Do only overlap.
  do iMltpl=iLow,mMltpl
    write(Label,'(A,I2)') 'dMP/dB',iMltpl
    mComp = (iMltpl+1)*(iMltpl+2)/2
    nComp = mComp*nB
    call DCopy_(3,Coor_MpM(1,iMltpl+1),1,Ccoor,1)
    call Allocate_Auxiliary()

    iComp = 0
    do ix=iMltpl,0,-1

      !---- Pick up which irrep each of the cartesian components of
      !     the operator belongs to. If the operator is associated
      !     with a center other than origin add the total symmetric
      !     irrep.

      if (mod(ix,2) == 0) then
        iSymX = 1
      else
        ixyz = 1
        iSymX = 2**IrrFnc(ixyz)
        if (Ccoor(1) /= Zero) iSymX = ibset(iSymX,0)
      end if
      do iy=iMltpl-ix,0,-1
        if (mod(iy,2) == 0) then
          iSymY = 1
        else
          ixyz = 2
          iSymY = 2**IrrFnc(ixyz)
          if (Ccoor(2) /= Zero) iSymY = ibset(iSymY,0)
        end if
        iz = iMltpl-ix-iy
        if (mod(iz,2) == 0) then
          iSymZ = 1
        else
          ixyz = 4
          iSymZ = 2**IrrFnc(ixyz)
          if (Ccoor(3) /= Zero) iSymZ = ibset(iSymZ,0)
        end if

        !----- Multiply cartesian components to generate which irreps
        !      the current element of the multipole moment operator
        !      belong to. The lowest significant bit if set indicate that
        !      a particular irrep is included.

        iTemp = MltLbl(iSymX,MltLbl(iSymY,iSymZ))

        iChO = mod(ix,2)*iChBas(2)+mod(iy,2)*iChBas(3)+mod(iz,2)*iChBas(4)

        !----- Now combine with the character of the first derivative
        !      with respect to the magnetic field.

        iSymRx = 2**IrrFnc(1)
        iSymRy = 2**IrrFnc(2)
        iSymRz = 2**IrrFnc(4)

        iB = 1
        iChOx = mod(ix,2)*iChBas(2)+mod(iy+1,2)*iChBas(3)+mod(iz+1,2)*iChBas(4)
        OperC(1+(iB-1)*mComp+iComp) = iChOx
        iSymBx = MltLbl(iSymRy,iSymRz)
        OperI(1+(iB-1)*mComp+iComp) = MltLbl(iTemp,iSymBx)
        call DCopy_(3,Coor_MPM(1,iMltpl+1),1,CoorO(1+((iB-1)*mComp+iComp)*3),1)

        iB = 2
        iChOy = mod(ix+1,2)*iChBas(2)+mod(iy,2)*iChBas(3)+mod(iz+1,2)*iChBas(4)
        OperC(1+(iB-1)*mComp+iComp) = iChOy
        iSymBy = MltLbl(iSymRz,iSymRx)
        OperI(1+(iB-1)*mComp+iComp) = MltLbl(iTemp,iSymBy)
        call DCopy_(3,Coor_MPM(1,iMltpl+1),1,CoorO(1+((iB-1)*mComp+iComp)*3),1)

        iB = 3
        iChOz = mod(ix+1,2)*iChBas(2)+mod(iy+1,2)*iChBas(3)+mod(iz,2)*iChBas(4)
        OperC(1+(iB-1)*mComp+iComp) = iChOz
        iSymBz = MltLbl(iSymRx,iSymRy)
        OperI(1+(iB-1)*mComp+iComp) = MltLbl(iTemp,iSymBz)
        call DCopy_(3,Coor_MPM(1,iMltpl+1),1,CoorO(1+((iB-1)*mComp+iComp)*3),1)

        iComp = iComp+1
      end do
    end do

    ! Zero nuclear contribution.
    call dcopy_(nComp,[Zero],0,Nuc,1)
    call OneEl(MltInt_GIAO,MltMem_GIAO,Label,ipList,OperI,nComp,CoorO,iMltpl,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

    call Deallocate_Auxiliary()

  end do

end if
!***********************************************************************
!***********************************************************************
!19)                                                                   *
!                                                                      *
!     Atomic mean field integrals for spin-orbit calculations          *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (lAMFI .and. (.not. Prprt) .and. (.not. Primitive_Pass)) then
  PLabel = ' '
  rHrmt = -One
  Label = 'AMFI    '
  nComp = 3
  call Allocate_Auxiliary()
  OperI(1) = 2**IrrFnc(6)
  OperC(1) = iChBas(7)
  OperI(1+1) = 2**IrrFnc(5)
  OperC(1+1) = iChBas(6)
  OperI(1+2) = 2**IrrFnc(3)
  OperC(1+2) = iChBas(4)

  ! BP - Turn off AMFI integrals for certain atom types
  !      as requested by the PAMF keyword
  !
  !      Note that Drv_AMFI will only process dbsc of valence type.
  !      Hence the restriction on the loop.

  call mma_allocate(iAtmNr2,nCnttp,Label='iAtmNr2')
  call mma_allocate(Charge2,nCnttp,Label='Charge2')
  do i=1,nCnttp
    iAtmNr2(i) = dbsc(i)%AtmNr
    Charge2(i) = dbsc(i)%Charge

    iAtom_Number = dbsc(i)%AtmNr
    if ((iAtom_Number < 0) .or. (iAtom_Number > size(No_AMFI))) then
      write(u6,*) 'Illegal atom number.'
      write(u6,*) 'Atom number=',iAtom_Number
      call Abend()
    end if
    if (No_AMFI(iAtom_Number)) then
      write(u6,*) 'Disabling AMFI for atom type ',dbsc(i)%AtmNr
      iAtmNr2(i) = 0
      Charge2(i) = Zero
    end if
  end do

  call Drv_AMFI(Label,OperI,nComp,iAtmNr2,Charge2)

  call mma_deallocate(iAtmNr2)
  call mma_deallocate(Charge2)
  call Deallocate_Auxiliary()
end if
!***********************************************************************
!***********************************************************************
!KAMAL,Rulin Update)                                                   *
!                                                                      *
!              GEN1INT                                                 *
!                                                                      *
!***********************************************************************
!***********************************************************************
!!!MXTC
if (lMXTC .and. DKroll .and. Primitive_Pass) then
# ifdef _GEN1INT_
  nOrdOp = 0
  ! Assume symmetric
  rHrmt = One
  nComp = 9
  call Get_nAtoms_All(nAtoms)
  do iCnt=1,nAtoms
    do jCnt=1,2
      if (jCnt == 1) then
        ! Label for lower triangular portion
        write(Label,'(A,I3)') 'MAGXP',iCnt
        write(PLabel,'(A6)') 'MagInt'
      else
        ! Label for upper triangular portion
        write(Label,'(A,I3)') 'MAGPX',iCnt
        write(PLabel,'(A6)') 'MagInt'
      end if
      call Allocate_Auxiliary()
      ! Dummy symmetry indices
      do i=1,nComp
        OperI(i) = 255
        OperC(i) = 0
      end do
      ! Zero nuclear contribution
      call dcopy_(nComp,[Zero],0,Nuc,1)
      ! Compute one electron integrals
      call OneEl(DumInt,DumMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
      call Deallocate_Auxiliary()
    end do
  end do
# else
  call WarningMessage(2,'Drv1El: NO Gen1int interface available!')
  call Abend()
# endif
end if ! lMXTC
!***********************************************************************
!***********************************************************************
!20)                                                                   *
!                                                                      *
!     The gradient of the kinetic energy and the nuclear attraction    *
!     energy with respect to the magnetic field.                       *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (GIAO .and. (.not. Primitive_Pass)) then
  PLabel = ' '
  rHrmt = -One
  nOrdOp = 0
  nComp = 3
  call Allocate_Auxiliary()
  call dcopy_(3*nComp,[Zero],0,CoorO,1)
  ixyz = 1
  OperI(1) = 2**IrrFnc(ixyz)
  OperC(1) = iChBas(2)
  ixyz = 2
  OperI(1+1) = 2**IrrFnc(ixyz)
  OperC(1+1) = iChBas(3)
  ixyz = 4
  OperI(1+2) = 2**IrrFnc(ixyz)
  OperC(1+2) = iChBas(4)

  call dcopy_(3,[Zero],0,Nuc,1)

  Label = 'dT/dB   '
  call OneEl(KneInt_GIAO,KneMem_GIAO,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  nOrdOp = 1
  Label = 'dV/dB   '
  call OneEl(NAInt_GIAO,NAMem_GIAO,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()

  ! Differentiate the generalized kinetic energy operator with respect
  ! to the magnetic moment at the centers.

  nOrdOp = 1
  nComp = 3
  ixyz = 1
  iSymX = 2**IrrFnc(ixyz)
  ixyz = 2
  iSymY = 2**IrrFnc(ixyz)
  ixyz = 4
  iSymZ = 2**IrrFnc(ixyz)
  ixyz = 3
  iSymXY = 2**IrrFnc(ixyz)
  ixyz = 5
  iSymXZ = 2**IrrFnc(ixyz)
  ixyz = 6
  iSymYZ = 2**IrrFnc(ixyz)
  ixyz = 7
  iSyXYZ = 2**IrrFnc(ixyz)

  iEF = 0
  do iCnttp=1,nCnttp
    do iCnt=1,dbsc(iCnttp)%nCntr
      iEF = iEF+1
      write(Label,'(A,I2)') 'dT/dmu',iEF
      CCoor(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
      call Allocate_Auxiliary()
      iSymC = 1
      if (Ccoor(1) /= Zero) iSymC = ior(iSymC,iSymX)
      if (Ccoor(2) /= Zero) iSymC = ior(iSymC,iSymY)
      if (Ccoor(3) /= Zero) iSymC = ior(iSymC,iSymZ)
      if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero)) iSymC = ior(iSymC,iSymXY)
      if ((Ccoor(1) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSymXZ)
      if ((Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSymYZ)
      if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSyXYZ)

      iComp = 0
      do ix=nOrdOp,0,-1
        do iy=nOrdOp-ix,0,-1
          iComp = iComp+1
          iz = nOrdOp-ix-iy
          ixyz = 0
          if (mod(ix,2) /= 0) ixyz = ibset(ixyz,0)
          if (mod(iy,2) /= 0) ixyz = ibset(ixyz,1)
          if (mod(iz,2) /= 0) ixyz = ibset(ixyz,2)
          iSym = 2**IrrFnc(ixyz)
          if (Ccoor(iComp) /= Zero) iSym = ibset(iSym,0)
          OperI(1+(iComp-1)) = MltLbl(iSymC,iSym)
          OperC(1+(iComp-1)) = iChBas(iComp+1)
          call dcopy_(3,Ccoor,1,CoorO(1+(iComp-1)*3),1)
        end do
      end do

      !call EFNuc(CoorO,Chrg,Centr,S%kCentr,Nuc,nOrdOp)

      call OneEl(dTdmu_Int,dTdmu_Mem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

      call Deallocate_Auxiliary()
    end do
  end do
end if
!***********************************************************************
!***********************************************************************
!21)                                                                   *
!                                                                      *
!     Atomic Fock matrix                                               *
!                                                                      *
!***********************************************************************
!***********************************************************************
PLabel = ' '
rHrmt = One
nComp = 1
nOrdOp = 0
if ((.not. Prprt) .and. (.not. Primitive_Pass) .and. Do_FckInt) then
  call Allocate_Auxiliary()
  call dcopy_(3,[Zero],0,CoorO,1)
  OperI(1) = 1
  OperC(1) = iChBas(1)

  Label = 'FckInt  '
  call Drv_Fck(Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)

  call Deallocate_Auxiliary()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory for property calculation.

if (Prprt) then
  if (Short) call mma_deallocate(Den)
  call mma_deallocate(Vec)
  call mma_deallocate(Occ)
  call CollapseOutput(0,'   Molecular properties:')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Is this the second pass in a relativistic calculation?
! In that case: close ONEINT, open ONEREL, read the pVp
! integrals, close ONEREL, re-open ONEINT, calculate the
! DK integrals in the contracted basis set, transform and
! replace the one-electron integrals on ONEINT.

if (DKroll .and. (.not. Primitive_Pass)) then
  if (BSS) then
    call BSSint()
  else
    call DKRelInt_DP()
  end if
end if
!***********************************************************************
!***********************************************************************
!20)                                                                   *
!                                                                      *
!     Produce information for the NEMO interface.                      *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (NEMO) then
  if (Primitive_Pass) then

    ! Compute p-matrix and put it temporarily on ONEREL.
    PLabel = ' '
    rHrmt = One
    nComp = 3
    nOrdOp = 0
    call Allocate_Auxiliary()
    do iComp=1,nComp
      call dcopy_(3,[Zero],0,CoorO(1+(iComp-1)*3),1)
      OperI(1+(iComp-1)) = 1
      OperC(1+(iComp-1)) = iChBas(1)
    end do

    Label = 'P_matrix'
    call OneEl(P_Int,P_Mem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero,Zero,Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
    call Deallocate_Auxiliary()
  else

    !-------Assemble transformation matrix between the contracted and
    !       the primitive basis.

    call NEMO_Opt1()

  end if
end if
!***********************************************************************
!***********************************************************************
!21)                                                                   *
!                                                                      *
!     Fragment AIEMP integrals: projection and 2-electron interaction  *
!                               integrals contracted with the          *
!                               fragment's density matrices            *
!                                                                      *
!***********************************************************************
!***********************************************************************
if (lFAIEMP .and. (.not. Primitive_Pass)) then
  ! projection integrals
  PLabel = ' '
  rHrmt = One
  nComp = 1
  nOrdOp = 0
  call Allocate_Auxiliary()
  call dcopy_(3,[Zero],0,CoorO,1)
  OperI(1) = 1
  OperC(1) = iChBas(1)
  Label = 'FragProj'
  call OneEl(FragPInt,FragPMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,[Zero],rHrmt,OperC,dum,1,dum,idum,0,0,dum,1,0)
  call Deallocate_Auxiliary()
  ! add the results to the one-electron hamiltonian
  iOpt = 0
  lOper = 0
  iRC = -1
  iCmp = 1
  call mma_allocate(FragP,n2Tri(1)+4,label='FragP')
  call RdOne(iRC,iOpt,Label,iCmp,FragP,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_READ_)
  end if
  Label = 'OneHam  '
  lOper = 0
  iRC = -1
  call mma_allocate(OneHam,n2Tri(1)+4,label='OneHam')
  call RdOne(iRC,iOpt,Label,iCmp,OneHam,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_READ_)
  end if
  call DaXpY_(n2Tri(1)+4,One,FragP,1,OneHam,1)
  iRC = -1
  call WrOne(iRC,iOpt,Label,1,OneHam,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error writing ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_WRITE_)
  end if
  iRC = -1
  Label = 'OneHam 0'
  call WrOne(iRC,iOpt,Label,1,OneHam,lOper)
  if (iRC /= 0) then
    call WarningMessage(2,'Drv1El: Error writing ONEINT;Label='//Label)
    call Quit(_RC_IO_ERROR_WRITE_)
  end if
  call mma_deallocate(FragP)
  call mma_deallocate(OneHam)
  ! 2-electron interaction integrals (are added to the one-electron
  ! hamiltonian locally)
  call Drv2El_FAIEMP()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
return

contains

subroutine Allocate_Auxiliary()
  implicit none

  call mma_Allocate(ipList,nComp,label='ipList')
  call mma_Allocate(OperI,nComp,label='OperI')
  call mma_Allocate(OperC,nComp,label='OperC')
  call mma_Allocate(CoorO,3*nComp,label='CoorO')
  call mma_Allocate(Nuc,nComp,label='Nuc')

  return
end subroutine Allocate_Auxiliary

subroutine Deallocate_Auxiliary()
  implicit none

  call mma_Deallocate(OperC)
  call mma_Deallocate(OperI)
  call mma_Deallocate(ipList)
  call mma_Deallocate(CoorO)
  call mma_Deallocate(Nuc)

  return
end subroutine Deallocate_Auxiliary

end subroutine Drv1el
