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
! Copyright (C) 1999, Roland Lindh                                     *
!***********************************************************************

subroutine Setup_NQ(Maps2p,nShell,nSym,nNQ,Do_Grad,On_Top,Pck_Old,PMode_old,R_Min,nR_Min)
!***********************************************************************
!                                                                      *
! Object: to set up information for calculation of integrals via a     *
!         numerical quadrature.                                        *
! Warning: The exponents of each shell are reordered diffuse to compact*
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept of Chemical Physics,                                *
!             University of Lund, Sweden                               *
!             August 1999                                              *
!***********************************************************************

use iSD_data, only: iSD, nskal_iSD
use Basis_Info, only: dbsc, MolWgh, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iOper
use nq_Grid, only: Angular, Coor, Fact, Mem, nGridMax, nR_Eff, Pax
use nq_structure, only: NQ_Data, Open_NQ_Data
use nq_Info, only: Angular_Pruning, Block_size, Crowding, Fade, Functional_Type, GGA_type, L_Quad, LDA_type, MBC, meta_GGA_type1, &
                   meta_GGA_type2, mRad, nAngularGrids, nAtoms, ndc, nR, ntotgp, number_of_subblocks, nx, ny, nz, Off, On, &
                   Rotational_Invariance, T_Y, Threshold, x_min, y_min, z_min
use Grid_On_Disk, only: Final_Grid, G_S, Grid_Status, GridInfo, iDisk_Grid, iDisk_Set, iGrid_Set, Intermediate, Lu_Grid, &
                        Not_Specified, Old_Functional_Type, Regenerate, Use_Old
use Index_Functions, only: nTri_Elem1
use Pack_mod, only: isPack, PkThrs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nShell, nSym, nR_Min
integer(kind=iwp), intent(out) :: Maps2p(nShell,0:nSym-1), nNQ
logical(kind=iwp), intent(in) :: Do_Grad, On_Top
real(kind=wp), intent(out) :: Pck_Old, R_Min(0:nR_Min)
logical(kind=iwp), intent(out) :: PMode_old
#include "status.fh"
integer(kind=iwp) :: iAng, iAng_, iANr, iAt, iBas, iCar, iCmp, iCnt, iCnttp, iDCRR(0:7), iDrv, iDum(1), iIrrep, iNQ, iNQ_, &
                     iNQ_MBC, iPrim, iReset, iS, iSet, ish, iShell, iShll, iSym, iuv, kAO, lAng, lAngular, LmbdR, lSO, mAO, mdci, &
                     mdcj, mExp, nAngular, nCntrc, nDCRR, nDegi, nDegj, nDrv, nFOrd, nForm, nMem, nR_tmp, nRad, nRadial, NrExp, &
                     nSO, nTerm, nxyz
real(kind=wp) :: A_high, A_low, Alpha(2), Box_Size, C(3), Crowding_tmp, Dummy(1), dx, dy, dz, Fct, R_BS, rm(2), Threshold_tmp, &
                 ValExp, x_max, XYZ(3), y_max, z_max
logical(kind=iwp) :: PMode
real(kind=wp), allocatable :: Crd(:,:), dOdx(:,:,:,:), TempC(:,:), ZA(:)
real(kind=wp), external :: Bragg_Slater, Eval_RMin
logical(kind=iwp), external :: EQ

!                                                                      *
!***********************************************************************
!                                                                      *
Maps2p(:,:) = -99999999
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
!write(u6,*) '********** Setup_NQ ***********'
ntotgp = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Check if NQ environment has been activated

if ((NQ_Status /= Active) .and. (NQ_Status /= Inactive)) then
  call WarningMessage(2,'Setup_NQ: NQ_Status not initialized')
  call Quit_OnUserError()
end if
if (NQ_Status == Active) return
NQ_Status = Active
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the coordinates to the centers of all Voronoi polyhedra
!
! Note that this will be all centers with valence basis sets on
! them. Hence this will also include any pseudo centers!

call mma_allocate(TempC,3,nShell*nSym,Label='TempC')
nAtoms = 0
if (nShell > nskal_iSD) then
  write(u6,*) 'nShell > nSkal_iSD'
  write(u6,*) 'nShell=',nShell
  write(u6,*) 'nSkal_iSD=',nSkal_iSD
  call AbEnd()
end if
do iShell=1,nShell
  iCnttp = iSD(13,iShell)
  iCnt = iSD(14,iShell)
  XYZ(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  call Process_Coor(XYZ,TempC,nAtoms,nSym,iOper)
end do
call mma_allocate(Coor,3,nAtoms,Label='Coor')
Coor(:,:) = TempC(:,1:nAtoms)
call mma_deallocate(TempC)
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the symmetry unique coordinates

nNQ = nAtoms
call Open_NQ_Data(Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up the requested accuracy.

!Dr = -Log10(Thr)
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over each unique center, and find the highest and lowest
! Gaussians exponents associated with this center. The latter
! information will be used to design the radial grid associated
! with this center.

! Loop over the shells

do iShell=1,nShell

  ! Get the Atom number
  iANr = dbsc(iSD(13,iShell))%AtmNr

  iShll = iSD(0,iShell)  ! Get the angular momentum of ishell
  iAng = iSD(1,iShell)   ! Get the angular momentum of ishell
  nCntrc = iSD(3,iShell) ! Get the # of contracted functions for iShell
  mExp = iSD(5,iShell)   ! Get the number of exponents of ishell

  !*********************************************************************
  !                                                                    *
  ! Order the exponents diffuse to compact for the active shell        *
  !                                                                    *
  !*********************************************************************
  call OrdExpD2C(mExp,Shells(iShll)%Exp,nCntrc,Shells(iShll)%pCff)

  ! Get the extreme exponents for the active shell.
  A_low = Shells(iShll)%Exp(1)
  A_high = Shells(iShll)%Exp(mExp)

  iCnttp = iSD(13,iShell)
  iCnt = iSD(14,iShell)
  C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  outer: do iIrrep=0,nIrrep-1
    call OA(iOper(iIrrep),C,XYZ)
    do iNQ=1,nNQ

      if (EQ(NQ_data(iNQ)%Coor,XYZ)) then

        NQ_Data(iNQ)%Atom_Nr = iANR

        ! Assign the BS radius to the center
        NQ_Data(iNQ)%R_RS = Bragg_Slater(iANr)

        ! What is the maximum angular momentum for the active center ?
        NQ_Data(iNQ)%l_Max = max(NQ_Data(iNQ)%l_max,iAng)

        ! Get the extreme exponents for the atom
        NQ_Data(iNQ)%A_high = max(NQ_Data(iNQ)%A_high,A_High)
        NQ_Data(iNQ)%A_low = min(NQ_Data(iNQ)%A_low,A_low)

        Maps2p(iShell,iIrrep) = iNQ
        cycle outer
      end if
    end do ! iNQ
    call WarningMessage(2,'Did not find a center associated with the shell!')
    call Abend()
  end do outer  ! iIrrep

end do          ! iShell

!***********************************************************************
!                                                                      *
!               END OF THE LOOP OVER THE SHELLS.                       *
! Now we have the number of unique center and their associated         *
! exponents and maximum angular momentum. And for each shell the       *
! exponents are ordered diffuse to compact.                            *
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over all the atoms to create the radial quadrature              *
!                                                                      *
!***********************************************************************
!
! Allocate memory to store the number of effective radial points for
! each center and the radius of this center.
call mma_Allocate(nR_Eff,nNQ,Label='nR_Eff')

iNQ_MBC = 0
iReset = 0
Threshold_tmp = Zero
nR_tmp = 0
do iNQ=1,nNQ
  ! Get the extreme exponents for the atom
  Alpha(1) = NQ_Data(iNQ)%A_low
  Alpha(2) = NQ_Data(iNQ)%A_high

  ! Get the coordinates of the atom
  XYZ(:) = NQ_Data(iNQ)%Coor

  ! For a special center we can increase the accuracy.

  if (MBC /= ' ') then
    do iS=1,nShell
      iCnttp = iSD(13,iS)
      iCnt = iSD(14,iS)
      C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
      if (EQ(NQ_Data(iNQ)%Coor,C)) then
        mdci = iSD(10,iS)
        if (dc(mdci)%LblCnt == MBC) then
          nR_tmp = nR
          nR = int(real(nR,kind=wp)*Two)
          Threshold_tmp = Threshold
          Threshold = Threshold*1.0e-6_wp

          iReset = 1
          iNQ_MBC = iNQ
          exit
        end if
      end if
    end do
  end if

  ! Max angular momentum for the atom -> rm(1)
  ! Max Relative Error -> rm(2)
  rm(1) = real(NQ_Data(iNQ)%l_Max,kind=wp)
  rm(2) = Threshold

  call GenVoronoi(nR_Eff(iNQ),Alpha,rm,iNQ)

  if (iReset == 1) then
    nR = nR_tmp
    Threshold = Threshold_tmp
    iReset = 0
  end if

end do  ! iNQ
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the principal axis system and optionally to
! compute derivatives of the principal axis. Needed in order to
! compute the gradient of the rotationally invariant DFT energy.

call mma_allocate(Pax,3,3,Label='Pax')
call mma_allocate(dOdx,3,3,nNQ,3,Label='dOdx')
dOdx(:,:,:,:) = Zero
call mma_Allocate(ZA,nNQ,Label='ZA')
call mma_Allocate(Crd,3,nNQ)

! Collect coordinates and charges of the nuclei

do iNQ=1,nNQ
  ZA(iNQ) = real(NQ_Data(iNQ)%Atom_Nr,kind=wp)
  Crd(:,iNQ) = NQ_data(iNQ)%Coor
end do

call RotGrd(Crd,ZA,Pax,dOdx,Dummy,nNQ,Do_Grad,.false.)

! Distribute derivative of the principal axis system

if (Do_Grad) then
  do iNQ=1,nNQ
    call mma_allocate(NQ_Data(iNQ)%dOdx,3,3,3,Label='dOdx')
    do iCar=1,3
      NQ_Data(iNQ)%dOdx(:,:,iCar) = dOdx(:,:,iNQ,iCar)
    end do
  end do
end if

call mma_deallocate(dOdX)
call mma_deallocate(Crd)
call mma_deallocate(ZA)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Rotational_Invariance == Off) then
  Pax(:,:) = reshape([One,Zero,Zero,Zero,One,Zero,Zero,Zero,One],[3,3])
  do iNQ=1,nNQ
    if (.not. allocated(NQ_Data(iNQ)%dOdx)) call mma_allocate(NQ_Data(iNQ)%dOdx,3,3,3,Label='dOdx')
    NQ_Data(iNQ)%dOdx(:,:,:) = Zero
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the angular grid

call Angular_grid()

Crowding_tmp = Zero
do iNQ=1,nNQ

  ! Assign the angular grid to be used with each radial grid point

  call mma_allocate(NQ_Data(iNQ)%Angular,nR_Eff(iNQ),Label='Angular')
  NQ_Data(iNQ)%Angular(:) = nAngularGrids

  ! Prune the angular grid

  if (Angular_Pruning == On) then

    ! Find the R_min values of each angular shell

    lAng = NQ_Data(iNQ)%l_max
    do iAng=0,lAng
      R_Min(iAng) = Zero
      ValExp = -One
      iSet = 0
      do iShell=1,nShell
        iShll = iSD(0,iShell)
        iAng_ = iSD(1,iShell)
        NrExp = iSD(5,iShell)
        !write(u6,*) 'iAng_,iAng=',iAng_,iAng
        if ((iAng_ == iAng) .and. (NrExp >= 1)) then
          do iSym=0,nSym-1
            iNQ_ = Maps2p(iShell,iSym)
            !write(u6,*) 'iNQ_,iNQ=',iNQ_,iNQ
            if (iNQ_ == iNQ) then
              ValExp = Shells(iShll)%Exp(NrExp)
              iSet = 1
            end if
          end do
        end if
      end do
      if ((ValExp < Zero) .and. (iSet == 1)) then
        call WarningMessage(2,'ValExp < Zero')
        call Abend()
      end if
      if (iSet == 1) then
        R_Min(iAng) = Eval_RMin(ValExp,iAng,Threshold)
        if (iAng == 0) R_Min(iAng) = Zero
      end if
    end do

    R_BS = NQ_Data(iNQ)%R_RS

    if (iNQ == iNQ_MBC) then
      Crowding_tmp = Crowding
      Crowding = One+(Crowding-One)*Quart
      iReset = 1
    end if

    call Angular_Prune(NQ_Data(iNQ)%R_Quad,nR_Eff(iNQ),NQ_Data(iNQ)%Angular,Crowding,Fade,R_BS,L_Quad,R_Min,lAng,nAngularGrids)

    if (iReset == 1) then
      Crowding = Crowding_tmp
      iReset = 0
    end if

  end if

end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,'(A)') ' =================================='
write(u6,'(A)') ' =        Grid information        ='
write(u6,'(A)') ' =================================='
write(u6,'(A)') ' Legend             '
write(u6,'(A)') ' ----------------------------------'
write(u6,'(A)') ' ANr: element number'
write(u6,'(A)') ' nR : number of radial grid points'
write(u6,'(A)') ' iNQ: grid index'
write(u6,'(A)') ' ----------------------------------'
write(u6,*)
write(u6,'(A)') ' iNQ ANr  nR'
do iNQ=1,nNQ
  iANr = NQ_Data(iNQ)%Atom_Nr
  kR = nR_Eff(iNQ)
  write(u6,'(3I4)') iNQ,iANr,kR
end do
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Determine the spatial extension of the molecular system

!Box_Size = Four      ! Angstrom
Box_Size = Two        ! Angstrom
!Box_Size = Half ! Angstrom
Block_size = Box_Size
x_min = huge(x_min)
y_min = huge(y_min)
z_min = huge(z_min)
x_max = -huge(x_max)
y_max = -huge(y_max)
z_max = -huge(z_max)
do iAt=1,nAtoms
  x_min = min(x_min,Coor(1,iAt))
  y_min = min(y_min,Coor(2,iAt))
  z_min = min(z_min,Coor(3,iAt))
  x_max = max(x_max,Coor(1,iAt))
  y_max = max(y_max,Coor(2,iAt))
  z_max = max(z_max,Coor(3,iAt))
end do

! Add half a box size around the whole molecule

x_min = x_min-Box_Size/Two
y_min = y_min-Box_Size/Two
z_min = z_min-Box_Size/Two
x_max = x_max+Box_Size/Two
y_max = y_max+Box_Size/Two
z_max = z_max+Box_Size/Two

! At least one finite box. Adjust to an even number of boxes.

nx = int((x_max-x_min+Box_Size)/Box_Size)
nx = 2*((nx+1)/2)
ny = int((y_max-y_min+Box_Size)/Box_Size)
ny = 2*((ny+1)/2)
nz = int((z_max-z_min+Box_Size)/Box_Size)
nz = 2*((nz+1)/2)

! Adjust extremal values to fit exactly with the box size.

dx = (real(nx,kind=wp)*Box_Size-(x_max-x_min))/Two
dy = (real(ny,kind=wp)*Box_Size-(y_max-y_min))/Two
dz = (real(nz,kind=wp)*Box_Size-(z_max-z_min))/Two

x_min = x_min-dx
y_min = y_min-dy
z_min = z_min-dz
x_max = x_max+dx
y_max = y_max+dy
z_max = z_max+dz

! Add the infinite edge boxes

nx = nx+2
ny = ny+2
nz = nz+2
#ifdef _DEBUGPRINT_
write(u6,*) 'x_min=',x_min,dx
write(u6,*) 'y_min=',y_min,dy
write(u6,*) 'z_min=',z_min,dz
write(u6,*) 'x_max=',x_max
write(u6,*) 'y_max=',y_max
write(u6,*) 'z_max=',z_max
write(u6,*) 'nx,ny,nz=',nx,ny,nz
write(u6,*) 'Total number of blocks=',nx*ny*nz
#endif
number_of_subblocks = nx*ny*nz
!                                                                      *
!***********************************************************************
!                                                                      *
! nFOrd: the order of the functional. nFOrd-1 is the number of times
!        the basis functions has to be differentiated to compute the
!        energy contribution.
! mRad: number of different radial functions associated with a
!       basis function. This number depends on the type of
!       functional and the number of times the basis function has
!       to be differentiated in order to produce the values of the
!       parameters which the functional depends on (rho, grad rho,
!       and nabla rho).
! mAO: number of elements a basis function generates upon
!      differentiation (1,4,10,20, etc.)

select case (Functional_type)

  case (LDA_type)
    nFOrd = 1
    mAO = (nFOrd*(nFOrd+1)*(nFOrd+2))/6
    if (do_grad) mAO = 4 !AMS - GRADIENTS?
    if (.not. Do_Grad) then
      mRad = nFOrd
    else
      mRad = nFOrd+1
    end if

  case (GGA_type)
    nFOrd = 2
    mAO = (nFOrd*(nFOrd+1)*(nFOrd+2))/6
    if (do_grad) mAO = 10
    if (.not. Do_Grad) then
      mRad = nFOrd
    else
      mRad = nFOrd+1
    end if

  case (meta_GGA_type1)
    nFOrd = 2
    mAO = (nFOrd*(nFOrd+1)*(nFOrd+2))/6
    if (.not. Do_Grad) then
      mRad = nFOrd
    else
      mRad = NFOrd+1
    end if

  case (meta_GGA_type2)
    nFOrd = 3
    mAO = (nFOrd*(nFOrd+1)*(nFOrd+2))/6
    if (.not. Do_Grad) then
      mRad = nFOrd
    else
      mRad = NFOrd+1
    end if

  case default
    mRad = 0 ! Dummy initialize
    mAO = 0  ! Dummy initialize
    call WarningMessage(2,'Functional_type == Other_type')
    call Abend()
end select
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate scratch for processing AO's on the grid

nMem = 0
nSO = 0
lSO = 0
lAngular = 0
do ish=1,nShell
  iAng = iSD(1,iSh)
  iCmp = iSD(2,iSh)
  iBas = iSD(3,iSh)
  iPrim = iSD(5,iSh)

  nxyz = nGridMax*3*(iAng+mRad)
  nDrv = mRad-1
  nForm = 0
  do iDrv=0,nDrv
    nForm = nForm+nTri_Elem1(iDrv)
  end do
  nTerm = 2**nDrv
  nAngular = 5*nForm*nTerm
  nRad = iPrim*nGridMax*mRad
  nRadial = iBas*nGridMax*mRad
  if (On_Top) then
    mdci = iSD(10,iSh)
    kAO = iCmp*iBas*nGridMax
    nSO = kAO*nSym/dc(mdci)%nStab*mAO
  end if
  nMem = max(nMem,nxyz+nRad+nRadial)
  lSO = max(lSO,nSO)
  lAngular = max(lAngular,nAngular)
end do

call mma_allocate(Angular,lAngular,Label='Angular')
call mma_allocate(Mem,nMem,Label='Mem')
!                                                                      *
!***********************************************************************
!                                                                      *
! Access the file with Grid points and weights.

! Open the file.
Lu_Grid = 88
call DaName_MF_WA(Lu_Grid,'NQGRID')

if (iGrid_Set == Not_Specified) iGrid_Set = Final_Grid

! Read the status flag.
iDisk_Grid = 0
call iDaFile(Lu_Grid,2,G_S,2,iDisk_Grid)
call iDaFile(Lu_Grid,2,iDisk_Set,2,iDisk_Grid)
call iDaFile(Lu_Grid,2,iDum,1,iDisk_Grid)
Old_Functional_Type = iDum(1)

Grid_Status = G_S(iGrid_Set)
if (Old_Functional_Type /= Functional_Type) then
  G_S(Final_Grid) = Regenerate
  G_S(Intermediate) = Regenerate
  Grid_Status = Regenerate
end if
iDisk_Grid = iDisk_Set(iGrid_Set)

! Allocate memory for the master TOC.
call mma_Allocate(GridInfo,2,number_of_subblocks,Label='GridInfo')

! Retrieve the TOC or regenerate it.

! The table contains two data items per subblock.
! 1) disk address and 2) number of batches.

if (Grid_Status == Regenerate) then
  !write(u6,*) 'Grid_Status == Regenerate'
  Grid_Status = Regenerate
  GridInfo(:,:) = 0
  call iDaFile(Lu_Grid,1,GridInfo,2*number_of_subblocks,iDisk_Grid)
  Old_Functional_Type = Functional_Type
else if (Grid_Status == Use_Old) then
  !write(u6,*) 'Grid_Status == Use_Old'
  call iDaFile(Lu_Grid,2,GridInfo,2*number_of_subblocks,iDisk_Grid)
else
  call WarningMessage(2,'Illegal Grid Status!')
  call Abend()
end if

Pck_Old = PkThrs
PMode_old = isPack
PMode = .true.
call IniPkR8(T_Y,PMode)
!                                                                      *
!***********************************************************************
!                                                                      *
! Setup some symmetry stuff outside the loop

ndc = 0
do iSh=1,nShell
  ndc = max(ndc,iSD(10,iSh))
end do
call mma_allocate(Fact,ndc,ndc,Label='Fact')
do mdci=1,ndc
  nDegi = nIrrep/dc(mdci)%nStab
  do mdcj=1,ndc
    nDegj = nIrrep/dc(mdcj)%nStab

    call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)

    iuv = dc(mdci)%nStab*dc(mdcj)%nStab
    if (MolWgh == 1) then
      Fct = real(nIrrep,kind=wp)/real(LmbdR,kind=wp)
    else if (MolWgh == 0) then
      Fct = real(iuv,kind=wp)/real(nIrrep*LmbdR,kind=wp)
    else
      Fct = sqrt(real(iuv,kind=wp))/real(LmbdR,kind=wp)
    end if
    Fct = Fct*real(nDCRR,kind=wp)/real(nDegi*nDegj,kind=wp)

    ! Save: Fact

    Fact(mdci,mdcj) = Fct

  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Setup_NQ
