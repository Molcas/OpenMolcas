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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module Slapaf_Info

implicit none
private

public :: Cx, Gx, Gx0, NAC, Q_nuclear, dMass, Coor, Grd, ANr, Weights, Shift, GNrm, Lambda, Energy, Energy0, DipM, MF, qInt, &
          dqInt, nSup, Atom, RefGeo, BMx, Degen, jStab, nStab, iCoSet, BM, dBM, iBM, idBM, nqBM, R12, GradRef, KtB, AtomLbl, &
          Smmtrc, Lbl, mRowH, RootMap, Free_Slapaf, Get_Slapaf, Dmp_Slapaf, dqInt_Aux, BMx_kriging

! Arrays always allocated

real*8, allocatable :: Cx(:,:,:)     ! list of Cartesian coordinates
real*8, allocatable :: Gx(:,:,:)     ! list of Cartesian Gradients, State 1
real*8, allocatable :: Gx0(:,:,:)    ! list of Cartesian Gradients, State 2 for optimization of conical intersections
real*8, allocatable :: NAC(:,:,:)    ! list of Cartesian non-adiabatic coupling vectors
real*8, allocatable :: Q_nuclear(:)  ! list nuclear charges
real*8, allocatable :: dmass(:)      ! list atomic mass in units of (C=12)
real*8, allocatable :: Coor(:,:)     ! Cartesian coordinates of the last iteraction
real*8, allocatable :: Grd(:,:)      ! gradient of the last iteraction in Cartesian coordinates
real*8, allocatable :: Weights(:)    ! list of weights of ALL centers, however, the symmetry unique are first.
real*8, allocatable :: Shift(:,:)    ! list of displacements in Cartesian coordinates
real*8, allocatable :: GNrm(:)       ! list of the gradient norm for each iteration
real*8, allocatable :: Energy(:)     ! list of the energies of each iteration, State 1
real*8, allocatable :: Energy0(:)    ! list of the energies of each iteration, State 2 for optimization of conical intersections
real*8, allocatable :: MF(:,:)       ! list of Cartesian mode following vectors for each iteration
real*8, allocatable :: DipM(:,:)     ! list of dipole moments for each iteration
real*8, allocatable :: qInt(:,:)     ! internal coordinates for each iteration
real*8, allocatable :: dqInt(:,:)    ! derivatives of internal coordinates for each iteration
real*8, allocatable :: dqInt_Aux(:,:,:)      ! derivatives of internal coordinates for each iteration, of sets > 1
real*8, allocatable, target :: RefGeo(:,:)   ! Reference geometry in Cartesian coordinates
real*8, allocatable, target :: R12(:,:)      ! Reference geometry in R-P calculation (not used right now)
real*8, allocatable, target :: GradRef(:,:)  ! Reference gradient
real*8, allocatable, target :: Bmx(:,:)      ! the B matrix
real*8, allocatable :: BMx_kriging(:,:)      ! the updated B matrix during the Kriging procedure
real*8, allocatable :: Degen(:,:)    ! list of degeneracy numbers of the unique atoms (three identical entries)
integer, allocatable :: jStab(:,:), iCoset(:,:), nStab(:) ! stabilizer and cosets information for the indivudual centers
#include "LenIn.fh"
character(len=LenIn), allocatable :: AtomLbl(:) ! atomic labels
logical, allocatable :: Smmtrc(:,:)    ! Array with logical symmetry information on if a Cartesian is symmetric or not.
character(len=8), allocatable :: Lbl(:) ! Labels for internal coordinates and constraints

! Arrays for automatic internal coordinates
real*8, allocatable :: BM(:)         ! ...
real*8, allocatable :: dBM(:)        ! ...
integer, allocatable :: iBM(:)       ! ...
integer, allocatable :: idBM(:)      ! ...
integer, allocatable :: nqBM(:)      ! ...

integer, allocatable :: ANr(:)       ! list of atomic numbers

! Arrays optionally allocated

real*8, allocatable :: Lambda(:,:)   ! list of the Lagrange multipiers
integer, allocatable :: mRowH(:)     ! rows of the Hessian to be explicitly computed
integer, allocatable :: RootMap(:)   ! Array to map the roots between iterations

! Utility arrays with explicit deallocation, i.e. not via Free_Slapaf()

integer, allocatable :: Atom(:)      ! Temporary arrays for the super symmetry case
integer, allocatable :: NSup(:)      ! Temporary arrays for the super symmetry case
real*8, allocatable :: KtB(:,:)     ! KtB array for the BMtrx family of subroutines

logical :: Initiated = .false.
integer nsAtom

contains

subroutine Free_Slapaf()

# include "stdalloc.fh"

  if (allocated(Energy)) call mma_deallocate(Energy)
  if (allocated(Energy0)) call mma_deallocate(Energy0)
  if (allocated(DipM)) call mma_deallocate(DipM)
  if (allocated(GNrm)) call mma_deallocate(GNrm)
  if (allocated(Cx)) call mma_deallocate(Cx)
  if (allocated(Gx)) call mma_deallocate(Gx)
  if (allocated(Gx0)) call mma_deallocate(Gx0)
  if (allocated(NAC)) call mma_deallocate(NAC)
  if (allocated(MF)) call mma_deallocate(MF)
  if (allocated(Lambda)) call mma_deallocate(Lambda)
  if (allocated(Degen)) call mma_deallocate(Degen)
  if (allocated(jStab)) call mma_deallocate(jStab)
  if (allocated(iCoSet)) call mma_deallocate(iCoSet)
  if (allocated(nStab)) call mma_deallocate(nStab)
  if (allocated(AtomLbl)) call mma_deallocate(AtomLbl)
  if (allocated(Smmtrc)) call mma_deallocate(Smmtrc)
  if (allocated(Lbl)) call mma_deallocate(Lbl)

  if (allocated(Q_nuclear)) call mma_deallocate(Q_nuclear)
  if (allocated(dMass)) call mma_deallocate(dMass)
  if (allocated(Coor)) call mma_deallocate(Coor)
  if (allocated(Grd)) call mma_deallocate(Grd)
  if (allocated(ANr)) call mma_deallocate(ANr)
  if (allocated(Weights)) call mma_deallocate(Weights)
  if (allocated(Shift)) call mma_deallocate(Shift)
  if (allocated(BMx)) call mma_deallocate(BMx)
  if (allocated(BMx_kriging)) call mma_deallocate(BMx_kriging)

  if (allocated(BM)) call mma_deallocate(BM)
  if (allocated(dBM)) call mma_deallocate(dBM)
  if (allocated(iBM)) call mma_deallocate(iBM)
  if (allocated(idBM)) call mma_deallocate(idBM)
  if (allocated(nqBM)) call mma_deallocate(nqBM)

  if (allocated(RefGeo)) call mma_deallocate(RefGeo)
  if (allocated(R12)) call mma_deallocate(R12)
  if (allocated(R12)) call mma_deallocate(GradRef)

  if (allocated(qInt)) call mma_deallocate(qInt)
  if (allocated(dqInt)) call mma_deallocate(dqInt)
  if (allocated(dqInt_Aux)) call mma_deallocate(dqInt_Aux)
  if (allocated(mRowH)) call mma_deallocate(mRowH)
  if (allocated(RootMap)) call mma_deallocate(RootMap)
end subroutine Free_Slapaf

subroutine Get_Slapaf(iter,MaxItr,mTROld,lOld_Implicit,nsAtom_In,mLambda)

  use UnixInfo, only: SuperName

  integer iter, MaxItr, mTROld, nsAtom_In, mLambda
  logical lOld_Implicit
# include "real.fh"
# include "stdalloc.fh"
  logical Exist
  integer itmp, iOff, Lngth
  integer, allocatable :: Information(:)
  real*8, allocatable :: Relax(:)

  Initiated = .true.

  nsAtom = nsAtom_In

  call mma_allocate(Information,7,Label='Information')

  call qpg_iArray('Slapaf Info 1',Exist,itmp)
  if (Exist) call Get_iArray('Slapaf Info 1',Information,7)

  if ((.not. Exist) .or. (Information(1) == -99)) then
    !write(6,*) 'Reinitiate Slapaf fields on runfile'
    Information(:) = 0
    Information(3) = -99
    call Put_iArray('Slapaf Info 1',Information,7)
  end if

  iter = Information(2)+1
  if (iter >= MaxItr+1) then
    write(6,*) 'Increase MaxItr in slapaf_info.f90'
    call WarningMessage(2,'iter >= MaxItr+1')
    call Abend()
  end if
  mTROld = Information(3)
  lOld_Implicit = Information(4) == 1

  call mma_deallocate(Information)

  if (.not. allocated(Energy)) then

    call mma_allocate(Energy,MaxItr+1,Label='Energy')
    Energy(:) = Zero
    call mma_allocate(Energy0,MaxItr+1,Label='Energy0')
    Energy0(:) = Zero
    call mma_allocate(DipM,3,MaxItr+1,Label='DipM')
    DipM(:,:) = Zero
    call mma_allocate(GNrm,MaxItr+1,Label='GNrm')
    GNrm(:) = Zero
    call mma_allocate(Cx,3,nsAtom,MaxItr+1,Label='Cx')
    Cx(:,:,:) = Zero
    call mma_allocate(Gx,3,nsAtom,MaxItr+1,Label='Gx')
    Gx(:,:,:) = Zero
    call mma_allocate(Gx0,3,nsAtom,MaxItr+1,Label='Gx0')
    Gx0(:,:,:) = Zero
    call mma_allocate(NAC,3,nsAtom,MaxItr+1,Label='NAC')
    NAC(:,:,:) = Zero
    call mma_allocate(MF,3,nsAtom,Label='MF')
    MF(:,:) = Zero
    if (mLambda > 0) then
      call mma_allocate(Lambda,mLambda,MaxItr+1,Label='Lambda')
      Lambda(:,:) = Zero
    end if

  end if

  if (iter == 1) return

  if (SuperName /= 'numerical_gradient') then

    Lngth = size(Energy)+size(Energy0)+size(DipM)+size(GNrm)+size(Cx)+size(Gx)+size(Gx0)+size(NAC)+size(MF)
    if (allocated(Lambda)) then
      Lngth = Lngth+size(Lambda)
    end if
    call mma_allocate(Relax,Lngth,Label='Relax')
    call Get_dArray('Slapaf Info 2',Relax,Lngth)

    iOff = 1
    call DCopy_(size(Energy),Relax(iOff),1,Energy,1)
    iOff = iOff+size(Energy)
    call DCopy_(size(Energy0),Relax(iOff),1,Energy0,1)
    iOff = iOff+size(Energy0)
    call DCopy_(size(DipM),Relax(iOff),1,DipM,1)
    iOff = iOff+size(DipM)
    call DCopy_(size(GNrm),Relax(iOff),1,GNrm,1)
    iOff = iOff+size(GNrm)
    call DCopy_(size(Cx),Relax(iOff),1,Cx,1)
    iOff = iOff+size(Cx)
    call DCopy_(size(Gx),Relax(iOff),1,Gx,1)
    iOff = iOff+size(Gx)
    call DCopy_(size(Gx0),Relax(iOff),1,Gx0,1)
    iOff = iOff+size(Gx0)
    call DCopy_(size(NAC),Relax(iOff),1,NAC,1)
    iOff = iOff+size(NAC)
    call DCopy_(size(MF),Relax(iOff),1,MF,1)
    iOff = iOff+size(MF)
    if (allocated(Lambda)) then
      call DCopy_(size(Lambda),Relax(iOff),1,Lambda,1)
      iOff = iOff+size(Lambda)
    end if
    call mma_deallocate(Relax)
  else
    iter = 1
  end if

end subroutine Get_Slapaf

subroutine Dmp_Slapaf(stop,Just_Frequencies,Energy_In,Iter,MaxItr,mTROld,lOld_Implicit,nsAtom)

  use UnixInfo, only: SuperName

  logical stop, Just_Frequencies, lOld_Implicit
  real*8 Energy_In
  integer Iter, MaxItr, mTROld, nsAtom
# include "stdalloc.fh"
  integer, allocatable :: Information(:)
  real*8, allocatable :: Relax(:)
  real*8, allocatable :: GxFix(:,:)
  integer iOff_Iter, nSlap, iOff, Lngth
  logical Found

  if (.not. Initiated) then
    write(6,*) 'Dmp_Slapaf: Slapaf not initiated!'
    call Abend()
  else
    Initiated = .false.
  end if

  ! Write information of this iteration to the RLXITR file
  call mma_allocate(Information,7,Label='Information')
  if (stop) then
    Information(1) = -99     ! Deactivate the record
    iOff_Iter = 0
    call Put_iScalar('iOff_Iter',iOff_Iter)

    ! Restore the runfile data as if the computation was analytic
    ! (note the gradient sign must be changed back)

    if (Just_Frequencies) then
      call Put_dScalar('Last Energy',Energy_In)
      call mma_allocate(GxFix,3,nsAtom,Label='GxFix')
      call dcopy_(3*nsAtom,Gx,1,GxFix,1)
      GxFix(:,:) = -GxFix(:,:)
      call Put_dArray('GRAD',GxFix,3*nsAtom)
      call mma_deallocate(GxFix)
      call Put_dArray('Unique Coordinates',Cx,3*nsAtom)
      call Put_Coord_New(Cx,nsAtom)
    end if
  else
    call qpg_iArray('Slapaf Info 1',Found,nSlap)
    if (Found) then
      call Get_iArray('Slapaf Info 1',Information,7)
      if (Information(1) /= -99) Information(1) = MaxItr
    else
      Information(1) = MaxItr
    end if
  end if

  if (SuperName /= 'numerical_gradient') then
    Information(2) = Iter
    Information(3) = mTROld ! # symm. transl /rot.
    if (lOld_Implicit) then
      Information(4) = 1
    else
      Information(4) = 0
    end if
    Information(5) = 0
    Information(6) = size(Energy)+size(Energy0)+size(DipM)+size(GNrm)
    Information(7) = size(Energy)+size(Energy0)+size(DipM)+size(GNrm)+size(Cx)
    call Put_iArray('Slapaf Info 1',Information,7)

    Lngth = size(Energy)+size(Energy0)+size(DipM)+size(GNrm)+size(Cx)+size(Gx)+size(Gx0)+size(NAC)+size(MF)
    if (allocated(Lambda)) then
      Lngth = Lngth+size(Lambda)
    end if
    call mma_allocate(Relax,Lngth,Label='Relax')
    iOff = 1
    call DCopy_(size(Energy),Energy,1,Relax(iOff),1)
    iOff = iOff+size(Energy)
    call DCopy_(size(Energy0),Energy0,1,Relax(iOff),1)
    iOff = iOff+size(Energy0)
    call DCopy_(size(DipM),DipM,1,Relax(iOff),1)
    iOff = iOff+size(DipM)
    call DCopy_(size(GNrm),GNrm,1,Relax(iOff),1)
    iOff = iOff+size(GNrm)
    call DCopy_(size(Cx),Cx,1,Relax(iOff),1)
    iOff = iOff+size(Cx)
    call DCopy_(size(Gx),Gx,1,Relax(iOff),1)
    iOff = iOff+size(Gx)
    call DCopy_(size(Gx0),Gx0,1,Relax(iOff),1)
    iOff = iOff+size(Gx0)
    call DCopy_(size(NAC),NAC,1,Relax(iOff),1)
    iOff = iOff+size(NAC)
    call DCopy_(size(MF),MF,1,Relax(iOff),1)
    iOff = iOff+size(MF)
    if (allocated(Lambda)) then
      call DCopy_(size(Lambda),Lambda,1,Relax(iOff),1)
      iOff = iOff+size(Lambda)
    end if
    call Put_dArray('Slapaf Info 2',Relax,Lngth)
    call mma_deallocate(Relax)
  end if
  call mma_deallocate(Information)

end subroutine Dmp_Slapaf

end module Slapaf_Info
