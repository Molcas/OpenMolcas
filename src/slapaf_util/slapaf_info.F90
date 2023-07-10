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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

! Arrays always allocated
!
! Cx                   list of Cartesian coordinates
! Gx                   list of Cartesian Gradients, State 1
! Gx0                  list of Cartesian Gradients, State 2 for optimization of conical intersections
! NAC                  list of Cartesian non-adiabatic coupling vectors
! Q_nuclear            list nuclear charges
! dmass                list atomic mass in units of (C=12)
! Coor                 Cartesian coordinates of the last iteration
! Grd                  gradient of the last iteraction in Cartesian coordinates
! Weights              list of weights of ALL centers, however, the symmetry unique are first.
! Shift                list of displacements in Cartesian coordinates
! GNrm                 list of the gradient norm for each iteration
! Energy               list of the energies of each iteration, State 1
! Energy0              list of the energies of each iteration, State 2 for optimization of conical intersections
! MF                   list of Cartesian mode following vectors for each iteration
! DipM                 list of dipole moments for each iteration
! qInt                 internal coordinates for each iteration
! dqInt                derivatives of internal coordinates for each iteration
! dqInt_Aux            derivatives of internal coordinates for each iteration, of sets > 1
! RefGeo               Reference geometry in Cartesian coordinates
! R12                  Reference geometry in R-P calculation (not used right now)
! GradRef              Reference gradient
! Bmx                  the B matrix
! BMx_kriging          the updated B matrix during the Kriging procedure
! Degen                list of degeneracy numbers of the unique atoms (three identical entries)
! Stab, iCoset, nStab  stabilizer and cosets information for the indivudual centers
! AtomLbl              atomic labels
! Smmtrc               Array with logical symmetry information on if a Cartesian is symmetric or not.
! Lbl                  Labels for internal coordinates and constraints
!Arrays for automatic internal coordinates
! BM                   ...
! dBM                  ...
! iBM                  ...
! idBM                 ...
! nqBM                 ...
! ANr                  list of atomic numbers
!
! Arrays optionally allocated
!
! Lambda               list of the Lagrange multipiers
! mRowH                rows of the Hessian to be explicitly computed
! RootMap              Array to map the roots between iterations
!
! Utility arrays with explicit deallocation, i.e. not via Free_Slapaf()
!
! Atom                  Temporary arrays for the super symmetry case
! NSup                  Temporary arrays for the super symmetry case
! KtB                   KtB array for the BMtrx family of subroutines

#include "LenIn.fh"
logical(kind=iwp) :: Initiated = .false.
integer(kind=iwp), allocatable :: ANr(:), Atom(:), iBM(:), idBM(:), iCoset(:,:), jStab(:,:), mRowH(:), nqBM(:), nStab(:), NSup(:), &
                                  RootMap(:)
real(kind=wp), allocatable :: BM(:), BMx_kriging(:,:), Coor(:,:), Cx(:,:,:), dBM(:), Degen(:,:), DipM(:,:), dmass(:), dqInt(:,:), &
                              dqInt_Aux(:,:,:), Energy(:), Energy0(:), GNrm(:), Grd(:,:), Gx(:,:,:), Gx0(:,:,:), KtB(:,:), &
                              Lambda(:,:), MF(:,:), NAC(:,:,:), Q_nuclear(:), qInt(:,:), Shift(:,:), Weights(:)
real(kind=wp), allocatable, target :: Bmx(:,:), GradRef(:,:), R12(:,:), RefGeo(:,:)
logical(kind=iwp), allocatable :: Smmtrc(:,:)
character(len=LenIn), allocatable :: AtomLbl(:)
character(len=8), allocatable :: Lbl(:)

public :: ANr, Atom, AtomLbl, BM, Bmx, BMx_kriging, Coor, Cx, dBM, Degen, DipM, dmass, Dmp_Slapaf, dqInt, dqInt_Aux, Energy, &
          Energy0, Free_Slapaf, Get_Slapaf, GNrm, GradRef, Grd, Gx, Gx0, iBM, iCoset, idBM, jStab, KtB, Lambda, Lbl, MF, mRowH, &
          NAC, nqBM, nStab, NSup, Q_nuclear, qInt, R12, RefGeo, RootMap, Shift, Smmtrc, Weights

contains

subroutine Free_Slapaf()

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

subroutine Get_Slapaf(iter,MaxItr,mTROld,lOld_Implicit,nsAtom,mLambda)

  use UnixInfo, only: SuperName
  use Constants, only: Zero
  use Definitions, only: u6

  integer(kind=iwp) :: iter, MaxItr, mTROld, nsAtom, mLambda
  logical(kind=iwp) :: lOld_Implicit
  integer(kind=iwp) :: iLn, iOff, itmp, Lngth
  logical(kind=iwp) :: Exists
  integer(kind=iwp), allocatable :: Information(:)
  real(kind=wp), allocatable :: Relax(:)

  Initiated = .true.

  call mma_allocate(Information,7,Label='Information')

  call qpg_iArray('Slapaf Info 1',Exists,itmp)
  if (Exists) call Get_iArray('Slapaf Info 1',Information,7)

  if ((.not. Exists) .or. (Information(1) == -99)) then
    !write(u6,*) 'Reinitiate Slapaf fields on runfile'
    Information(:) = 0
    Information(3) = -99
    call Put_iArray('Slapaf Info 1',Information,7)
  end if

  iter = Information(2)+1
  if (iter >= MaxItr+1) then
    write(u6,*) 'Increase MaxItr in slapaf_info.f90'
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
    if (allocated(Lambda)) Lngth = Lngth+size(Lambda)
    call mma_allocate(Relax,Lngth,Label='Relax')
    call Get_dArray('Slapaf Info 2',Relax,Lngth)

    iOff = 1
    iLn = size(Energy)
    Energy(:) = Relax(iOff:iOff+iLn-1)
    iOff = iOff+iLn
    iLn = size(Energy0)
    Energy0(:) = Relax(iOff:iOff+iLn-1)
    iOff = iOff+iLn
    iLn = size(DipM)
    DipM(:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(DipM))
    iOff = iOff+iLn
    iLn = size(GNrm)
    GNrm(:) = Relax(iOff:iOff+iLn-1)
    iOff = iOff+iLn
    iLn = size(Cx)
    Cx(:,:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(Cx))
    iOff = iOff+iLn
    iLn = size(Gx)
    Gx(:,:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(Gx))
    iOff = iOff+iLn
    iLn = size(Gx0)
    Gx0(:,:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(Gx0))
    iOff = iOff+iLn
    iLn = size(NAC)
    NAC(:,:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(NAC))
    iOff = iOff+iLn
    iLn = size(MF)
    MF(:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(MF))
    iOff = iOff+iLn
    if (allocated(Lambda)) then
      iLn = size(Lambda)
      Lambda(:,:) = reshape(Relax(iOff:iOff+iLn-1),shape(Lambda))
      iOff = iOff+iLn
    end if
    call mma_deallocate(Relax)
  else
    iter = 1
  end if

end subroutine Get_Slapaf

subroutine Dmp_Slapaf(lStop,Just_Frequencies,Energy_In,Iter,MaxItr,mTROld,lOld_Implicit,nsAtom)

  use UnixInfo, only: SuperName
  use Definitions, only: u6

  logical(kind=iwp) :: lStop, Just_Frequencies, lOld_Implicit
  real(kind=wp) :: Energy_In
  integer(kind=iwp) :: Iter, MaxItr, mTROld, nsAtom
  integer(kind=iwp) :: iLn, iOff, iOff_Iter, Lngth, nSlap
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: Information(:)
  real(kind=wp), allocatable :: GxFix(:,:), Relax(:)

  if (.not. Initiated) then
    write(u6,*) 'Dmp_Slapaf: Slapaf not initiated!'
    call Abend()
  else
    Initiated = .false.
  end if

  ! Write information of this iteration to the RLXITR file
  call mma_allocate(Information,7,Label='Information')
  if (lStop) then
    Information(1) = -99     ! Deactivate the record
    iOff_Iter = 0
    call Put_iScalar('iOff_Iter',iOff_Iter)

    ! Restore the runfile data as if the computation was analytic
    ! (note the gradient sign must be changed back)

    if (Just_Frequencies) then
      call Put_dScalar('Last Energy',Energy_In)
      call mma_allocate(GxFix,3,nsAtom,Label='GxFix')
      GxFix(:,:) = -Gx(:,:,1)
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
    if (allocated(Lambda)) Lngth = Lngth+size(Lambda)
    call mma_allocate(Relax,Lngth,Label='Relax')
    iOff = 1
    iLn = size(Energy)
    Relax(iOff:iOff+iLn-1) = pack(Energy,.true.)
    iOff = iOff+iLn
    iLn = size(Energy0)
    Relax(iOff:iOff+iLn-1) = pack(Energy0,.true.)
    iOff = iOff+iLn
    iLn = size(DipM)
    Relax(iOff:iOff+iLn-1) = pack(DipM,.true.)
    iOff = iOff+iLn
    iLn = size(GNrm)
    Relax(iOff:iOff+iLn-1) = pack(GNrm,.true.)
    iOff = iOff+iLn
    iLn = size(Cx)
    Relax(iOff:iOff+iLn-1) = pack(Cx,.true.)
    iOff = iOff+iLn
    iLn = size(Gx)
    Relax(iOff:iOff+iLn-1) = pack(Gx,.true.)
    iOff = iOff+iLn
    iLn = size(Gx0)
    Relax(iOff:iOff+iLn-1) = pack(Gx0,.true.)
    iOff = iOff+iLn
    iLn = size(NAC)
    Relax(iOff:iOff+iLn-1) = pack(NAC,.true.)
    iOff = iOff+iLn
    iLn = size(MF)
    Relax(iOff:iOff+iLn-1) = pack(MF,.true.)
    iOff = iOff+iLn
    if (allocated(Lambda)) then
      iLn = size(Lambda)
      Relax(iOff:iOff+iLn-1) = pack(Lambda,.true.)
      iOff = iOff+iLn
    end if
    call Put_dArray('Slapaf Info 2',Relax,Lngth)
    call mma_deallocate(Relax)
  end if
  call mma_deallocate(Information)

end subroutine Dmp_Slapaf

end module Slapaf_Info
