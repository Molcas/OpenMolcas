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
use Constants, only: Zero, One, Three, Half
use Definitions, only: wp, iwp

implicit none
private

! Baker      Convergence a la Baker
! Beta       The threshold for restricted step optimization.
! Beta_Disp  The threshold for restricted variance optimization.
!
!     Hessian update
! 0   iOptH=0000001 ( 1) Meyer (disabled)
! 1   iOptH=0000010 ( 2) BP (disabled)
! 2   iOptH=0000100 ( 4) BFGS
! 3   iOptH=0001000 ( 8) None
! 4   iOptH=0010000 (16) MPS, for TS search
! 5   iOptH=0100000 (32) EU, for TS search
! 6   iOptH=1000000 (64) TS-BFGS, for TS search
!
! Optimization method. DO NOT EVER GO BEYOND BIT 30!!!
!
!     iOptC=00000000000000  (  0) No optimization
!  0  iOptC=00000000000001  (  1) Quasi Newton-Raphson
!  1  iOptC=00000000000010  (  2) c1-DIIS
!  2  iOptC=00000000000100  (  4) c2-DIIS
!  3  iOptC=00000000001000  (  8) RS-RFO
!  4  iOptC=0000000001....  ( 16) DIIS, <dx|dx>
!  5  iOptC=0000000010....  ( 32) DIIS, <dx|g>
!  6  iOptC=0000000100....  ( 64) DIIS, <g|g>
!  7  iOptC=0000001.......  (128) Minimum, if not set TS search
!  8  iOptC=0000010.......  (256) Optimization with constraint
!  9  iOptC=00001.........  (512) set: RS-I-RFO, unset: RS-P-RFO
! 10  iOptC=0001.......... (1024) HMF augmented with weak interactions
! 11  iOptC=0010.......... (2048) augmented HMF used for selection of internal coordinates
! 12  iOptC=01............ (4096) set if FindTS
! 13  iOptC=10............ (8192) set if FindTS and in TS regime

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
integer(kind=iwp), parameter :: MaxItr = 2000
integer(kind=iwp) :: iInt = 0, iNeg(2) = 0, iOptC = int(b'111011001000'), iOptH = int(b'100'), IRC = 0, iRef = 0, iRow = 0, &
                     iRow_c = 0, iState(2) = 0, iter = 0, Max_Center = 15, mB_Tot = 0, mdB_Tot = 0, MEPNum = 0, Mode = -1, &
                     mq = 0, mTROld = 0, mTtAtm = 0, MxItr = 0, nBVec = 0, nDimBC = 0, nFix = 0, nLambda = 0, nMEP = MaxItr, &
                     NmIter = 0, nsRot = 0, nUserPT = 0, nWndw = 5
real(kind=wp) :: Beta = 0.3_wp, Beta_Disp = 0.3_wp, CnstWght = One, Delta = 1.0e-2_wp, dMEPStep = 0.1_wp, E_Delta = Zero, &
                 GNrm_Threshold = 0.2_wp, GrdMax = Zero, StpMax = Zero, rFuzz = Half, rHidden = Zero, RtRnc = Three, &
                 ThrCons = Zero, ThrMEP = Zero, ThrEne = Zero, ThrGrd = Zero, UserP = Zero, UserT(64) = Zero
logical(kind=iwp) :: Analytic_Hessian = .false., ApproxNADC = .false., Baker = .false., BSet = .false., CallLast = .true., &
                     Cubic = .false., Curvilinear = .true., ddV_Schlegel = .false., EDiffZero = .false., eMEPTest = .true., &
                     Fallback = .true., FindTS = .false., Force_dB = .false., HrmFrq_Show = .false., HSet = .false., &
                     HWRS = .true., Initiated = .false., isFalcon = .false., lCtoF = .false., lDoubleIso = .false., &
                     Line_Search = .true., lNmHss = .false., lOld = .false., lOld_Implicit = .false., lRP = .false., &
                     lSoft = .false., lTherm = .false., MEP = .false., MEPCons = .false., NADC = .false., Numerical = .false., &
                     PrQ = .false., Redundant = .false., Request_Alaska = .false., Request_RASSI = .false., rMEP = .false., &
                     SlStop = .false., Track = .false., TSConstraints = .false., TwoRunFiles = .false., User_Def = .false., &
                     WeightedConstraints = .false.
character(len=10) :: MEP_TYPE = 'SPHERE'
character(len=8) :: GrdLbl = '', StpLbl = ''
character(len=6) :: HUpMet = ' None ', UpMeth = '  RF  '
character(len=2) :: MEP_Algo = 'GS'
character :: Header(144) = ''
integer(kind=iwp), allocatable :: ANr(:), Atom(:), iBM(:), idBM(:), iCoset(:,:), jStab(:,:), mRowH(:), nqBM(:), nStab(:), NSup(:), &
                                  RootMap(:)
real(kind=wp), allocatable :: BM(:), BMx_kriging(:,:), Coor(:,:), Cx(:,:,:), dBM(:), Degen(:,:), DipM(:,:), dmass(:), dqInt(:,:), &
                              dqInt_Aux(:,:,:), Energy(:), Energy0(:), GNrm(:), Grd(:,:), Gx(:,:,:), Gx0(:,:,:), KtB(:,:), &
                              Lambda(:,:), MF(:,:), NAC(:,:,:), Q_nuclear(:), qInt(:,:), Shift(:,:), Weights(:)
real(kind=wp), allocatable, target :: Bmx(:,:), GradRef(:,:), R12(:,:), RefGeo(:,:)
logical(kind=iwp), allocatable :: Smmtrc(:,:)
character(len=LenIn), allocatable :: AtomLbl(:)
character(len=8), allocatable :: Lbl(:)

integer(kind=iwp), parameter :: Covalent_Bond = 0, &
                                vdW_Bond = 1, &
                                Fragments_Bond = 2, &
                                Magic_Bond = 3
character(len=*), parameter :: BondType(0:3) = ['Covalent_Bond ', &
                                                'vdW_Bond      ', &
                                                'Fragments_Bond', &
                                                'Magic_Bond    ']

public :: Analytic_Hessian, ANr, ApproxNADC, Atom, AtomLbl, Baker, Beta, Beta_Disp, BM, Bmx, BMx_kriging, BondType, BSet, &
          CallLast, CnstWght, Coor, Covalent_Bond, Cubic, Curvilinear, Cx, dBM, ddV_Schlegel, Degen, Delta, DipM, dmass, dMEPStep, &
          Dmp_Slapaf, dqInt, dqInt_Aux, E_Delta, EDiffZero, eMEPTest, Energy, Energy0, Fallback, FindTS, Force_dB, Fragments_Bond, &
          Free_Slapaf, Get_Slapaf, GNrm, GNrm_Threshold, GradRef, Grd, GrdLbl, GrdMax, Gx, Gx0, Header, HrmFrq_Show, HSet, HUpMet, &
          HWRS, iBM, iCoset, idBM, iInt, iNeg, iOptC, iOptH, IRC, iRef, iRow, iRow_c, isFalcon, iState, iter, jStab, KtB, Lambda, &
          Lbl, lCtoF, lDoubleIso, Line_Search, lNmHss, lOld, lOld_Implicit, lRP, lSoft, lTherm, Magic_Bond, Max_Center, MaxItr, &
          mB_Tot, mdB_Tot, MEP, MEP_Algo, MEP_TYPE, MEPCons, MEPNum, MF, Mode, mq, mRowH, mTROld, mTtAtm, MxItr, NAC, NADC, nBVec, &
          nDimBC, nFix, nLambda, nMEP, NmIter, nqBM, nsRot, nStab, NSup, Numerical, nUserPT, nWndw, PrQ, Q_nuclear, qInt, R12, &
          Redundant, RefGeo, Request_Alaska, Request_RASSI, rFuzz, rHidden, rMEP, RootMap, RtRnc, Shift, SlStop, Smmtrc, StpLbl, &
          StpMax, ThrCons, ThrEne, ThrGrd, ThrMEP, Track, TSConstraints, TwoRunFiles, UpMeth, User_Def, UserP, UserT, vdW_Bond, &
          WeightedConstraints, Weights

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

  integer(kind=iwp), intent(out) :: iter, mTROld
  integer(kind=iwp), intent(in) :: MaxItr, nsAtom, mLambda
  logical(kind=iwp), intent(out) :: lOld_Implicit
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
  if (iter > MaxItr) then
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

  logical(kind=iwp), intent(in) :: lStop, Just_Frequencies, lOld_Implicit
  real(kind=wp), intent(in) :: Energy_In
  integer(kind=iwp), intent(in) :: Iter, MaxItr, mTROld, nsAtom
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
    Information(4) = merge(1,0,lOld_Implicit)
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
