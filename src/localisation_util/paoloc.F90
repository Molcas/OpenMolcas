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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  PAOLoc
!
!> @brief
!>   Compute projected atomic orbitals (PAOs)
!> @author Thomas Bondo Pedersen
!>
!> @details
!> A set of linearly independent nonorthonormal
!> Projected Atomic Orbitals (PAOs) are computed by projecting
!> the AOs onto the virtual space. Subsequently, the linear
!> dependence may be removed through Cholesky decomposition of a
!> density-type matrix constructed from the linearly dependent
!> PAOs according to \f$ D_{a,b} = \sum_c \mathit{PAO}_{a,c} \mathit{PAO}_{b,c} \f$. Finally,
!> orthonormal linearly independent PAOs may be obtained by
!> multiplying with the inverse square root of the PAO overlap
!> matrix.
!>
!> The Mode input string controls which set is returned in array
!> PAO:
!>
!> - \p Mode = ``'RAW'``: return linearly dependent nonorthonormal PAOs. In this case, array PAO
!>                        should be dimensioned as \p nBas &timens \p nBas (in symmetry blocks).
!> - \p Mode = ``'CHO'``: return linearly independent nonorthonormal PAOs obtained by Cholesky
!>                        decomposition of the density-type matrix. In this case, array PAO should
!>                        be dimensioned as \p nBas &times; \p nVir (in symmetry blocks).
!> - \p Mode = ``'ORT'``: return linearly independent orthonormal PAOs. In this case, array PAO
!>                        should be dimensioned as \p nBas &times; \p nVir (in symmetry blocks).
!>
!> If Mode is anything else, \p Mode = ``'ORT'`` is assumed.
!> Note that if \p nOcc+nVir < \p nBas the remaining \p nBas-nOcc-nVir
!> orbitals are considered part of the orthogonal complement of
!> the virtual space (and are projected out of the AOs when
!> obtaining the raw PAOs). Thus, one may obtain PAOs for any
!> subset of orbitals by specifying \p nOcc and \p nVir arrays
!> appropriately.
!>
!> @note
!> Needs the AO overlap matrix on disk.
!>
!> @param[out] irc  return code
!> @param[in]  CMO  MO coefficients
!> @param[out] PAO  PAOs
!> @param[in]  Thr  Cholesky decomposition threshold
!> @param[in]  nBas Number of basis functions/irrep
!> @param[in]  nOrb Number of orbitals/irrep
!> @param[in]  nOcc Number of occupied orbs/irrep
!> @param[in]  nVir Number of virtual orbs/irrep
!> @param[in]  nSym Number of irreps
!> @param[in]  Mode Mode of calculation
!***********************************************************************

subroutine PAOLoc(irc,CMO,PAO,Thr,nBas,nOrb,nOcc,nVir,nSym,Mode)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: CMO(*), Thr
real(kind=wp), intent(_OUT_) :: PAO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), nOcc(nSym), nVir(nSym)
character(len=*), intent(in) :: Mode
integer(kind=iwp) :: DefLevel, i, iSym, kOff, l_D, Level, lMode, nOrthPs
real(kind=wp) :: ThrLoc, xNrm
logical(kind=iwp) :: Normalize, TestOrth
character(len=3) :: myMode
type(DSBA_Type) :: P, R
real(kind=wp), allocatable :: D(:)
real(kind=wp), parameter :: DefThr = 1.0e-12_wp
character(len=*), parameter :: DefMode = 'ORT', LevMode(3) = ['RAW','CHO','ORT'], SecNam = 'PAOLoc'

#if defined (_DEBUGPRINT_)
TestOrth = .true.
#else
TestOrth = .false.
#endif

! Set return code.
! ----------------

irc = 0

! Set DefLevel from DefMode.
! --------------------------

DefLevel = 0
do i=1,size(LevMode)
  if (DefMode == LevMode(i)) then
    DefLevel = i
    exit
  end if
end do
if (DefLevel == 0) then
  call SysAbendMsg(SecNam,'DefMode not recognized!','Note: this is a programming error...')
  DefLevel = -999999
end if

! Interpret Mode input.
! ---------------------

lMode = len(Mode)
if (lMode < 3) then
  Level = DefLevel
else
  myMode = Mode(1:3)
  call UpCase(myMode)
  Level = DefLevel
  do i=1,size(LevMode)
    if (myMode == LevMode(i)) then
      Level = i
      exit
    end if
  end do
end if

! Compute raw projected AOs.
! --------------------------

call Allocate_DT(R,nBas,nBas,nSym,label='R')

Normalize = .true.
call GetRawPAOs(R%A0,CMO,nBas,nOrb,nOcc,nVir,nSym,Normalize)

if (Level == 1) then
  PAO(1:size(R%A0)) = R%A0
  call FreeMem()
  return
end if

call Allocate_DT(P,nBas,nVir,nSym,label='P',Ref=PAO)

! Use Cholesky decomposition to compute a linearly independent set
! of nonorthonormal PAOs.
! ----------------------------------------------------------------

l_D = nBas(1)**2
do iSym=2,nSym
  l_D = max(l_D,nBas(iSym)**2)
end do
call mma_allocate(D,l_D,label='PAOL_D')

if (Thr <= Zero) then
  ThrLoc = DefThr
else
  ThrLoc = Thr
end if

do iSym=1,nSym
  if (nVir(iSym) > 0) then
    call GetDens_Localisation(D,R%SB(iSym)%A2,nBas(iSym),nBas(iSym))
    call ChoLoc(irc,D,P%SB(iSym)%A2,ThrLoc,xNrm,nBas(iSym),nVir(iSym))
    if (irc /= 0) then
      call FreeMem()
      return
    end if
  end if
end do

if (Level == 2) then
  call FreeMem()
  return
end if

! Orthonormalize the PAOs.
! ------------------------

do iSym=1,nSym
  kOff = nOcc(iSym)+1
  R%SB(iSym)%A2(:,kOff:kOff+nVir(iSym)-1) = P%SB(iSym)%A2(:,1:nVir(iSym))
end do
nOrthPs = 2 ! orthonormalization passes to ensure num. accuracy
call OrthoPAO_Localisation(R%A0,nBas,nOcc,nVir,nSym,nOrthPs,TestOrth)
do iSym=1,nSym
  kOff = nOcc(iSym)+1
  P%SB(iSym)%A2(:,1:nVir(iSym)) = R%SB(iSym)%A2(:,kOff:kOff+nVir(iSym)-1)
end do

call FreeMem()

return

contains

subroutine FreeMem()

  call Deallocate_DT(R)
  call Deallocate_DT(P)
  if (allocated(D)) call mma_deallocate(D)

end subroutine FreeMem

end subroutine PAOLoc
