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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: irc, nSym, nBas(nSym), nOrb(nSym), nOcc(nSym), nVir(nSym)
real(kind=wp) :: CMO(*), PAO(*), Thr
character(len=*) :: Mode
#include "WrkSpc.fh"
integer(kind=iwp) :: DefLevel, ip_D, ip_Dum, ip_R, iSym, kOff1, kOffP, kOffR, l_D, l_Dum, l_R, Level, lMode, nOrthPs
real(kind=wp) :: ThrLoc, xNrm
logical(kind=iwp) :: Normalize, TestOrth
character(len=3) :: myMode
real(kind=wp), parameter :: DefThr = 1.0e-12_wp
character(len=*), parameter :: DefMode = 'ORT', SecNam = 'PAOLoc'

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

if (DefMode == 'RAW') then
  DefLevel = 1
else if (DefMode == 'CHO') then
  DefLevel = 2
else if (DefMode == 'ORT') then
  DefLevel = 3
else
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
  if (myMode == 'RAW') then
    Level = 1
  else if (myMode == 'CHO') then
    Level = 2
  else if (myMode == 'ORT') then
    Level = 3
  else
    Level = DefLevel
  end if
end if

! Make a dummy allocation to enable de-allocation by flushing.
! ------------------------------------------------------------

l_Dum = 1
call GetMem('PAOL_Dummy','Allo','Real',ip_Dum,l_Dum)

! Compute raw projected AOs.
! --------------------------

l_R = nBas(1)**2
do iSym=2,nSym
  l_R = l_R+nBas(iSym)**2
end do
call GetMem('PAOL_R','Allo','Real',ip_R,l_R)

Normalize = .true.
call GetRawPAOs(Work(ip_R),CMO,nBas,nOrb,nOcc,nVir,nSym,Normalize)

if (Level == 1) then
  call dCopy_(l_R,Work(ip_R),1,PAO,1)
  Go To 1 ! return after de-allocation
end if

! Use Cholesky decomposition to compute a linearly independent set
! of nonorthonormal PAOs.
! ----------------------------------------------------------------

l_D = nBas(1)**2
do iSym=2,nSym
  l_D = max(l_D,nBas(iSym)**2)
end do
call GetMem('PAOL_D','Allo','Real',ip_D,l_D)

if (Thr <= Zero) then
  ThrLoc = DefThr
else
  ThrLoc = Thr
end if

kOffR = ip_R
kOffP = 1
do iSym=1,nSym
  if (nVir(iSym) > 0) then
    call GetDens_Localisation(Work(ip_D),Work(kOffR),nBas(iSym),nBas(iSym))
    call ChoLoc(irc,Work(ip_D),PAO(kOffP),ThrLoc,xNrm,nBas(iSym),nVir(iSym))
    if (irc /= 0) Go To 1 ! return after de-allocation
  end if
  kOffR = kOffR+nBas(iSym)**2
  kOffP = kOffP+nBas(iSym)*nVir(iSym)
end do

if (Level == 2) then
  Go To 1 ! return after de-allocation
end if

! Orthonormalize the PAOs.
! ------------------------

kOffP = 1
kOffR = ip_R
do iSym=1,nSym
  kOff1 = kOffR+nBas(iSym)*nOcc(iSym)
  call dCopy_(nBas(iSym)*nVir(iSym),PAO(kOffP),1,Work(kOff1),1)
  kOffP = kOffP+nBas(iSym)*nVir(iSym)
  kOffR = kOffR+nBas(iSym)**2
end do
nOrthPs = 2 ! orthonormalization passes to ensure num. accuracy
call OrthoPAO_Localisation(Work(ip_R),nBas,nOcc,nVir,nSym,nOrthPs,TestOrth)
kOffP = 1
kOffR = ip_R
do iSym=1,nSym
  kOff1 = kOffR+nBas(iSym)*nOcc(iSym)
  call dCopy_(nBas(iSym)*nVir(iSym),Work(kOff1),1,PAO(kOffP),1)
  kOffP = kOffP+nBas(iSym)*nVir(iSym)
  kOffR = kOffR+nBas(iSym)**2
end do

if (Level == 3) then
  Go To 1 ! return after de-allocation
end if

! De-allocation by flushing.
! --------------------------

1 call GetMem('PAOL_Dummy','Flus','Real',ip_Dum,l_Dum)
call GetMem('PAOL_Dummy','Free','Real',ip_Dum,l_Dum)

end subroutine PAOLoc
