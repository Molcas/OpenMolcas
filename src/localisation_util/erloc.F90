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
!  ERLoc
!
!> @brief
!>   Edmiston--Ruedenberg localization of occupied orbitals
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Performs iterative (eta sweeps) Edmiston--Ruedenberg
!> localization of
!> molecular orbitals. CMO is assumed to be the full set of
!> orbitals, but only the first \p nOcc will be localized. On exit,
!> the localized MOs are returned in CMO.
!> Note that symmetry is not allowed.
!>
!> If successful, \p irc = ``0`` is returned. If \p irc = ``-1``
!> is returned, an error was found in the input and nothing has
!> been done. If \p irc = ``1`` is returned, the iterative procedure did
!> not converge to within the threshold (\p Thr, \p ThrGrad) in the
!> requested max. number of iterations (\p MxIter). If the user
!> specifies negative thresholds, \p Thr = ``1.0e-6_wp`` and \p ThrGrad = ``1.0e-3_wp``
!> will be used.
!> If the user specifies a negative screening threshold,
!> \p ThrRot = ``1.0e-10_wp`` is used.
!>
!> @param[out]    irc     Return code
!> @param[in,out] CMO     Molecular orbital coefficients
!> @param[in]     Thr     Threshold for functional
!> @param[in]     ThrGrad Threshold for gradient
!> @param[in]     ThrRot  Screening threshold
!> @param[in]     MxIter  Max. number of iterations
!> @param[in]     nBas    Number of basis functions per irrep
!> @param[in]     nOcc    Number of occupied orbitals to localize
!> @param[in]     nFro    Number of frozen occupied orbitals
!> @param[in]     nSym    Number of irreps
!> @param[in]     Silent  Flag to avoid printing
!***********************************************************************

subroutine ERLoc(irc,CMO,Thr,ThrGrad,ThrRot,MxIter,nBas,nOcc,nFro,nSym,Silent)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(inout) :: CMO(*)
real(kind=wp), intent(in) :: Thr, ThrGrad, ThrRot
integer(kind=iwp), intent(in) :: MxIter, nSym, nBas(nSym), nOcc(nSym), nFro(nSym)
logical(kind=iwp), intent(in) :: Silent
#include "Molcas.fh"
integer(kind=iwp) :: iSym, nBasT, nOccT
real(kind=wp) :: Functional, ThrGLoc, ThrLoc, ThrRotLoc
logical(kind=iwp) :: Converged, Debug, Maximization
character(len=*), parameter :: SecNam = 'ERLoc'

! Initialization.
! ---------------

irc = 0
if (MxIter < 1) return

nBasT = nBas(1)
do iSym=2,nSym
  nBasT = nBasT+nBas(iSym)
end do
if (nBasT < 1) return

nOccT = nOcc(1)
do iSym=2,nSym
  nOccT = nOccT+nOcc(iSym)
end do
if (nOccT < 1) return

! Edmiston-Ruedenberg does not work with symmetry.
! TODO/FIXME: we might consider de-symmetrization.
! ------------------------------------------------

if (nSym /= 1) then
  irc = -1
  return
end if

! Localize.
! ---------

Functional = -huge(Functional)
if (Thr <= Zero) then
  ThrLoc = 1.0e-6_wp
else
  ThrLoc = Thr
end if
if (ThrGrad <= Zero) then
  ThrGLoc = 1.0e-3_wp
else
  ThrGLoc = ThrGrad
end if
if (ThrRot < Zero) then
  ThrRotLoc = 1.0e-10_wp
else
  ThrRotLoc = ThrRot
end if
Maximization = .true.
Converged = .false.
Debug = .false.
call EdmistonRuedenberg(Functional,CMO,ThrLoc,ThrRotLoc,ThrGLoc,nBas,nOcc,nFro,nSym,MxIter,Maximization,Converged,Debug,Silent)

! Check convergence.
! ------------------

if (.not. Converged) irc = 1

end subroutine ERLoc
