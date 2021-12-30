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
!  PMLoc
!
!> @brief
!>   Pipek--Mezey localization of occupied orbitals
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Performs iterative (Jacobi sweeps) Pipek--Mezey localization of
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
!> specifies negative thresholds, \p Thr = ``1.0d-6`` and \p ThrGrad = ``1.0d-3``
!> will be used.
!> If the user specifies a negative screening threshold,
!> \p ThrRot = ``1.0d-10`` is used.
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

subroutine PMLoc(irc,CMO,Thr,ThrGrad,ThrRot,MxIter,nBas,nOcc,nFro,nSym,Silent)

implicit real*8(a-h,o-z)
real*8 CMO(*)
integer nBas(nSym), nOcc(nSym), nFro(nSym)
logical Silent
#include "Molcas.fh"
character*5 SecNam
parameter(SecNam='PMLoc')
character*(LENIN8) myName(MxAtom)
character*80 Txt
logical Maximization, Converged, Debug

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

! Pipek-Mezey does not work with symmetry.
! TODO/FIXME: we might consider de-symmetrization.
! ------------------------------------------------

if (nSym /= 1) then
  irc = -1
  return
end if

! Read number of atoms, atomic labels, and basis functions labels
! from runfile.
! ---------------------------------------------------------------

call Get_nAtoms_All(nAtoms)
if ((nAtoms < 1) .or. (nAtoms > MxAtom)) then
  write(Txt,'(A,I9)') 'nAtoms =',nAtoms
  call SysAbendMsg(SecNam,'Atom limit exceeded!',Txt)
end if
call Get_cArray('Unique Basis Names',myName,LENIN8*nBasT)

! Localize.
! ---------

Functional = -9.9d9
if (Thr <= 0.0d0) then
  ThrLoc = 1.0d-6
else
  ThrLoc = Thr
end if
if (ThrGrad <= 0.0d0) then
  ThrGLoc = 1.0d-3
else
  ThrGLoc = ThrGrad
end if
if (ThrRot < 0.0d0) then
  ThrRotLoc = 1.0d-10
else
  ThrRotLoc = ThrRot
end if
Maximization = .true.
Converged = .false.
Debug = .false.
call PipekMezey(Functional,CMO,ThrLoc,ThrRotLoc,ThrGLoc,myName,nBas,nOcc,nFro,nSym,nAtoms,MxIter,Maximization,Converged,Debug, &
                Silent)

! Check convergence.
! ------------------

if (.not. Converged) irc = 1

end subroutine PMLoc
