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

subroutine EdmistonRuedenberg(Functional,CMO,Thrs,ThrRot,ThrGrad,nBas,nOrb2Loc,nFro,nSym,nMxIter,Maximisation,Converged,Debug, &
                              Silent)
! Author: T.B. Pedersen
!
! Purpose: Edmiston-Ruedenberg localisation of occupied orbitals.

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Functional
real(kind=wp), intent(inout) :: CMO(*)
real(kind=wp), intent(in) :: Thrs, ThrRot, ThrGrad
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb2Loc(nSym), nFro(nSym), nMxIter
logical(kind=iwp), intent(in) :: Maximisation, Debug, Silent
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: irc, kOffC, nBasT, nFroT, nOrb2LocT
real(kind=wp) :: FracMem
character(len=80) :: Txt
character(len=*), parameter :: SecNam = 'EdmistonRuedenberg'

! Symmetry is NOT allowed.
! ------------------------

if (nSym /= 1) then
  call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
end if

! Initializations.
! ----------------

Functional = -huge(Functional)

nBasT = nBas(1)
nOrb2LocT = nOrb2Loc(1)
nFroT = nFro(1)

Converged = .false.

irc = -1
FracMem = 0.3_wp ! 30 percent of memory used as vector buffer
call Cho_X_Init(irc,FracMem)
if (irc /= 0) then
  write(Txt,'(A,I6)') 'Cho_X_Init returned',irc
  call SysAbendMsg(SecNam,'Cholesky initialization error:',Txt)
end if

! Localise orbitals.
! ------------------

kOffC = nBasT*nFroT+1
call EdmistonRuedenberg_Iter(Functional,CMO(kOffC),Thrs,ThrRot,ThrGrad,nBasT,nOrb2LocT,nMxIter,Maximisation,Converged,Debug,Silent)

! Finalizations.
! --------------

irc = -1
call Cho_X_Final(irc)
if (irc /= 0) then
  write(Txt,'(A,I6)') 'Cho_X_Final returned',irc
  call SysAbendMsg(SecNam,'Cholesky finalization error:',Txt)
end if

end subroutine EdmistonRuedenberg
