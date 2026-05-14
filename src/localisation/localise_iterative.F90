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

subroutine Localise_Iterative(irc,Model,Functional)
! Author: T.B. Pedersen
!
! Purpose: Iterative localisation of orbitals.
!          Models implemented:
!            Pipek-Mezey         [MODEL='PIPE']
!            Boys                [MODEL='BOYS']
!            Edmiston-Ruedenberg [MODEL='EDMI']

use Localisation_globals, only: ChoStart, CMO, Maximisation, BName, nAtoms, nMxIter, nBas, nFro, nOrb2Loc, nSym, Silent, ThrGrad, &
                                ThrRot, Thrs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
character(len=4), intent(in) :: Model
real(kind=wp), intent(out) :: Functional
integer(kind=iwp) :: iSym
real(kind=wp) :: Thrs_Save, xNrm
character(len=80) :: Txt
character(len=4) :: myModel
logical(kind=iwp) :: Converged
logical(kind=iwp), parameter :: debug = .false.
character(len=*), parameter :: SecNam = 'Localise_Iterative'

irc = 0
Functional = -huge(Functional)
Converged = .false.

! Generate Cholesky start guess, if requested.
! --------------------------------------------

if (ChoStart) then
  Thrs_Save = Thrs
  Thrs = 1.0e-12_wp
  call Localise_Noniterative(irc,'Chol',xNrm)
  if (irc /= 0) then
    write(Txt,'(A,I4)') 'Return code:',irc
    call SysAbendMsg(SecNam,'Localise_Noniterative failed!',Txt)
  end if
  Thrs = Thrs_Save
end if

! Localise.
! ---------

myModel = Model
call UpCase(myModel)
if (myModel == 'PIPE') then
  !if (.not. Silent) then
  write(u6,'(//,1X,A)') 'Pipek-Mezey localisation'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call PipekMezey(Functional,CMO,Thrs,ThrRot,ThrGrad,BName,nBas,nOrb2Loc,nFro,nSym,nAtoms,nMxIter,Maximisation,Converged, &
                  Debug,Silent)
else if (myModel == 'BOYS') then
  !if (.not. Silent) then
  write(u6,'(/,1X,A)') 'Boys localisation'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call Boys(Functional,CMO,Thrs,ThrRot,ThrGrad,nBas,nOrb2Loc,nFro,nSym,nMxIter,Maximisation,Converged,Debug,Silent)
else if (myModel == 'EDMI') then
  !if (.not. Silent) then
  write(u6,'(/,1X,A)') 'Edmiston-Ruedenberg localisation'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  !write(u6,'(1X,A,1X,ES12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call EdmistonRuedenberg(Functional,CMO,Thrs,ThrRot,ThrGrad,nBas,nOrb2Loc,nFro,nSym,nMxIter,Maximisation,Converged,Debug,Silent)
else
  write(Txt,'(A,A4)') 'Model = ',Model
  call SysAbendMsg(SecNam,'Unknown model',Txt)
end if

if (.not. Converged) then
  irc = 1
  return
end if

end subroutine Localise_Iterative
