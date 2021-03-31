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

implicit real*8(a-h,o-z)
character*4 Model
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"
#include "debug.fh"

character*18 SecNam
parameter(SecNam='Localise_Iterative')

character*4 myModel
character*80 Txt
logical Converged

irc = 0
Functional = -9.9d9
Converged = .false.

! Generate Cholesky start guess, if requested.
! --------------------------------------------

if (ChoStart) then
  Thrs_Save = Thrs
  Thrs = 1.0d-12
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
  write(6,'(//,1X,A)') 'Pipek-Mezey localisation'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  write(6,'(1X,A,1X,D12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call PipekMezey(Functional,Work(ipCMO),Thrs,ThrRot,ThrGrad,Name,nBas,nOrb2Loc,nFro,nSym,nAtoms,nMxIter,Maximisation,Converged, &
                  Debug,Silent)
else if (myModel == 'BOYS') then
  !if (.not. Silent) then
  write(6,'(/,1X,A)') 'Boys localisation'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  write(6,'(1X,A,1X,D12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call Boys(Functional,Work(ipCMO),Thrs,ThrRot,ThrGrad,nBas,nOrb2Loc,nFro,nSym,nMxIter,Maximisation,Converged,Debug,Silent)
else if (myModel == 'EDMI') then
  !if (.not. Silent) then
  write(6,'(/,1X,A)') 'Edmiston-Ruedenberg localisation'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  !write(6,'(1X,A,1X,D12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call EdmistonRuedenberg(Functional,Work(ipCMO),Thrs,ThrRot,ThrGrad,nBas,nOrb2Loc,nFro,nSym,nMxIter,Maximisation,Converged,Debug, &
                          Silent)
else
  write(Txt,'(A,A4)') 'Model = ',Model
  call SysAbendMsg(SecNam,'Unknown model',Txt)
end if

if (.not. Converged) then
  irc = 1
  return
end if

end subroutine Localise_Iterative
