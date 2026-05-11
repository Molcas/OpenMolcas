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
!#define _DEBUGPRINT_
!#define _SCR_

subroutine Localise_Iterative(irc,Model,Functional)
! Author: T.B. Pedersen
!
! Purpose: Iterative localisation of orbitals.
!          Models implemented:
!            Pipek-Mezey         [MODEL='PIPE']
!            Boys                [MODEL='BOYS']
!            Edmiston-Ruedenberg [MODEL='EDMI']

use Localisation_globals, only: ChoStart, CMO, nBas, nFro, nOrb2Loc, nSym, ThrGrad, ThrRot, Thrs, ScrFac,OptMeth, ChargeType
use Definitions, only: wp, iwp, u6
use Constants, only: Zero

implicit none
integer(kind=iwp), intent(out) :: irc
character(len=4), intent(in) :: Model
real(kind=wp), intent(out) :: Functional
integer(kind=iwp) :: iSym
real(kind=wp) :: Thrs_Save, xNrm
character(len=80) :: Txt
character(len=4) :: myModel
logical(kind=iwp) :: Converged
character(len=*), parameter :: SecNam = 'Localise_Iterative'

#ifdef _SCR_
if (ScrFac == Zero .and. OptMeth == 4 .or. OptMeth == 5 .or. OptMeth == 2) ScrFac = 0.5
#endif

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
#ifdef _DEBUGPRINT_
call recprt("CMO before localisation"," ",CMO,nBas,nOrb2Loc)
#endif

myModel = Model
call UpCase(myModel)
if (myModel == 'PIPE') then
  !if (.not. Silent) then
  write(u6,'(//,1X,A)') 'Pipek-Mezey localisation'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  If (ScrFac/=Zero) write(u6,'(1X,A,1X,ES12.4)')'Scrambling factor    :',ScrFac
  write(u6,'(1X,A,1X,ES12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  If (OptMeth == 1) then
    write(u6,'(1X,A)') 'Optimization Method  : Jacobi Sweeps'
  else if (OptMeth == 2) then
    write(u6,'(1X,A)') 'Optimization Method  : Newton Raphson'
  else if (OptMeth == 3) then
    write(u6,'(1X,A)') 'Optimization Method  : Gradient Ascent'
  else if (OptMeth == 4) then
    write(u6,'(1X,A)') 'Optimization Method  : GEK (fullspace)'
  else if (OptMeth == 5) then
    write(u6,'(1X,A)') 'Optimization Method  : S-GEK'
  else if (OptMeth == 6) then
    write(u6,'(1X,A)') 'Optimization Method  : S-GEK with Jacobi Sweep start'
  end if
  If (ChargeType == 1) then
    write(u6,'(1X,A)') 'Framework for PMLoc  : Mulliken charges'
  else if (OptMeth == 2) then
    write(u6,'(1X,A)') 'Framework for PMLoc  : Loewdin charges'
  end if
  !end if
  call PipekMezey        (Functional,CMO,nBas,nOrb2Loc,nFro,nSym,Converged)
else if (myModel == 'BOYS') then
  !if (.not. Silent) then
  write(u6,'(/,1X,A)') 'Boys localisation'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  If (ScrFac/=Zero) write(u6,'(1X,A,1X,ES12.4)') 'Scrambling factor     ',ScrFac
  write(u6,'(1X,A,1X,ES12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call Boys              (Functional,CMO,nBas,nOrb2Loc,nFro,nSym,Converged)
else if (myModel == 'EDMI') then
  !if (.not. Silent) then
  write(u6,'(/,1X,A)') 'Edmiston-Ruedenberg localisation'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',Thrs,' (functional)'
  write(u6,'(1X,A,1X,ES12.4,A)') 'Convergence threshold:',ThrGrad,' (gradient)'
  If (ScrFac/=Zero) write(u6,'(1X,A,1X,ES12.4)') 'Scrambling factor     ',ScrFac
  !write(u6,'(1X,A,1X,ES12.4,A)') 'Screening threshold  :',ThrRot,' (orbital rotations)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  call EdmistonRuedenberg(Functional,CMO,nBas,nOrb2Loc,nFro,nSym,Converged)
else
  write(Txt,'(A,A4)') 'Model = ',Model
  call SysAbendMsg(SecNam,'Unknown model',Txt)
end if

#ifdef _DEBUGPRINT_
call recprt("CMO after localisation"," ",CMO,nBas,nOrb2Loc)
#endif


if (.not. Converged) then
  irc = 1
  return
end if

end subroutine Localise_Iterative
