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

subroutine Chk_Input(irc)

! Author: T.B. Pedersen
! Purpose: Check input.
!
! Return codes:
!   irc = 0: all OK
!   irc < 0: all OK, but nothing to do
!   irc > 0: input error

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"

character*9 SecNam
parameter(SecNam='Chk_Input')

logical DoCholesky

irc = 0
doCholesky = .false.

nOrb2LocT = 0
do iSym=1,nSym
  n = nFro(iSym)+nOrb2Loc(iSym)
  if ((n < 0) .or. (n > nOrb(iSym))) then
    irc = irc+1
    write(6,*) SecNam,': nFro + nOrb2Loc out of bounds:'
    write(6,*) '    iSym     = ',iSym
    write(6,*) '    nFro     = ',nFro(iSym)
    write(6,*) '    nOrb2Loc = ',nOrb2Loc(iSym)
    write(6,*) '    nOrb     = ',nOrb(iSym)
  end if
  if (n > nBas(iSym)) then
    irc = irc+1
    write(6,*) SecNam,': nFro + nOrb2Loc > nBas:'
    write(6,*) '    iSym     = ',iSym
    write(6,*) '    nFro     = ',nFro(iSym)
    write(6,*) '    nOrb2Loc = ',nOrb2Loc(iSym)
    write(6,*) '    nBas     = ',nBas(iSym)
  end if
  nOrb2LocT = nOrb2LocT+nOrb2Loc(iSym)
end do
if (nOrb2LocT == 0) then
  irc = -1
  return
end if

if ((LocModel < 0) .or. (LocModel > nLocModel)) then
  write(6,*) SecNam,': LocModel must satisfy 0 <= LocModel <= ',nLocModel
  write(6,*) '    LocModel = ',LocModel
  irc = irc+1
end if

if (LocModel == 4) then
  call DecideOnCholesky(doCholesky)
  if (.not. doCholesky) then
    call SysAbendMsg(SecNam,'Edmiston-Ruedenberg localisation not possible:','Cholesky integrals required!')
  end if
end if

if (EvalER) then
  call DecideOnCholesky(doCholesky)
  if (.not. doCholesky) then
    write(6,*) SecNam,': evaluation of ER functional requires',' Cholesky decomposition of ERIs!'
    write(6,*) 'Evaluation of ER functional is cancelled...'
    EvalER = .false.
  end if
end if

if (Analysis .and. (.not. Test_Localisation)) then
  Test_Localisation = .true.
end if

end subroutine Chk_Input
