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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine validate_wfn()

use rasscf_global, only: iRoot, lRoots, NAC, NIN, nRoots
use general_data, only: ispin, nactel, nash, nbas, nbas, ndel, nelec3, nfro, nhole1, nish, norb, nrs1, nrs1t, nrs2, nrs2t, nrs3, &
                        nrs3t, nssh, nsym, ntot, stsym
use Molcas, only: MxAct, MxBas, MxIna, MxOrb, MxRoot
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ia0, ib0, ic0
logical(kind=iwp) :: is_negative, is_too_large

if ((NSYM /= 1) .and. (NSYM /= 2) .and. (NSYM /= 4) .and. (NSYM /= 8)) then
  call WarningMessage(2,'Wrong nr of symmetries.')
  write(u6,*) '*************** ERROR ***************'
  write(u6,'(1X,A,I8)') 'Nr of symmetries NSYM=',NSYM
  write(u6,*) ' Only possible values are 1,2,4 or 8.'
  write(u6,*) '*************************************'
  call Quit_OnUserError()
end if
if (STSYM > NSYM) then
  call WarningMessage(2,'Wrong symmetry.')
  write(u6,*) '***************** ERROR *****************'
  write(u6,'(1X,A,I8)') 'CHKINP Error: Wrong symmetry.'
  write(u6,'(1X,A,I8)') 'State symmetry   STSYM=',STSYM
  write(u6,'(1X,A,I8)') 'Point group order NSYM=',NSYM
  write(u6,*) '************************************************'
  call Quit_OnUserError()
end if
if (NTOT > mxOrb) then
  call WarningMessage(2,'Too many orbitals.')
  write(u6,*) '************ ERROR ******************'
  write(u6,'(1X,A,I8)') 'Too many orbitals NTOT=',NTOT
  write(u6,'(1X,A,I8)') 'Limit is MXORB=',MXORB
  write(u6,*) '*************************************'
  call Quit_OnUserError()
end if
if (NAC > mxAct) then
  call WarningMessage(2,'Too many active orbitals.')
  write(u6,*) '*************** ERROR ***************'
  write(u6,'(1X,A,I8)') 'Too many active orbitals NAC=',NAC
  write(u6,'(1X,A,I8)') 'Limit is MXACT=',MXACT
  write(u6,*) '*************************************'
  call Quit_OnUserError()
end if
if (NIN > mxIna) then
  call WarningMessage(2,'Too many inactive orbitals.')
  write(u6,*) '*************** ERROR ***************'
  write(u6,'(1X,A,I8)') 'Too many inactive orbitals NIN=',NIN
  write(u6,'(1X,A,I8)') 'Limit is MXINA=',MXINA
  write(u6,*) '*************************************'
  call Quit_OnUserError()
end if
if (NACTEL > 2*NAC) then
  call WarningMessage(2,'Too many active electrons.')
  write(u6,*) '********************* ERROR **********************'
  write(u6,'(1X,A,I6)') 'Too many active electrons NACTEL=',NACTEL
  write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of active orbitals=',2*NAC
  write(u6,*) '**************************************************'
  call Quit_OnUserError()
end if
if (NHOLE1 > 2*NRS1T) then
  call WarningMessage(2,'Too many holes in Ras1.')
  write(u6,*) '******************** WARNING *********************'
  write(u6,'(1X,A,I6)') 'You allow too many holes in Ras1 NHOLE1=',NHOLE1
  write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras1 orbitals=',2*NRS1T
  write(u6,*) '**************************************************'
  call Quit_OnUserError()
end if
if (NELEC3 > 2*NRS3T) then
  call WarningMessage(1,'Too many electrons in Ras3.')
  write(u6,*) '******************** WARNING *********************'
  write(u6,'(1X,A,I6)') 'You allow too many electrons in Ras3 NELEC3=',NELEC3
  write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras3 orbitals=',2*NRS3T
  write(u6,*) '**************************************************'
  call Quit_OnUserError()
end if
if (NACTEL > 2*(NRS1T+NRS2T)+NELEC3) then
  call WarningMessage(2,'Too many active electrons.')
  write(u6,*) '********************* ERROR **********************'
  write(u6,'(1X,A,I8)') 'Too many active electrons NACTEL=',NACTEL
  write(u6,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
  write(u6,*) '**************************************************'
  call Quit_OnUserError()
end if
if (NACTEL < 2*NRS1T-NHOLE1) then
  call WarningMessage(2,'Too few active electrons.')
  write(u6,*) '********************* ERROR **********************'
  write(u6,'(1X,A,I8)') 'Too few active electrons NACTEL=',NACTEL
  write(u6,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
  write(u6,*) '**************************************************'
  call Quit_OnUserError()
end if

is_negative = .false.
if (any(nBas(1:nSym) < 0)) is_negative = .true.
if (any(nFro(1:nSym) < 0)) is_negative = .true.
if (any(nDel(1:nSym) < 0)) is_negative = .true.
if (any(nOrb(1:nSym) < 0)) is_negative = .true.
if (any(nISh(1:nSym) < 0)) is_negative = .true.
if (any(nASh(1:nSym) < 0)) is_negative = .true.
if (any(nSSh(1:nSym) < 0)) is_negative = .true.
if (any(nRS1(1:nSym) < 0)) is_negative = .true.
if (any(nRS2(1:nSym) < 0)) is_negative = .true.
if (any(nRS3(1:nSym) < 0)) is_negative = .true.
is_too_large = .false.
if (any(nBas(1:nSym) > mxBas)) is_too_large = .true.
if (any(nFro(1:nSym) > mxBas)) is_too_large = .true.
if (any(nDel(1:nSym) > mxBas)) is_too_large = .true.
if (any(nOrb(1:nSym) > mxBas)) is_too_large = .true.
if (any(nISh(1:nSym) > mxBas)) is_too_large = .true.
if (any(nASh(1:nSym) > mxBas)) is_too_large = .true.
if (any(nSSh(1:nSym) > mxBas)) is_too_large = .true.
if (any(nRS1(1:nSym) > mxBas)) is_too_large = .true.
if (any(nRS2(1:nSym) > mxBas)) is_too_large = .true.
if (any(nRS3(1:nSym) > mxBas)) is_too_large = .true.
if (is_negative .or. is_too_large) then
  write(u6,*)
  write(u6,*) '****************** ERROR *******************'
  call WarningMessage(2,'Erroneous nr of orbitals.')
  write(u6,*) 'Inappropriate nr of orbitals. One or more of'
  write(u6,*) 'these orbital counts is wrong or too large.'
  if (is_negative) write(u6,*) ' Negative values.'
  if (is_too_large) write(u6,*) ' Extremely large values.'
  write(u6,'(1X,A,8I4)') '   All orbitals:',NBAS(1:NSYM)
  write(u6,'(1X,A,8I4)') '         Frozen:',NFRO(1:NSYM)
  write(u6,'(1X,A,8I4)') '       Inactive:',NISH(1:NSYM)
  write(u6,'(1X,A,8I4)') '         Active:',NASH(1:NSYM)
  write(u6,'(1X,A,8I4)') '          RAS-1:',NRS1(1:NSYM)
  write(u6,'(1X,A,8I4)') '          RAS-2:',NRS2(1:NSYM)
  write(u6,'(1X,A,8I4)') '          RAS-3:',NRS3(1:NSYM)
  write(u6,'(1X,A,8I4)') '      Secondary:',NSSH(1:NSYM)
  write(u6,'(1X,A,8I4)') '        Deleted:',NDEL(1:NSYM)
  write(u6,'(1X,A,8I4)') 'Basis functions:',NBAS(1:NSYM)
  if (is_negative) write(u6,*) ' Have you used a too small basis set?'
  write(u6,*) '********************************************'
  call Quit_OnUserError()
end if

IB0 = ISPIN-1
IA0 = (NACTEL-IB0)/2
IC0 = (NAC-IA0-IB0)
if ((2*IA0+IB0 /= NACTEL) .or. (IA0 < 0) .or. (IB0 < 0) .or. (IC0 < 0)) then
  call WarningMessage(2,'No such wave function.')
  write(u6,*) '******************** ERROR *********************'
  write(u6,*) 'The following combined specifications are wrong.'
  write(u6,'(1X,A,8I4)') 'Nr of active electrons NACTEL=',NACTEL
  write(u6,'(1X,A,8I4)') 'Nr of active orbitals     NAC=',NAC
  write(u6,'(1X,A,8I4)') 'Spin degeneracy         ISPIN=',ISPIN
  write(u6,*) 'There can be no such wave function.'
  write(u6,*) '************************************************'
  call Quit_OnUserError()
end if

if ((NROOTS > mxRoot) .or. (LROOTS > mxROOT)) then
  call WarningMessage(2,'Max roots exceeded.')
  write(u6,*) '***************** ERROR *****************'
  write(u6,'(1X,A,I6)') 'Input Error: Max roots exceeded.',mxRoot
  write(u6,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
  write(u6,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
  write(u6,*) '************************************************'
  call Quit_OnUserError()
end if
if ((minval(IROOT(1:NROOTS)) < 0) .or. (maxval(IROOT(1:NROOTS)) > LROOTS)) then
  call WarningMessage(2,'Root specifications are wrong.')
  write(u6,*) '***************** ERROR *****************'
  write(u6,*) 'CHKINP Error: Root specifications are wrong.'
  write(u6,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
  write(u6,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
  write(u6,*) '************************************************'
  call Quit_OnUserError()
end if

if (((NHOLE1 < 0) .and. (NRS1T /= 0)) .or. ((NELEC3 < 0) .and. (NRS3T /= 0)) .or. ((NACTEL < 0) .and. (NRS2T /= 0))) then
  call WarningMessage(2,'Error in RAS specification.')
  write(u6,*) '***************** ERROR *****************'
  write(u6,*) 'Error in RAS specification.'
  write(u6,'(1X,A,I8)') 'Max holes in Ras1,     NHOLE1=',NHOLE1
  write(u6,'(1X,A,I8)') 'Max electrons in Ras3, NELEC3=',NELEC3
  write(u6,'(1X,A,I8)') 'Nr of active electrons NACTEL=',NACTEL
  write(u6,*) '*****************************************'
  call Quit_OnUserError()
end if

end subroutine validate_wfn
