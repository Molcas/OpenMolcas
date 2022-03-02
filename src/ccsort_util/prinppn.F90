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
! Copyright (C) 1994, Markus P. Fuelscher                              *
!               1994, Per Ake Malmqvist                                *
!               Pavel Neogrady                                         *
!***********************************************************************

subroutine PrInpPN()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     - echo the input parameters                                      *
!                                                                      *
!     calling parameters: none                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P.-AA. Malmqvist                              *
!     University of Lund, Sweden, 1994                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     Modified by P.N.                                                 *
!                                                                      *
!***********************************************************************

use ccsort_global, only: cckey, clopkey, fullprint, iokey, ISCF, ISPIN, LROOT, LSYM, NACTEL, NASH, NASHT, NBAS, NCONF, NDEL, &
                         ndelr, NELE3, NFRO, nfror, NHOLE1, NISH, NISHT, noop, NROOTS, NSSH, NSSHT, NSYM, t3key, zrkey
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iSym
!LD character(len=120) :: BlLine, Line, StLine
!LD character(len=8) :: Fmt1, Fmt2

!----------------------------------------------------------------------*
!     Start and define the paper width                                 *
!----------------------------------------------------------------------*
!lPaper = 132
!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*
!lLine = Len(Line)
!do i=1,lLine
!  BlLine(i:i) = ' '
!  StLine(i:i) = '*'
!end do
!lPaper = 132
!left = (lPaper-lLine)/2
!write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
!write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
!     Print the project title                                          *
!----------------------------------------------------------------------*
!if (nTit > 0) then
!  write(u6,*)
!  nLine = nTit+5
!  do i=1,nLine
!    Line = BlLine
!    If ((i == 1) .or. (i == nLine)) Line = StLine
!    If (i == 3) Line = 'Project:'
!    If ((i >= 4) .and. (i <= nLine-2)) write(Line,'(18A4)') (Title(i-3,j),j=1,18)
!    call Center_Text(Line)
!    write(u6,Fmt1) '*'//Line//'*'
!  end do
!  write(u6,*)
!end if
!----------------------------------------------------------------------*
!     Stop if NOOPeration key is used                                  *
!----------------------------------------------------------------------*
if (noop == 1) then
  write(u6,'(6X,A)') ' No operation is required'
  write(u6,'(6X,A)') ' Happy Landing '
  call Finish(0)
end if
!----------------------------------------------------------------------*
!     Print iokey                                                      *
!----------------------------------------------------------------------*
if (iokey == 1) write(u6,'(6X,A)') 'Standard Fortran IO handling used '

if (iokey == 2) write(u6,'(6X,A)') 'MOLCAS DA IO handling used '
!----------------------------------------------------------------------*
!     Print zrkey                                                      *
!----------------------------------------------------------------------*
if (fullprint == 2) then
  if (zrkey == 1) write(u6,'(6X,A)') 'Separate V and Ind IO'

  if (zrkey == 0) write(u6,'(6X,A)') 'Simultanneous V and Ind IO'
end if
!----------------------------------------------------------------------*
!     Print cckey and t3key                                            *
!----------------------------------------------------------------------*
if (cckey == 1) write(u6,'(6X,A)') 'Integrals for CCSD will be produced'

if (t3key == 1) write(u6,'(6X,A)') 'Integrals for Noniterative T3 will be produced'

if (clopkey == 1) then
  write(u6,'(6X,A)') 'ROHF open shell reference function'
else
  write(u6,'(6X,A)') 'RHF closed shell reference function'
end if
!----------------------------------------------------------------------*
!     Print allocation and printing parameters                         *
!----------------------------------------------------------------------*
!if (maxspace == 0) then
!  write(u6,'(6X,A)') ' Allocatable work space   : Unlimited'
!else
!  write(u6,'(6X,A,I10)') ' Allocatable work space   : ',maxspace
!end if
!if (fullprint == 0) then
!  write(u6,'(6X,A)') ' Level of output printing : Minimal'
!else if (fullprint == 1) then
!  write(u6,'(6X,A)') ' Level of output printing : Medium '
!else if (fullprint == 2) then
!  write(u6,'(6X,A)') ' Level of output printing : Full'
!end if
!----------------------------------------------------------------------*
!     Print actual frozen and deleted orbitals                         *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'Actual numbers of frozen and deleted orbitals :'
write(u6,'(6X,A)') '-----------------------------------------------'
write(u6,*)
write(u6,'(6X,A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Frozen orbitals',(nFror(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Deleted orbitals',(nDelr(iSym),iSym=1,nSym)
write(u6,*)
!----------------------------------------------------------------------*
!     Print orbital and wavefunction specifications                    *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'Wave function specifications from previous RASSCF:'
write(u6,'(6X,A)') '--------------------------------------------------'
write(u6,*)
write(u6,'(6X,A,T45,I6)') 'Number of closed shell electrons',2*NISHT
write(u6,'(6X,A,T45,I6)') 'Number of electrons in active shells',NACTEL
write(u6,'(6X,A,T45,I6)') 'Max number of holes in RAS1 space',NHOLE1
write(u6,'(6X,A,T45,I6)') 'Max number of electrons in RAS3 space',NELE3
write(u6,'(6X,A,T45,I6)') 'Number of inactive orbitals',NISHT
write(u6,'(6X,A,T45,I6)') 'Number of active orbitals',NASHT
write(u6,'(6X,A,T45,I6)') 'Number of secondary orbitals',NSSHT
write(u6,'(6X,A,T45,F6.1)') 'Spin quantum number',Half*real(ISPIN-1,kind=wp)
write(u6,'(6X,A,T45,I6)') 'State symmetry',LSYM
write(u6,'(6X,A,T45,I6)') 'Number of configuration state fnc.',NCONF
write(u6,'(6X,A,T45,I6)') 'Number of root(s) available',NROOTS
write(u6,'(6X,A,T45,5I6)') 'CI root used',LROOT
if (ISCF == 0) then
  write(u6,'(6X,A)') 'This is a CASSCF reference function'
else if (ISCF == 1) then
  write(u6,'(6X,A)') 'This is a closed shell RHF reference function'
else
  write(u6,'(6X,A)') 'This is a high spin open shell RHF reference function'
end if
write(u6,*)
write(u6,*)
write(u6,'(6X,A)') 'Orbital specifications from previous RASSCF:'
write(u6,'(6X,A)') '--------------------------------------------'
write(u6,*)
write(u6,'(6X,A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Secondary orbitals',(nSsh(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
write(u6,*)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine PrInpPN
