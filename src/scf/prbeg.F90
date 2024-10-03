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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine PrBeg(Meth)
!***********************************************************************
!                                                                      *
!     purpose: Print the header to iterations                          *
!                                                                      *
!***********************************************************************

use InfSCF, only: iDummy_run, InVec, jPrint, nD, nIter, nIterP, SCF_FileOrb
use Definitions, only: u6

implicit none
character(len=*), intent(in) :: Meth
character(len=10) :: Label
character(len=4) :: cUHF

if (jPrint < 2) return

write(u6,*)
call CollapseOutput(1,'Convergence information')
iDummy_run = 0
cUHF = '    '
if (nD == 2) cUHF = 'UHF '
Label = Meth(1:10)
if (nIter(nIterP) > 0) then
  write(u6,'(31x,A,A,A)') cUHF,Label,' iterations: Energy and convergence statistics'
  write(u6,*)
  write(u6,'(A,A,A)') 'Iter     Tot. ',Label, &
                      ' One-elec.       Two-elec.     Energy      Max Dij or  Max Fij      DNorm      TNorm      AccCon     Time'
  write(u6,'(A)') '         Energy          Energy          Energy        Change      Delta Norm'// &
                  '                                                in Sec.'
else
  iDummy_run = 1
  write(u6,'(45x,A)') 'No optimization is performed'
  if (InVec == 1) then
    write(u6,'(29x,A)') 'Results refer to orbitals obtained from core diagonalization'
  else if (InVec == 2) then
    write(u6,'(34x,A,A)') 'Results refer to input orbitals read from ',trim(SCF_FileOrb)
  else if (InVec == 3) then
    write(u6,'(34x,A)') 'Results refer to density matrix read from COMOLD'
  end if
end if

end subroutine PrBeg
