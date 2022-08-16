!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine OpnFls_McKinley()

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep, lIrrep

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "disp.fh"
#include "disp2.fh"
#include "etwas.fh"
character*8 Method, MckLbl
character*288 Header

iOpt = 1
iRC = -1
MckLbl = 'Title'
call cWrMck(iRC,iOpt,MckLbl,1,Header,iDummer)
if (iRC /= 0) then
  write(6,*) 'OpnFls: Error writing to MCKINT'
  write(6,'(A,A)') 'MckLbl=',MckLbl
  call Abend()
end if
iOpt = 1
iRC = -1
MckLbl = 'nSym'
call iWrMck(iRC,iOpt,MckLbl,1,[nIrrep],iDummer)
if (iRC /= 0) then
  write(6,*) 'OpnFls: Error writing to MCKINT'
  write(6,'(A,A)') 'MckLbl=',MckLbl
  call Abend()
end if
iOpt = 0
iRC = -1
MckLbl = 'nBas'
call iWrMck(iRC,iOpt,MckLbl,1,nBas,iDummer)
if (iRC /= 0) then
  write(6,*) 'OpnFls: Error writing to MCKINT'
  write(6,'(A,A)') 'MckLbl=',MckLbl
  call Abend()
end if
iOpt = 0
iRC = -1
MckLbl = 'SymOp'
call cWrMck(iRC,iOpt,MckLbl,1,lirrep(0),iDummer)
if (iRC /= 0) then
  write(6,*) 'OpnFls: Error writing to MCKINT'
  write(6,'(A,A)') 'MckLbl=',MckLbl
  call Abend()
end if
iOpt = 0
iRC = -1
MckLbl = 'ldisp'
call iWrMck(iRC,iOpt,MckLbl,1,ldisp,iDummer)
if (iRC /= 0) then
  write(6,*) 'OpnFls: Error writing to MCKINT'
  write(6,'(A,A)') 'MckLbl=',MckLbl
  call Abend()
end if
ngrad = 0
do i=0,nIrrep-1
  nGrad = nGrad+ldisp(i)
end do
iOpt = 0
iRC = -1
MckLbl = 'chdisp'
call cWrMck(iRC,iOpt,MckLbl,1,chdisp(1),iDummer)
if (iRC /= 0) then
  write(6,*) 'OpnFls: Error writing to MCKINT'
  write(6,'(A,A)') 'MckLbl=',MckLbl
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the method label

call Get_cArray('Relax Method',Method,8)
if (Method == 'RHF-SCF ') then
  nMethod = SCF
else if (Method == 'CASSCF  ') then
  nMethod = RASSCF
else if (Method == 'CASSCFSA') then
  nMethod = RASSCF
  call Get_iScalar('SA ready',iGo)
  if (lHss .and. (iGo /= 2)) then
    write(6,*)
    write(6,*) ' Wavefunction type: RASSCF-SA'
    write(6,*)
    write(6,*) ' This option is not allowed when computing the Hessian. Use the RHS option!'
    call Quit_OnUserError()
  end if
else
  write(6,*) ' OpnFls: Wavefunction type:',Method
  write(6,*) '         Illegal type of wave function!'
  write(6,*) '         McKinley can not continue'
  write(6,*)
  call Quit_OnUserError()
end if

return

end subroutine OpnFls_McKinley
