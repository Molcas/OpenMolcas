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

subroutine Last_Energy(iReturn)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
!#include "timtra.fh"
integer(kind=iwp) :: lengthlast
character(len=8) :: Method
logical(kind=iwp) :: Do_ESPF, StandAlone, FoundLastEn
!                                                                      *
!***********************************************************************
!                                                                      *
iReturn = 99

! Get information regarding the last method used

call qpg_cArray('LastEnergyMethod',FoundLastEn,lengthlast)
if (FoundLastEn) then
  call Get_cArray('LastEnergyMethod',Method,8)
else
  call Get_cArray('Relax Method',Method,8)
end if

call DecideOnESPF(Do_ESPF)

if (Method(5:7) == 'SCF'     .or. &
    Method(1:6) == 'KS-DFT'  .or. &
    Method(1:6) == 'CASSCF'  .or. &
    Method(1:6) == 'RASSCF'  .or. &
    Method(1:6) == 'CASPT2'  .or. &
    Method(1:5) == 'MBPT2'   .or. &
    Method(1:5) == 'CCSDT'   .or. &
    Method(1:4) == 'CHCC'    .or. &
    Method(1:6) == 'MCPDFT'  .or. &
#   ifdef _DMRG_
    Method(1:7) == 'DMRGSCF' .or. &
#   endif
    Method(1:4) == 'CHT3'    .or. &
    Method(1:8) == 'EXTERNAL') then
  continue
else
  write(u6,'(A,A,A)') 'Last Energy for ',Method,' is not implemented yet.'
  call Abend()
end if

if (Method(1:6) == 'MCPDFT') Do_ESPF = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute integrals

call StartLight('seward')
call Disable_Spool()
call Seward(iReturn)
if (iReturn /= 0) then
  write(u6,*) 'Last_Energy failed ...'
  write(u6,*) 'SEWARD returned with return code, rc = ',iReturn
  call Abend()
end if

! Compute ESPF

if (Do_ESPF) then
  call StartLight('espf')
  call Disable_Spool()
  StandAlone = .true.
  call ESPF(iReturn,StandAlone)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'ESPF returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

! Compute the wave function

if ((Method(5:7) == 'SCF' .and. Method(1:4) /= 'DMRG') .or. &
    Method(1:6) == 'KS-DFT' .or. &
    Method(1:5) == 'MBPT2'  .or. &
    Method(1:4) == 'CHCC'   .or. &
    Method(1:4) == 'CHT3') then
  call StartLight('scf')
  call Disable_Spool()
  call xml_open('module',' ',' ',0,'scf')
  call SCF(iReturn)
  call xml_close('module')
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'SCF returned with return code, rc = ',iReturn
    call Abend()
  end if
else if (Method(1:6) == 'RASSCF' .or. &
         Method(1:6) == 'CASSCF' .or. &
         Method(1:6) == 'CASPT2' .or. &
         Method(1:6) == 'MCPDFT' .or. &
         Method(1:5) == 'CCSDT') then
  call StartLight('rasscf')
  call Disable_Spool()
  call RASSCF(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'RASSCF returned with return code, rc = ',iReturn
    call Abend()
  end if
#ifdef _DMRG_
else if (Method(1:7) == 'DMRGSCF') then
  call StartLight('dmrgscf')
  call Disable_Spool()
  call DMRGSCF(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'DMRGSCF returned with return code, rc = ',iReturn
    call Abend()
  end if
#endif
else if (Method(1:8) == 'EXTERNAL') then
  call StartLight('false')
  call Disable_Spool()
  call False_program(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'FALSE returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

if (Method(1:5) == 'MBPT2') then
  call StartLight('mbpt2')
  call Disable_Spool()
  call MP2_Driver(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'MBPT2 returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

if (Method(1:5) == 'CCSDT') then
  call StartLight('motra')
  call Disable_Spool()
  call MOTra(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'MOTra returned with return code, rc = ',iReturn
    call Abend()
  end if

  call StartLight('ccsdt')
  call Disable_Spool()
  call CCSDT(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'CCSDT returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

if (Method(1:4) == 'CHCC' .or. &
    Method(1:4) == 'CHT3') then
  call StartLight('chcc')
  call Disable_Spool()
  call CHCC(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'CHCC returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

if (Method(1:4) == 'CHT3') then
  call StartLight('cht3')
  call Disable_Spool()
  call CHT3(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'CHT3 returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

if (Method(1:6) == 'CASPT2') then
  call StartLight('caspt2')
  call Disable_Spool()
  call CASPT2(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'CASPT2 returned with return code, rc = ',iReturn
    call Abend()
  end if
end if

if (Method(1:6) == 'MCPDFT') then
  call StartLight('mcpdft')
  call Disable_Spool()
  call MCPDFT(iReturn)
  if (iReturn /= 0) then
    write(u6,*) 'Last_Energy failed ...'
    write(u6,*) 'MCPDFT returned with return code, rc = ',iReturn
    call Abend()
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Since Last_Energy itself cleans up after the modules, it's necessary
! to force nfld_tim and nfld_stat to zero, or finish will scream.
!nfld_tim = 0
!nfld_stat = 0
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Last_Energy
