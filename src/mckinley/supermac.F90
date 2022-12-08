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

subroutine SuperMac()

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
#include "warnings.h"
#include "temperatures.fh"
integer(kind=iwp) :: i, iDNG, iErr, irlxroot, iSigma, LuInput, nData, nXF
character(len=16) :: StdIn
character(len=8) :: Method
logical(kind=iwp) :: Do_Cholesky, Do_ESPF, Numerical, Found
integer(kind=iwp), allocatable :: Scr1(:)
integer(kind=iwp), external :: IsFreeUnit

call Get_cArray('Relax Method',Method,8)

Numerical = (Method(1:6) == 'RASSCF') .or. (Method(1:6) == 'GASSCF') .or. (Method == 'CASSCFSA') .or. (Method == 'DMRGSCFS') .or. &
            (Method == 'CASPT2') .or. (Method == 'UHF-SCF') .or. (Method == 'MBPT2') .or. (Method == 'CCSDT') .or. &
            (Method == 'KS-DFT') .or. (Method == 'UKS-DFT') .or. (Method == 'MCPDFT') .or. (Method == 'MSPDFT')

if (Method == 'CASSCF') then
  call Get_iScalar('NumGradRoot',irlxroot)
  Numerical = irlxroot /= 1
end if

call DecideOnCholesky(Do_Cholesky)
if (Do_Cholesky) Numerical = .true.

call Qpg_iScalar('nXF',Found)
if (Found) then
  call Get_iScalar('nXF',nXF)
  Numerical = Numerical .or. (nXF > 0)
end if
call DecideOnESPF(Do_ESPF)
Numerical = Numerical .or. Do_ESPF

! Analytical PCM frequencies are currently not implemented

call Qpg_dArray('PCM info',Found,nData)
Numerical = Numerical .or. (Found .and. (nData > 0))

call Qpg_iScalar('DNG',Found)
if (Found) then
  call Get_iScalar('DNG',iDNG)
  Numerical = Numerical .or. (iDNG == 1)
end if

if (.not. Numerical) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Create a backup runfile before running the numerical differentiation

call fCopy('RUNFILE','RUNBACK',iErr)
if (iErr /= 0) call Abend()

call mma_allocate(Scr1,7,Label='Scr1')
Scr1(1) = 0
Scr1(2) = 0
Scr1(3) = -99
Scr1(4) = 0
call Put_iArray('Slapaf Info 1',Scr1,7)
call mma_deallocate(Scr1)
LuInput = IsFreeUnit(11)
call StdIn_Name(StdIn)
call Molcas_open(LuInput,StdIn)
!                                                                      *
!***********************************************************************
!                                                                      *
write(LuInput,'(A)') '>ECHO OFF'
write(LuInput,'(A)') '>export MCK_OLD_TRAP=$MOLCAS_TRAP'
write(LuInput,'(A)') '>export MCK_OLD_MAXITER=$MOLCAS_MAXITER'
write(LuInput,'(A)') '> export MOLCAS_TRAP=ON'
write(LuInput,'(A)') '> export MOLCAS_MAXITER=500'

! If SA-CASSCF run the MCLR code so that the reference dipole moment is variational.

if ((Method == 'RASSCFSA') .or. (Method == 'CASSCFSA')) write(LuInput,'(A)') '&MCLR'

write(LuInput,'(A)') '> DO WHILE <'
write(LuInput,'(A)') '> IF (ITER NE 1) <'

call Lu2Lu('SEWARINP',LuInput)
write(LuInput,*)

if (Do_ESPF) then
  call Lu2Lu('ESPFINP',LuInput)
end if

if ((Method == 'RASSCFSA') .or. (Method == 'CASSCFSA') .or. (Method == 'CASSCF')) then
  call Lu2Lu('RASSCINP',LuInput)
else if ((Method == 'MCPDFT') .or. (Method == 'MSPDFT')) then
  call Lu2Lu('RASSCINP',LuInput)
  write(LuInput,'(A)')
  call Lu2Lu('MCPDFINP',LuInput)
else if (Method == 'CASPT2') then
  call Lu2Lu('RASSCINP',LuInput)
  write(LuInput,'(A)')
  call Lu2Lu('CASPTINP',LuInput)
else if (Method == 'MBPT2') then
  call Lu2Lu('SCFINP',LuInput)
else if (Method == 'CCSDT') then
  call Lu2Lu('SCFINP',LuInput)
  write(LuInput,'(A)')
  call Lu2Lu('CCSDTINP',LuInput)
else if ((Method == 'KS-DFT') .or. (Method == 'RHF-SCF') .or. (Method == 'UKS-DFT') .or. (Method == 'UHF-SCF')) then
  call Lu2Lu('SCFINP',LuInput)
end if

write(LuInput,'(A)') '> END IF <'

! To make sure MBPT2 is run with the Grdt option, run it always (including first iteration, i.e. outside the IF)

if (Method == 'MBPT2') then
  write(LuInput,'(A)')
  call Lu2Lu('MBPT2INP',LuInput)
end if

write(LuInput,'(A)')
write(LuInput,'(A)') '&Slapaf &End'
write(LuInput,'(A)') 'Numerical'
write(LuInput,'(A)') 'Iterations'
write(LuInput,'(A)') '0'
write(LUInput,'(A)') 'THERmochemistry'

call Get_iScalar('Rotational Symmetry Number',iSigma)
write(LUInput,'(I3)') iSigma
write(LUInput,'(A)') '1.0'
do i=1,NDefTemp
  write(LUInput,'(F7.2)') DefTemp(i)
end do
write(LUInput,'(A)') 'End of PT'
write(LuInput,'(A)') 'End of Input'
write(LuInput,'(A)') '> END DO <'
write(LuInput,'(A)') '> export MOLCAS_TRAP=$MCK_OLD_TRAP'
write(LuInput,'(A)') '> export MOLCAS_MAXITER=$MCK_OLD_MAXITER'
write(LuInput,'(A)') '>ECHO ON'
close(LuInput)
!                                                                      *
!***********************************************************************
!                                                                      *
call Finish(_RC_INVOKED_OTHER_MODULE_)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine SuperMac
