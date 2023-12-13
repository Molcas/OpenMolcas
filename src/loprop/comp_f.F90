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

subroutine Comp_F(h0,Ei,nBas,Delta_i,Energy,S,Refx,Originx)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas
real(kind=wp), intent(in) :: h0(nBas*(nBas+1)/2+4), Ei(nBas*(nBas+1)/2+4), Delta_i, S(nBas*(nBas+1)/2+4), Refx, Originx
real(kind=wp), intent(out) :: Energy
integer(kind=iwp) :: i, iComp, iOpt, ireturn, iRc, iSyLbl, mBas(8), nInts, nsize
real(kind=wp) :: PotNuc_Save
character(len=8) :: Method, Label
real(kind=wp), allocatable :: h0_temp(:)

nInts = nBas*(nBas+1)/2
call mma_allocate(h0_temp,nInts+4,label='h0_temp')

call Get_cArray('Relax Method',Method,8)
call Get_iScalar('nSym',i)
call Get_iArray('nBas',mBas,i)
nsize = nInts+4

! Add perturbations to h0

call dcopy_(nSize,h0,1,h0_temp,1)
call DaXpY_(nSize,Delta_i,Ei,1,h0_temp,1)
call DaXpY_(nSize,Delta_i*(Originx-Refx),S,1,h0_temp,1)
call Get_dScalar('PotNuc',PotNuc_Save)
call Put_dScalar('PotNuc',h0_temp(nInts+4))

!---- Write the one-electron hamiltonian.

!call TriPrt('H0 before wrone',' ',h0,nBas)
iComp = 1
iSyLbl = 1
Label = 'OneHam  '
iRc = -1
iOpt = 0
call WrOne(iRc,iOpt,Label,iComp,h0_temp,iSyLbl)
!call TriPrt('H0_temp after wrone',' ',h0_temp,nBas)

if ((Method == 'RHF-SCF') .or. (Method == 'UHF-SCF') .or. (Method == 'KS-DFT')) then

  call StartLight('scf')
  call Disable_Spool()
  call xml_open('module',' ',' ',0,'scf')
  call SCF(ireturn)
  call xml_close('module')
  if (iReturn /= 0) call Error()

else if (Method(1:5) == 'MBPT2') then

  call StartLight('scf')
  call Disable_Spool()
  call xml_open('module',' ',' ',0,'scf')
  call SCF(ireturn)
  call xml_close('module')
  if (iReturn /= 0) call Error()
  call StartLight('mbpt2')
  call Disable_Spool()
  call mp2_driver(ireturn)
  if (iReturn /= 0) call Error()

else if ((Method == 'RASSCF') .or. (Method == 'CASSCF')) then

  call StartLight('rasscf')
  call Disable_Spool()
  call RASSCF(ireturn)
  if (iReturn /= 0) call Error()

else if (Method == 'CASPT2') then

  call StartLight('rasscf')
  call Disable_Spool()
  call RASSCF(ireturn)
  if (iReturn /= 0) call Error()
  call StartLight('caspt2')
  call Disable_Spool()
  call CASPT2(ireturn)
  if (iReturn /= 0) call Error()

else
  write(u6,*) 'Method=',Method
  write(u6,*) ' Oups!'
  call Abend()
end if

!call Get_Energy(Energy)
call Get_dScalar('Last energy',Energy)

call WrOne(iRc,0,Label,iComp,h0,iSyLbl)
call Put_dScalar('PotNuc',PotNuc_Save)

call mma_deallocate(h0_temp)

return

contains

subroutine Error()

  write(u6,*)
  write(u6,*) 'Comp_f: Wave function calculation failed!'
  write(u6,*)
  call Abend()

end subroutine Error

end subroutine Comp_F
