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

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
real*8 h0(*), Ei(*), S(*)

nInts = nBas*(nBas+1)/2
call Allocate_Work(ip_h0,nInts+4)
call Comp_F_(h0,Ei,nBas,Delta_i,Energy,Work(ip_h0),S,Refx,Originx,nInts)
call Free_Work(ip_h0)

return

end subroutine Comp_F

subroutine Comp_F_(h0,Ei,nBas,Delta_i,Energy,h0_temp,S,Refx,Originx,nInts)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
real*8 h0(nInts+4), h0_temp(nInts+4), Ei(nInts+4), S(nInts+4)
character*8 Method, Label
integer mBas(8)

call Get_cArray('Relax Method',Method,8)
call Allocate_Work(ipC,1)
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
call WrOne(iRc,0,Label,iComp,h0_temp,iSyLbl)
!call TriPrt('H0_temp after wrone',' ',h0_temp,nBas)

if (Method == 'RHF-SCF' .or. Method == 'UHF-SCF' .or. Method == 'KS-DFT') then

  call StartLight('scf')
  call Disable_Spool()
  call xml_open('module',' ',' ',0,'scf')
  call SCF(ireturn)
  call xml_close('module')
  if (iReturn /= 0) Go To 99
  call GetMem('PT2','Flush','Real',ipC,iDum)

  else if (Method(1:5) == 'MBPT2') then

  call StartLight('scf')
  call Disable_Spool()
  call xml_open('module',' ',' ',0,'scf')
  call SCF(ireturn)
  call xml_close('module')
  if (iReturn /= 0) Go To 99
  call GetMem('PT2','Flush','Real',ipC,iDum)
  call StartLight('mbpt2')
  call Disable_Spool()
  call mp2_driver(ireturn)
  if (iReturn /= 0) Go To 99
  call GetMem('PT2','Flush','Real',ipC,iDum)

else if (Method == 'RASSCF' .or. Method == 'CASSCF') then

  call StartLight('rasscf')
  call Disable_Spool()
  call RASSCF(ireturn)
  if (iReturn /= 0) Go To 99
  call GetMem('PT2','Flush','Real',ipC,iDum)

else if (Method == 'CASPT2') then

  call StartLight('rasscf')
  call Disable_Spool()
  call RASSCF(ireturn)
  if (iReturn /= 0) Go To 99
  call GetMem('PT2','Flush','Real',ipC,iDum)
  call StartLight('caspt2')
  call Disable_Spool()
  call CASPT2(ireturn)
  if (iReturn /= 0) Go To 99
  call GetMem('PT2','Flush','Real',ipC,iDum)

else
  write(6,*) 'Method=',Method
  write(6,*) ' Oups!'
  call Abend()
end if

!call Get_Energy(Energy)
call Get_dScalar('Last energy',Energy)
call Free_Work(ipC)

call WrOne(iRc,0,Label,iComp,h0,iSyLbl)
call Put_dScalar('PotNuc',PotNuc_Save)

return

99 continue
write(6,*)
write(6,*) 'Comp_f: Wave function calculation failed!'
write(6,*)
call Abend()
! Avoid unused argument warnings
if (.false.) call Unused_integer(nBas)

end subroutine Comp_F_
