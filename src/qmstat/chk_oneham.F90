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

subroutine Chk_OneHam(nBas)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "WrkSpc.fh"
dimension nBas(MxSym)
character Label_Read*8, Label_Pure*8

Lu_One = 49
Lu_One = IsFreeUnit(Lu_One)
Label_Read = 'OneHam  '
Label_Pure = 'OneHam 0'
nBT = nBas(1)*(nBas(1)+1)/2
call OpnOne(irc,0,'ONEINT',Lu_One)
call GetMem('Read','Allo','Real',iOneR,nBT+4)
call GetMem('Pure','Allo','Real',iOneP,nBT+4)

irc = -1
iopt = 0
iSmLbl = 0
call RdOne(irc,iopt,Label_Read,1,Work(iOneR),iSmLbl)
irc = -1
iopt = 0
iSmLbl = 0
call RdOne(irc,iopt,Label_Pure,1,Work(iOneP),iSmLbl)
call ClsOne(irc,Lu_One)

call DaxPy_(nBT,-1.0d0,Work(iOneR),1,Work(iOneP),1)

dNorm = dnrm2_(nBT,Work(iOneP),1)

if (dNorm > 1d-8) then
  write(6,*)
  write(6,*)
  write(6,*) ' WARNING!'
  write(6,*)
  write(6,*) '   Your one-electron hamiltonian is not purely vacuum. This means that the Hamiltonian'
  write(6,*) '   in QmStat can be contaminated. Is this intentional? If not, then make sure that the ONEINT'
  write(6,*) '   file comes directly from a Seward calculation without any calls from'
  write(6,*) '   FFPT (or similar) in between.'
  write(6,*)
  write(6,*)
end if

call GetMem('Read','Free','Real',iOneR,nBT+4)
call GetMem('Pure','Free','Real',iOneP,nBT+4)

return

end subroutine Chk_OneHam
