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

subroutine ExtractS(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nBas,xyzQu,lExtr,iExtr_Atm,ip_ExpVal,ip_ExpCento,ENR,ENP)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
dimension xyzMy(3), Hmat(*), xyzQu(6)
dimension iDt(3), iExtr_Atm(*)
logical lExtr(*)

write(iLu,*) '<<<<<<<Configuration ',i9,'>>>>>>>'
if (lExtr(1)) then
  write(iLu,*) 'Total Energy'
  write(iLu,'(F15.8)') Etot
end if
if (lExtr(2)) then
  write(iLu,*) 'QM-Dipole'
  write(iLu,'(3(F12.5))') (xyzMy(k),k=1,3)
end if
if (lExtr(3)) then
  write(iLu,*) 'QM-Quadrupole'
  write(iLu,'(6(F12.5))') (xyzQu(k),k=1,6)
end if
if (lExtr(6)) then
  write(iLu,*) 'Expectation values (T+H_nuc,V_el,V_pol,V_pp)'
  write(iLu,*) '  Nuc cont:',ENR
  write(iLu,'(4(F15.8))') (Work(ip_ExpVal+k),k=0,3)
  call GetMem('ExpVals','Free','Real',ip_ExpVal,4)
end if
if (lExtr(7)) then
  write(iLu,*) 'Expectation values partial V_el, V_pol'
  write(iLu,*) '  Nuc cont:',ENP
  write(iLu,'(2(F15.8))') (Work(ip_ExpCento+k),k=1,2)
  call GetMem('ExpVals','Free','Real',ip_ExpCento,4)
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Hmat)
  call Unused_integer(iC)
  call Unused_integer_array(iDt)
  call Unused_integer(nBas)
  call Unused_integer_array(iExtr_Atm)
end if

end subroutine ExtractS
