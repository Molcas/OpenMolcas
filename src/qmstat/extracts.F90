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

subroutine ExtractS(iLu,i9,Etot,xyzMy,xyzQu,lExtr,ExpVal,ExpCento,ENR,ENP)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iLu, i9
real(kind=wp), intent(in) :: Etot, xyzMy(3), xyzQu(6), ExpVal(4,1), ExpCento(4,1), ENR, ENP
logical(kind=iwp), intent(in) :: lExtr(*)

write(iLu,*) '<<<<<<<Configuration ',i9,'>>>>>>>'
if (lExtr(1)) then
  write(iLu,*) 'Total Energy'
  write(iLu,'(F15.8)') Etot
end if
if (lExtr(2)) then
  write(iLu,*) 'QM-Dipole'
  write(iLu,'(3(F12.5))') xyzMy(:)
end if
if (lExtr(3)) then
  write(iLu,*) 'QM-Quadrupole'
  write(iLu,'(6(F12.5))') xyzQu(:)
end if
if (lExtr(6)) then
  write(iLu,*) 'Expectation values (T+H_nuc,V_el,V_pol,V_pp)'
  write(iLu,*) '  Nuc cont:',ENR
  write(iLu,'(4(F15.8))') ExpVal(:,1)
end if
if (lExtr(7)) then
  write(iLu,*) 'Expectation values partial V_el, V_pol'
  write(iLu,*) '  Nuc cont:',ENP
  write(iLu,'(2(F15.8))') ExpCento(:,1)
end if

return

end subroutine ExtractS
