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

subroutine PrtCav(IOut,ITyp,NS,NOrd,Alpha,Rad)
! Print out sphere radii for Pauling or input cavitites

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut, ITyp, NS, NOrd(NS)
real(kind=wp) :: Alpha, Rad(NS)
integer(kind=iwp) :: IS

write(iOut,*)
write(iOut,*)
write(iOut,'(6X,A)') 'Polarized Continuum Model Cavity'
write(iOut,'(6X,A)') '================================'
if (ITyp == 2) write(iOut,'(6X,A)') 'Pauling radii'
if (ITyp == 3) write(iOut,'(6X,A)') 'Sphere radii from input'
write(iOut,*)
write(IOut,'(6X,A)') ' NOrd  Alpha  Radius'
do IS=1,NS
  write(IOut,'(6X,I5,2X,F5.2,2X,F6.3)') NOrd(IS),Alpha,Rad(IS)
end do
write(IOut,'(6X,1X,78("-"))')
write(IOut,*)

return

end subroutine PrtCav
