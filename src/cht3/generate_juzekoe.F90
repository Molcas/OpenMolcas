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

subroutine generate_juzekOE(oe,oeh,oep,no,nv)
! this routine does:
!
! modifies standard oe record to one used by DIRCC
!
! oeh = ( oe_occ(alpha) ... oe_occ(beta) )   (alpha=beta for this implementation)
!
! oep = ( oe_virt(alpha) ... oe_virt(beta) ) (alpha=beta for this implementation)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no, nv
real(kind=wp), intent(in) :: oe(no+nv)
real(kind=wp), intent(out) :: oeh(no,2), oep(nv,2)

oeh(:,1) = oe(1:no)
oeh(:,2) = oe(1:no)

oep(:,1) = oe(no+1:no+nv)
oep(:,2) = oe(no+1:no+nv)

return

end subroutine generate_juzekOE
