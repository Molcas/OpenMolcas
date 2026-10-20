!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

function energy_mcwfn(dm1,h1e,vj,energy_nuc,n)

use PrintLevel, only: DEBUG
use mcpdft_output, only: iPrGlb
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: energy_mcwfn
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: dm1(n), h1e(n), vj(n), energy_nuc
real(kind=wp) :: e_j, te_vne
real(kind=wp), external :: ddot_

te_vne = dDot_(n,h1e,1,dm1,1)
e_j = Half*dDot_(n,vj,1,dm1,1)

if (iPrGlb >= DEBUG) then
  write(u6,*) 'Nuclear Repulsion energy',energy_nuc
  write(u6,*) 'Te_Vne',te_vne
  write(u6,*) 'E_j',e_j
end if

energy_mcwfn = energy_nuc+te_vne+e_j

end function energy_mcwfn
