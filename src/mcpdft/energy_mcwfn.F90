!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************
function energy_mcwfn(dm1,h1e,vj,energy_nuc,dim)
  use definitions,only:iwp,wp,u6
  use printlevel,only:debug
  use constants,only:half
  use mcpdft_output,only:iPrGlb
  implicit none

  integer(kind=iwp),intent(in) :: dim
  real(kind=wp),intent(in) :: dm1(dim),h1e(dim),vj(dim),energy_nuc
  real(kind=wp) :: energy_mcwfn

  real(kind=wp),external :: ddot_

  real(kind=wp) :: te_vne,e_j

  te_vne = dDot_(dim,h1e,1,dm1,1)
  e_j = half*dDot_(dim,vj,1,dm1,1)

  if(iPrGlb >= debug) then
    write(u6,*) 'Nuclear Repulsion energy',energy_nuc
    write(u6,*) 'Te_Vne',te_vne
    write(u6,*) 'E_j',e_j
  endif

  energy_mcwfn = energy_nuc+te_vne+e_j

endfunction
