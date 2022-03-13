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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  HelState
!
!> @brief
!>   Couple the electrostatic part of the solvent with the QM-region.
!>   Only include the static part, no polarization at this moment
!> @author A. Ohrn
!>
!> @note
!> The quadrupoles are put in 'Buckingham-style'.
!>
!> @details
!> Rather easy to follow. This subroutine is a slightly modified
!> copy of hel.f. The interesting quantities are collected in
!> \p Vmat and are later to be added to the 'RASSI-matrix'.
!>
!> @param[in]  Eint    Field from static part of solvent on the Qm-molecule centers
!> @param[in]  nrstate Number of states in RASSI
!> @param[in]  ici     Number of MME-centers
!> @param[in]  Cha     Charges
!> @param[in]  Dip     Dipoles
!> @param[in]  Qua     Quadrupoles
!> @param[out] Vmat    The electrostatic part of the solute-solvent interaction matrix
!> @param[in]  iPrint  Print-level
!***********************************************************************

subroutine HelState(Eint,nrstate,ici,Cha,Dip,Qua,Vmat,iPrint)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
real(kind=wp) :: Eint(MxQCen,10), Cha(MxStOT,MxQCen), Dip(MxStOT,3,MxQCen), Qua(MxStOT,6,MxQCen), Vmat(MxOt)
integer(kind=iwp) :: nrstate, ici, iPrint
integer(kind=iwp) :: i, j, k, kaunt

kaunt = 0
do i=1,nrState
  do j=1,i
    kaunt = kaunt+1
    Vmat(kaunt) = Zero
  end do
end do

! The interaction between the distributed multipoles
! and the generalized field from the solvent.
kaunt = 0
do i=1,nrState
  do j=1,i
    kaunt = kaunt+1
    do k=1,ici
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,1)*Cha(kaunt,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,2)*Dip(kaunt,1,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,3)*Dip(kaunt,2,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,4)*Dip(kaunt,3,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,5)*Qua(kaunt,1,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,7)*Qua(kaunt,3,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,10)*Qua(kaunt,6,k)
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,6)*Qua(kaunt,2,k)*Two
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,8)*Qua(kaunt,4,k)*Two
      Vmat(kaunt) = Vmat(kaunt)+Eint(k,9)*Qua(kaunt,5,k)*Two
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iPrint)

end subroutine HelState
