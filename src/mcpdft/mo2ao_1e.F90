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
!  mo2act_1e
!
!> @brief extracts the active block from folded MO one-particle integrals
!>
!> @param[in] d_mo folded 1 particle integral matrix in MO basis
!> @param[in] nSym number of symmetry groups
!> @param[in] nOrb number of orbitals per symmetry group
!> @param[in] nAsh number of active orbitals per symmetry group
!> @param[in] nIsh number of inactive orbitals per symmetry group
!>
!> @param[out] d_act folded 1 particle integral matrix in active MO basis
!***********************************************************************

subroutine mo2act_1e(d_mo,d_act,nSym,nOrb,nAsh,nIsh)

use Index_Functions, only: iTri, nTri_Elem
use general_data, only: ntot1
use rasscf_global, only: nacpar
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: d_mo(ntot1)
real(kind=wp), intent(out) :: d_act(nacpar)
integer(kind=iwp), intent(in) :: nsym, norb(nsym), nash(nsym), nish(nsym)
integer(kind=iwp) :: iact, iactoff, imo, imoff, iash, iorb, isym, jact, jmo

d_act(:) = Zero
iActOff = 1
iMOff = 1
do iSym=1,nSym
  iOrb = nOrb(iSym)
  iAsh = nAsh(iSym)
  if ((iOrb == 0) .or. (iAsh == 0)) then
    iMOff = iMOff+nTri_Elem(iOrb)
    cycle
  end if

  do iAct=1,iAsh
    iMO = nIsh(iSym)+iAct
    do jAct=1,iAct
      jMO = nIsh(iSym)+jAct
      d_act(iActOff) = d_mo(iMOff+iTri(iMO,jMO)-1)
      iActOff = iActOff+1
    end do
  end do
  iMOff = iMOff+nTri_Elem(iOrb)
end do

end subroutine mo2act_1e
