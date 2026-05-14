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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

!>  @brief
!>  Transform the inactive one-body density from MO to AO basis
!>
!>  @author
!>  Markus P. Fuelscher
!>
!>  @details
!>  The underlying equation is the basis transformatin:
!>  \f[ D^{\text{AO}} = C D C^\dagger \f]
!>  For inactive orbitals it simplifies to:
!>  \f[D^{\text{AO}, I} = C^I D^I (C^I)^\dagger = 2 C^I \mathbf{1} (C^I)^\dagger = 2 C (C^I)^\dagger \f]
!>  Where (\f$ C^I, D^I \f$) are the coefficients and densities of the inactive MOs.
!>
!>  @param[in] CMO The MO-coefficients
!>  @param[out] D1I_AO The inactive one-body density matrix in AO-space
subroutine Get_D1I_RASSCF(CMO,D1I_AO)
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use general_data, only: nBas, nFro, nIsh, nSym
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(_OUT_) :: D1I_AO(*)
integer(kind=iwp) :: ista, iSym, nb, nbsq, nfi

ista = 1
do isym=1,nsym
  nb = nbas(isym)
  nbsq = nb**2
  nfi = nfro(isym)+nish(isym)
  if (nb > 0) then
    d1i_AO(ista:ista+nbsq-1) = Zero
    if (nfi > 0) call DGEMM_('N','T',nb,nb,nfi, &
                             Two,cmo(ista),nb, &
                             cmo(ista),nb, &
                             Zero,d1i_AO(ista),nb)
    ista = ista+nbsq
  end if
end do

end subroutine Get_D1I_RASSCF
