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
!>    Transform the active one-body density from MO to AO basis
!>
!>  @author
!>    Markus P. Fuelscher
!>
!>  @details
!>  The underlying equation is the basis transformatin:
!>  \f[ D^{\text{AO}} = C D C^\dagger \f]
!>  For the aktive orbitals this becomes:
!>  \f[ D^{\text{AO}, A} = C^A D^A (C^A)^\dagger \f]
!>  Where (\f$ C^A, D^A \f$) are the coefficients and densities of the active MOs.
!>
!>  @param[in] CMO The MO-coefficients
!>  @param[in] D1A_MO The active one-body density matrix in MO-space
!>  @param[out] D1A_AO The active one-body density matrix in AO-space
subroutine Get_D1A_RASSCF(CMO,D1A_MO,D1A_AO)
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

use Index_Functions, only: nTri_Elem
use general_data, only: nAsh, nBas, nFro, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*), D1A_MO(*)
real(kind=wp), intent(_OUT_) :: D1A_AO(*)
integer(kind=iwp) :: iAsh, iBas, iFro, iIsh, iOff1, iOff2, iSym
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)

iOff1 = 1
iOff2 = 1
do iSym=1,nSym
  iBas = nBas(iSym)
  iAsh = nAsh(iSym)
  iIsh = nIsh(iSym)
  iFro = nFro(iSym)
  D1A_AO(iOff2:iOff2+iBas*iBas-1) = Zero
  if (iAsh /= 0) then
    call mma_allocate(Tmp1,iAsh*iAsh,Label='Tmp1')
    call mma_allocate(Tmp2,iAsh*iBas,Label='Tmp2')
    call Square(D1A_MO(iOff1),Tmp1,1,iAsh,iAsh)
    call DGEMM_('N','T',iBas,iAsh,iAsh, &
                One,CMO(iOff2+(iFro+iIsh)*iBas),iBas, &
                Tmp1,iAsh, &
                Zero,Tmp2,iBas)
    call DGEMM_('N','T',iBas,iBas,iAsh, &
                One,Tmp2,iBas, &
                CMO(iOff2+(iFro+iIsh)*iBas),iBas, &
                Zero,D1A_AO(iOff2),iBas)
    call mma_deallocate(Tmp2)
    call mma_deallocate(Tmp1)
  end if
  iOff1 = iOff1+nTri_Elem(iAsh)
  iOff2 = iOff2+iBas*iBas
end do

end subroutine Get_D1A_RASSCF
