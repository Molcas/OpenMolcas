************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************

*>  @brief
*>    Transform the active one-body density from MO to AO basis
*>
*>  @author
*>    Markus P. Fuelscher
*>
*>  @details
*>  The underlying equation is the basis transformatin:
*>  \f[ D^{\text{AO}} = C D C^\dagger \f]
*>  For the aktive orbitals this becomes:
*>  \f[ D^{\text{AO}, A} = C^A D^A (C^A)^\dagger \f]
*>  Where (\f$ C^A, D^A \f$) are the coefficients and densities of the active MOs.
*>
*>  @param[in] CMO The MO-coefficients
*>  @param[in] D1A_MO The active one-body density matrix in MO-space
*>  @param[out] D1A_AO The active one-body density matrix in AO-space
      subroutine Get_D1A_RASSCF(CMO,D1A_MO,D1A_AO)
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use general_data, only : nBas, nSym, nFro, nIsh, nAsh
      use stdalloc, only: mma_allocate, mma_deallocate
      implicit none
      real*8, intent(in) :: CMO(*) , D1A_MO(*)
      real*8, intent(out) :: D1A_AO(*)
      real*8, parameter :: Zero = 0.0d0
      integer :: iOff1, iOff2, iOff3, iSym, iBas, iAsh, iIsh, iFro
      Real*8, Allocatable:: Tmp1(:), Tmp2(:)

      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iFro = nFro(iSym)
        Call dCopy_(iBas*iBas,[Zero],0,D1A_AO(iOff3),1)
        If ( iAsh /= 0 ) then
          Call mma_allocate(Tmp1,iAsh*iAsh,Label='Tmp1')
          Call mma_allocate(Tmp2,iAsh*iBas,Label='Tmp2')
          Call Square(D1A_MO(iOff1),Tmp1,1,iAsh,iAsh)
          Call DGEMM_('N','T',
     &                iBas,iAsh,iAsh,
     &                1.0d0,CMO(iOff2+(iFro+iIsh)*iBas),iBas,
     &                      Tmp1,iAsh,
     &                0.0d0,Tmp2,iBas)
          Call DGEMM_('N','T',
     &                iBas,iBas,iAsh,
     &                1.0d0,Tmp2,iBas,
     &                      CMO(iOff2+(iFro+iIsh)*iBas),iBas,
     &                0.0d0,D1A_AO(iOff3),iBas)
          Call mma_deallocate(Tmp2)
          Call mma_deallocate(Tmp1)
        End If
        iOff1 = iOff1 + (iAsh*iAsh+iAsh)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + iBas*iBas
      End Do

      End subroutine Get_D1A_RASSCF
