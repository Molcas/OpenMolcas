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
      Subroutine Get_D1I_RASSCF(CMO, D1I_AO)
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
      use general_data, only : nBas, nSym, nFro, nIsh
      implicit none
      real(8), intent(in) :: CMO(*)
      real(8), intent(out) :: D1I_AO(*)
      real(8), parameter :: Zero = 0.0d0, Two = 2.0d0
      integer :: ista, iSym, nb, nbsq, nfi

      Call qEnter('Get_D1I')

      ista=1
      do isym=1,nsym
        nb=nbas(isym)
        nbsq=nb**2
        nfi=nfro(isym)+nish(isym)
        if(nb.gt.0) then
          call dcopy_(nbsq,[0.0d0],0,d1i_AO(ista),1)
          if(nfi.gt.0) then
            call DGEMM_('n','t',nb,nb,nfi,two,
     &           cmo(ista),nb,cmo(ista),nb,zero,d1i_AO(ista),nb)
          end if
          ista=ista+nbsq
        end if
      end do

      Call qExit('Get_D1I')

      Return
      End
