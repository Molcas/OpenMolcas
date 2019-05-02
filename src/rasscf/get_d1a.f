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
!>    Transform the active one-body density from MO to AO basis
!>
!>  @author
!>    Markus P. Fuelscher
!>
!>  @param[in] CMO The occupation number vector in MO-space.
!>  @param[in] D1A_MO The active one-body density matrix in MO-space
!>  @param[out] D1A_AO The active one-body density matrix in AO-space
      Subroutine Get_D1A_RASSCF(CMO,D1A_MO,D1A_AO)
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
      implicit none
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
      real(8), intent(in) :: CMO(*) , D1A_MO(*)
      real(8), intent(out) :: D1A_AO(*)
      real(8), parameter :: Zero = 0.0d0
      integer :: iOff1, iOff2, iOff3, iSym, iBas, iAsh, iIsh, iFro,
     &    iTmp1, iTmp2

      Call qEnter('Get_D1A')

      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iFro = nFro(iSym)
        Call dCopy_(iBas*iBas,[Zero],0,D1A_AO(iOff3),1)
        If ( iAsh.ne.0 ) then
          Call GetMem('Scr1','Allo','Real',iTmp1,iAsh*iAsh)
          Call GetMem('Scr2','Allo','Real',iTmp2,iAsh*iBas)
          Call Square(D1A_MO(iOff1),Work(iTmp1),1,iAsh,iAsh)
          Call DGEMM_('N','T',
     &                iBas,iAsh,iAsh,
     &                1.0d0,CMO(iOff2+(iFro+iIsh)*iBas),iBas,
     &                Work(iTmp1),iAsh,
     &                0.0d0,Work(iTmp2),iBas)
          Call DGEMM_('N','T',
     &                iBas,iBas,iAsh,
     &                1.0d0,Work(iTmp2),iBas,
     &                CMO(iOff2+(iFro+iIsh)*iBas),iBas,
     &                0.0d0,D1A_AO(iOff3),iBas)
          Call GetMem('Scr2','Free','Real',iTmp2,iAsh*iBas)
          Call GetMem('Scr1','Free','Real',iTmp1,iAsh*iAsh)
        End If
        iOff1 = iOff1 + (iAsh*iAsh+iAsh)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + iBas*iBas
      End Do

      Call qExit('Get_D1A')

      Return
      End
