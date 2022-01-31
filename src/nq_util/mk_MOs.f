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
* Copyright (C) 2022, Roland Lindh                                     *
************************************************************************
      Subroutine mk_MOs(SOValue,mAO,nCoor,MOValue,nMOs,CMOs,nCMO)
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 SOValue(mAO*nCoor,nMOs),
     &       MOValue(mAO*nCoor,nMOs), CMOs(nCMO)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*80 Label
#endif
*
#ifdef _DEBUGPRINT_
      Write (6,*) 'mk_MOs: MO-Coefficients'
      iOff=1
      Do iIrrep = 0, nIrrep-1
         If (nBas(iIrrep).gt.0) Then
            Write (6,*) ' Symmetry Block',iIrrep
            Call RecPrt(' ',' ',CMOs(iOff),nBas(iIrrep),nBas(iIrrep))
         End If
         iOff=iOff+nBas(iIrrep)**2
      End Do
#endif
*
*---- Compute some offsets
*
      iSO=1
      iCMO=1
      Do iIrrep = 0, nIrrep-1
         If (nBas(iIrrep)==0) Cycle
         Call DGeMM_('N','N',
     &               mAO*nCoor,nBas(iIrrep),nBas(iIrrep),
     &               One,SOValue(:,iSO:),mAO*nCoor,
     &                   CMOs(iCMO:),nBas(iIrrep),
     &               Zero,MOValue(:,iSO:),mAO*nCoor)
         iSO =iSO +nBas(iIrrep)
         iCMO=iCMO+nBas(iIrrep)*nBas(iIrrep)
      End Do
*
#ifdef _DEBUGPRINT_
      Write (Label,'(A)')'mk_MOs: MOValue(mAO*nCoor,nMOs)'
      Call RecPrt(Label,' ',MOValue(1,1),mAO*nCoor,nMOs)
#endif
*
      Return
      End
