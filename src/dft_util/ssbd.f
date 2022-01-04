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
* Copyright (C) 2001, Roland Lindh                                     *
*               2010, Grigory A. Shamov                                *
*               2017, Giovanni Li Manni                                *
*               2017, Aron Cohen                                       *
************************************************************************
      Subroutine SSBD(mGrid,Rho,nRho,
     &                iSpin,F_xc)
************************************************************************
*                                                                      *
* Object:     Combination of exchange SSB-D and PBE correlation terms  *
*             ref (secondary): This is true SSB-D                      *
*                                                                      *
*      Author: G. Li Manni A. Cohen, Max Planck Institute Stuttgart    *
*              Summer 2017, edited in Cambridge (UK) & Palermo (Sicily)*
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),
     &       F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*---- SSB-D Exchange -- unlike OPTX, SSB-D has its LDA part included !
      Coeff=1.0d0*CoefX
      Call xSSBD(mGrid,
     &          Coeff,iSpin,F_xc)

*---- CSPBE Correlation
      Coeff=1.0d0*CoefR
      Call CSPBE(mGrid,
     &         Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_Integer(nRho)
         Call Unused_real_array(Rho)
      End If
      End
