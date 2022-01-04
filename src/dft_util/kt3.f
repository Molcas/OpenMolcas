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
* Copyright (C) 2009, Grigory A. Shamov                                *
************************************************************************
      Subroutine KT3(mGrid,Rho,nRho,P2_ontop,
     &                nP2_ontop,iSpin,F_xc,
     &                dF_dP2ontop,ndF_dP2ontop)
************************************************************************
*                                                                      *
* Object:     KT3 combination                                          *
*            ref: Keal, TW, Tozer, DJ. J.Chem.Phys. 121 (2004) 5654    *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
* Author:      Grigory A Shamov, U of Manitoba, 2009                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange with the non-UEG factor!
*
      Coeff= 1.092d0*CoefX
      Call Diracx(mGrid,iSpin,F_xc,Coeff)
*
*---- OPTX Exchange
*
      Coeff= 0.925452d0*CoefX
      Call xOPT(mGrid,
     &          Coeff,iSpin,F_xc)
*
*---- KT term Exchange
*
      Coeff= -0.0040d0*CoefX
      Call KealTozer(mGrid,
     &          Coeff,iSpin,F_xc)
*
*---- Lee-Yang-Parr Correlation
*
      Coeff=0.864409d0*CoefR
      Call LYP(mGrid,
     &         Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_Integer(nRho)
         Call Unused_real_array(Rho)
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
