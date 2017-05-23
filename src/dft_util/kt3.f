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
     &                dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
     &                T_X)
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
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
C      Call QEnter('OLYP')
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange with the non-UEG factor!
*
      Coeff= 1.092d0
      Call Diracx(mGrid,Rho,nRho,iSpin,F_xc,
     &            dF_dRho,ndF_dRho,Coeff,T_X)
*
*---- OPTX Exchange
*
      Coeff= 0.925452d0
      Call xOPT(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          Coeff,iSpin,F_xc,T_X)
*
*---- KT term Exchange
*
      Coeff= -0.0040d0
      Call KealTozer(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          Coeff,iSpin,F_xc,T_X)
*
*---- Lee-Yang-Parr Correlation
*
      Coeff=0.864409d0
      Call LYP(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &         Coeff,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
C      Call QExit('OLYP')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
