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
*               2009, Grigory A. Shamov                                *
************************************************************************
      Subroutine B2PLYP(mGrid,Rho,nRho,P2_ontop,
     &                  nP2_ontop,iSpin,F_xc,
     &                  dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
************************************************************************
*                                                                      *
* Object:   DFT part for the Grimme's double hybrud functional.        *
*         requires (manual?) scaling of the &mbpt2 correlation energy! *
*         Grimme S., J. Chem.Phys 124(2006) 034108                     *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author: Template by Roland Lindh, Department of Chemical Physi  *
*             of Lund, SWEDEN. March 2001                              *
*              Grigory A Shamov, University of Manitoba, 2009          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
C      Call QEnter('B2PLYP')
*                                                                      *
************************************************************************
*                                                                      *
      Coeff_A=0.470D0
      Coeff_B=0.0000
      Coeff_C=0.73D0
C    requires 0.27 times MP2 correlation energy from MBPT2
*                                                                      *
*---- Dirac Exchange Functional                                        *
*                                                                      *
      Call Diracx(mGrid,Rho,nRho,
     &            iSpin,F_xc,dF_dRho,
     &            ndF_dRho,Coeff_A,T_X)
*                                                                      *
*---- Becke 88 Exchange Functional                                     *
*                                                                      *
      Call xB88(Rho,nRho,mGrid,
     &          dF_dRho,ndF_dRho,
     &          Coeff_A,iSpin,F_xc,T_X)
*                                                                      *
*---- Vosko-Wilk-Nusair Correlation Functional V                       *
*                                                                      *
C--      Call VWN_V(mGrid,Rho,nRho,
C--     &           iSpin,F_xc,dF_dRho,
C--     &           ndF_dRho,One-Coeff_C,T_X)
*                                                                      *
*---- Lee-Yang-Parr Correlation Functional                             *
*                                                                      *
      Call LYP(Rho,nRho,mGrid,
     &         dF_dRho,ndF_dRho,
     &         Coeff_C,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
C      Call QExit('B2PLYP')
C
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
