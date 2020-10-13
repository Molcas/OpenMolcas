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
      Subroutine O2PLYP(mGrid,Rho,nRho,P2_ontop,
     &                  nP2_ontop,iSpin,F_xc,
     &                  dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
************************************************************************
*                                                                      *
* Object:   OPTX analog of Grimme's B2PLYP double hybrid               *
*           preserves non-UEG feature of the OPTX combinations         *
*           requires scaling of &mbpt2 correlation energy!             *
*                                                                      *
*      Author: Grigory A Shamov, U of Manitoba, 2009                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
C "physically motivated" values from Grimme's B2PLYP paper are 0.5 * GGA XC
* + 0.5* HF; since the parent GGA OLYP is non-UEG Grigory took 0.5* OLYP's
* LDA and GGA coefficients here
      Coeff_A=0.525755D0*CoefX
      Coeff_B=0.715845D0*CoefX
      Coeff_C=0.75D0*CoefR
C MP2 correlation energy to be scaled to 1 - Coeff_C = 0.25

*                                                                      *
*---- Dirac Exchange Functional     * 0.5 OPTX                         *
*                                                                      *
      Call Diracx(mGrid,Rho,nRho,
     &            iSpin,F_xc,dF_dRho,
     &            ndF_dRho,Coeff_A,T_X)
*                                                                      *
*---- OPTX Exchange Functional * 0.5  OPTX                             *
*                                                                      *
      Call xOPT(Rho,nRho,mGrid,
     &          dF_dRho,ndF_dRho,
     &          Coeff_B,iSpin,F_xc,T_X)
*                                                                      *
*  LYP correlation                                                     *
*                                                                      *
      Call LYP(Rho,nRho,mGrid,
     &         dF_dRho,ndF_dRho,
     &         Coeff_C,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
C
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
