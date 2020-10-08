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
      Subroutine O3LYP(mGrid,Rho,nRho,P2_ontop,
     &                 nP2_ontop,iSpin,F_xc,
     &                 dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
************************************************************************
*                                                                      *
* Object:  The O3LYP combination as defined in Hoe, Cohen  and Handy,  *
*            Chem.Phys.Lett 341 (2001) 319                             *
*          Note that it really is O4LYP, breaking UEG limit, and thus  *
*          is probably not the same as Gaussian's O3LYP implementation *
*          that, according to their docs, has (1-a) like Becke's B3LYP *
*                                                                      *
* Author:    Grigory A. Shamov, U. of Manitoba, 2009                   *
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
      Coeff_A=0.9262D0*CoefX
      Coeff_B=0.8133D0*CoefX
      Coeff_C=0.81D0*CoefR
*                                                                      *
*---- Dirac Exchange Functional                                        *
*                                                                      *
      Call Diracx(mGrid,Rho,nRho,
     &            iSpin,F_xc,dF_dRho,
     &            ndF_dRho,Coeff_A,T_X)
*                                                                      *
*---- OPTX Exchange Functional                                         *
*                                                                      *
      Call xOPT(Rho,nRho,mGrid,
     &          dF_dRho,ndF_dRho,
     &          Coeff_B,iSpin,F_xc,T_X)
*                                                                      *
*---- Vosko-Wilks-Nusair Correlation Functional III                    *
*                                                                      *
      Call VWN_III(mGrid,Rho,nRho,
     &             iSpin,
     &             F_xc,dF_dRho,ndF_dRho,CoefR - Coeff_C,T_X)
*                                                                      *
*---- Lee-Yang-Parr Correlation Functional                             *
*                                                                      *
      Call LYP(Rho,nRho,mGrid,
     &         dF_dRho,ndF_dRho,
     &         Coeff_C,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
