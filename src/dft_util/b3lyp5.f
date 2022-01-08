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
************************************************************************
      Subroutine B3LYP5(mGrid,iSpin)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. March 2001                              *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Coeff_A=0.80D0*CoefX
      Coeff_A=0.08D0*CoefX
      Coeff_B=0.72D0*CoefX
      Coeff_C=0.81D0*CoefR
*                                                                      *
*---- Dirac Exchange Functional                                        *
*                                                                      *
      Call Diracx(mGrid,
     &            iSpin,F_xc,Coeff_A)
*                                                                      *
*---- Becke 88 Exchange Functional                                     *
*                                                                      *
      Call xB88(mGrid,
     &          Coeff_B,iSpin,F_xc)
*                                                                      *
*---- Vosko-Wilk-Nusair Correlation Functional V                       *
*                                                                      *
      Call VWN_V(mGrid,iSpin,F_xc,CoefR-Coeff_C)
*                                                                      *
*---- Lee-Yang-Parr Correlation Functional                             *
*                                                                      *
      Call LYP(mGrid,
     &         Coeff_C,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
