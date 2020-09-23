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
* Copyright (C) 1990,1996, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine SetUp_RW(DoRys,nDiff)
************************************************************************
*                                                                      *
* Object: to setup tables for auxiliary functions to be used direct    *
*         in the recurrence relations of integral form or indirectly   *
*         to compute the Rys roots and weights which are used in the   *
*         recurrence relations of integrand form. For the lower order  *
*         Rys polynomials the roots and weight are computed from expa- *
*         nsion coefficients.                                          *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              SetHer                                                  *
*              SetUpR                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*                                                                      *
*             Unified version August '96, RL.                          *
************************************************************************
      use External_Centers, only: XF, nOrdEF
      use Sizes_of_Seward, only: S
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
      Logical DoRys
*
*     Compute max sum of angular momentum index
*
      iAng2 = 4*S%iAngMx
*
*     Set up roots and weights for Hermite polynomials.
*
      Call SetHer(nDiff)
*
*     Set up coefficients for Rys polynomials.
*
*     1) for two-electron integrals
*     2) for external field and nuclear attraction
*
      mRys =(iAng2+2+nDiff)/2
      If (Allocated(XF).or.(nOrdEF.eq.1).or.GIAO)
     &   mRys=Max(mRys,(2*S%iAngMx+1+2+nDiff)/2)
      If (nOrdEF.eq.2) mRys=Max(mRys,(2*S%iAngMx+2+2+nDiff)/2)
      If (DoRys) Call SetUpR(mRys)
*
      Return
      End
