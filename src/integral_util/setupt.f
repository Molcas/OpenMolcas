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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine SetUpT(DSCF,nDiff)
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
*              GetMem                                                  *
*              SetHer                                                  *
*              SetUpR                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Logical DSCF
*
      iRout=7
      iPrint = nPrint(iRout)
*     iQ = 1
      Call qEnter('SetUpT')
      If (iPrint.ge.99) Call GetMem(' In SetupT','CHECK','REAL',
     &                              iDum,iDum)
*
*     Compute max sum of angular momentum index
*
      iAng2 = 4*iAngMx
*
*     Set up roots and weights for Hermite polynomials.
*
      Call SetHer(nDiff)
*
*     Set up coefficients for Rys polynomials.
*
      mRys =(iAng2+2)/2
      If (DSCF) Call SetUpR(mRys)
*
*     Call GetMem('SetUpT','CHECK','REAL',iDum,iDum)
      Call qExit('SetUpT')
      Return
      End
