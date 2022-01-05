************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine KT2(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object:     KT2 combination, as described in Keal&Tozer paper        *
*               J.Chem.Phys 119 (2003) 3015                            *
*                                                                      *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange with non-UEG factor
*
      Coeff= 1.07173d0
      Call Diracx(mGrid,iSpin,F_xc,Coeff)
*
*---- KT term Exchange/Correlation
*
      Coeff= - 0.0060d0
      Call KealTozer(mGrid,
     &          Coeff,iSpin,F_xc)
*
*---- VWN-3 Correlation
*
      Coeff=0.576727d0
      Call VWN_III(mGrid,iSpin,F_xc,Coeff)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
