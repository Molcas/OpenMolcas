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
* Copyright (C) 2017, Giovanni Li Manni                                *
*               2017, Aron Cohen                                       *
************************************************************************
      Subroutine S12h(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object:     Combination of exchange S12h (25% HF exchange) and PBE   *
*             correlation terms                                        *
*             Reference:  Swart, Marcel                                *
*             Chemical Physics Letters 580 (2013) 166-171              *
*                                                                      *
*      Author: G. Li Manni and A. Cohen, Department of Electronic      *
*              Structure Theory, Max Planck Institute, Stuttgart       *
*              2017                                                    *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange
*
C      Coeff=1.079966d0
C
C      Call Diracx(mGrid,iSpin,F_xc,Coeff)
*
*---- S12h has its LDA part included!
*
!     Coeff=0.75d0
      Coeff=One
      Call xS12h(mGrid,Coeff,iSpin,F_xc)
*
*---- PBE Correlation
*
      Coeff=1.0d0
      Call CPBE(mGrid,Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
