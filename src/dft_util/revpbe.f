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
* Copyright (C) 2005, Per Ake Malmqvist                                *
************************************************************************
      Subroutine REVPBE(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object: To compute the sum of the functional c_pbe and x_pbe as      *
* described in the Density Functional Repository,                      *
* (http://www.cse.clrc.ac.uk/qcg/dft) and in the article               *
*    J.P. Perdew, K. Burke, and M. Ernzerhof,                          *
*  "Generalized gradient approximation made simple,"                   *
*      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. December 2005                *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *

      CoeffB=One*CoefX
      Call XrevPBE(mGrid,CoeffB,iSpin,F_xc)

      CoeffA=One*CoefR
      Call CPBE(mGrid,CoeffA,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
