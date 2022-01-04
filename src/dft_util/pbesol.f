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
      Subroutine PBEsol(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object: To compute the sum of the functional c_pbesol and x_pbesol   *
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
      CoeffA=One*CoefR
      Call CPBEsol(mGrid,CoeffA,iSpin,F_xc)

      CoeffB=One*CoefX
      Call XPBEsol(mGrid,CoeffB,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
