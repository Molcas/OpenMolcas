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
      Subroutine cBLYP_emb(mGrid,iSpin,F_xc)
************************************************************************
      use OFembed, only: dFMD
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
*
*---- Lee-Yang-Parr Correlation
*
      Coeff=dFMD
      Call LYP(mGrid,
     &         Coeff,iSpin,F_xc)
*
      Return
      End
