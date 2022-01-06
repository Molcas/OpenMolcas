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
      Subroutine VWN_III_emb(mGrid,iSpin)
************************************************************************
      use nq_Grid, only: F_xc => Exc
      use OFembed, only: dFMD
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "hflda.fh"
*
*---- Vosko-Wilk-Nusair correlation functional III
*
      Coeff=dFMD
      Call VWN_III(mGrid,iSpin,F_xc,Coeff)

      Return
      End
