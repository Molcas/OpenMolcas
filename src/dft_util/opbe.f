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
      Subroutine OPBE(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object:     OPTX + PBE combination                                   *
*                 OPBE                                                 *
*                                                                      *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange with the a1 OPTX factor!
*
      Coeff= 1.051510d0*CoefX
      Call Diracx(mGrid,iSpin,F_xc,Coeff)
*
*---- OPTX Exchange, the a2 coeff is here!
*
      Coeff= 1.431690d0*CoefX
      Call xOPT(mGrid,
     &          Coeff,iSpin,F_xc)
*
*---- PBEc
*
      Coeff=One*CoefR
      Call CPBE(mGrid,
     &         Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
