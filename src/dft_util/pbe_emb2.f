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
      Subroutine PBE_emb2(mGrid,nDmat)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      use OFembed, only: KEonly
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "hflda.fh"
*
************************************************************************
*
*---- Thomas-Fermi Kinetic energy functional
*---- NDSD potential
*
      Coeff=One
      Call ndsd_Ts(mGrid,nDmat,F_xc,Coeff)

      If (KEonly) Return
*
*---- PBE for exchange-correlation energy functional (no potential)
*
      Call PBE_(mGrid,nDmat,F_xc)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *

      Subroutine PBE_(mGrid,iSpin,F_xc)
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      CoeffA=0.0D0
      Call CPBE_ofe(mGrid,
     &              CoeffA,iSpin,F_xc)

      CoeffB=0.0D0
      Call XPBE_ofe(mGrid,
     &              CoeffB,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
