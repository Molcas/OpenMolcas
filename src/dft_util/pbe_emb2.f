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
      use libxc, only: Only_exc
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
      Only_exc=.True.
      Call PBE(mGrid,nDmat,F_xc)
      Only_exc=.False.
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
