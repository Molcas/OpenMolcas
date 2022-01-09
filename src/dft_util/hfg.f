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
* Copyright (C) 2001, Roland Lindh                                     *
************************************************************************
       Subroutine HFG(mGrid,iSpin)
************************************************************************
*                                                                      *
* Object:   Stand-alone Gill96 exchange                                *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. March 2001                              *
************************************************************************
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac (Slater) exchange                                          *
*                                                                      *
!     Coeff=One*CoefX
!     Call DiracX(mGrid,iSpin,F_xc,Coeff)
*                                                                      *
*---- Gill 96 Exchange Functional                                      *
*                                                                      *
      Coeff=One*CoefX
      Call xG96(mGrid,
     &          Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
