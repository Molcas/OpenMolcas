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
       Subroutine HFO(mGrid,iSpin,F_xc)
************************************************************************
*                                                                      *
* Object:            standalone OPTX exchange                          *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. March 2001                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac (Slater) exchange                                          *
*                                                                      *
      Coeff=1.051510d0*CoefX
      Call DiracX(mGrid,iSpin,F_xc,Coeff)
*                                                                      *
*---- OPTX Exchange Functional                                         *
*                                                                      *
      Coeff=1.431690d0*CoefX
      Call xOPT(mGrid,
     &          Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
