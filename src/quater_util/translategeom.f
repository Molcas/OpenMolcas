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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
*  translategeom
*
*> @brief
*>   Performs the translation of geometry 2 to geometry 1
*> @author Y. Carissan
*>
*> @details
*> Performs the translation of geometry 2 to
*> geometry 1 and stores the result in geometry 3.
************************************************************************
      subroutine translategeom(Vtrans, nat, Geom1, Geom2)
      implicit none
#include "real.fh"
      Integer nat,iat
      Real*8 Vtrans(3)
      Real*8 Geom1(3,nat)
      Real*8 Geom2(3,nat)

      Do iat=1,nat
        call dcopy_(3,Geom1(1,iat),1,Geom2(1,iat),1)
        Call Daxpy_(3,One,Vtrans,1,Geom2(1,iat),1)
      End do
      Return
      End
