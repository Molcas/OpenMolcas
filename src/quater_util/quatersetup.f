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
*  QuaterSetup
*
*> @brief
*>   Transforms \p V2 such that \f$ U_1 \cdot V_1 = U_2 \cdot V_2 \f$
*> @author Y. Carissan
*>
*> @details
*> Transforms \p V2 such that \f$ U_1 \cdot V_1 = U_2 \cdot V_2 \f$.
*>
*> @param[in]     U1 Input Vector
*> @param[in]     V1 Input Vector
*> @param[in]     U2 Input Vector
*> @param[in,out] V2 Input and output Vector
************************************************************************
      Subroutine QuaterSetup(U1,U2,V1,V2)
      Implicit none
#include "debug.fh"
#include "real.fh"
      Real*8 U1(3),V1(3)
      Real*8 U2(3),V2(3)
      Real*8 coeff,U1dU2,V1dV2
      real*8 ddot_
      external ddot_
*
      call normalizeVec(U1)
      call normalizeVec(V1)
      call normalizeVec(U2)
      call normalizeVec(V2)
*
      if (debug) then
        call RecPrt("IN QUATERSETUP normalized U1","",U1,3,1)
        call RecPrt("IN QUATERSETUP normalized V1","",V1,3,1)
        call RecPrt("IN QUATERSETUP normalized U2","",U2,3,1)
        call RecPrt("IN QUATERSETUP normalized V2","",V2,3,1)
      end if
*
      U1dU2=ddot_(3,U1,1,U2,1)
      V1dV2=ddot_(3,V1,1,V2,1)
      coeff = ( One-U1dU2**2 )/( One-V1dV2**2 )
      coeff = sqrt(coeff)
      V2(1) = ( U1dU2 - V1dV2 * coeff ) * V1(1) + coeff * V2(1)
      V2(2) = ( U1dU2 - V1dV2 * coeff ) * V1(2) + coeff * V2(2)
      V2(3) = ( U1dU2 - V1dV2 * coeff ) * V1(3) + coeff * V2(3)
      if (debug) call RecPrt("IN QUATERSETUP modified V2","",V2,3,1)
      Return
      End
