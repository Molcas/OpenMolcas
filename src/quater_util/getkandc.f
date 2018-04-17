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
*  GetKandC
*
*> @brief
*>   Utility routine for quaternion resolution
*> @author Y. Carissan
*>
*> @details
*> Performs the following operation:
*>
*> \f[ K = (V_1-U_1) \times (V_2-U_2) \\
*>     C = K . (U_1 \times U_2) \f]
*>
*> @param[in]  U1 Input vector \f$ U_1 \f$
*> @param[in]  U2 Input vector \f$ U_2 \f$
*> @param[in]  V1 Input vector \f$ V_1 \f$
*> @param[in]  V2 Input vector \f$ V_2 \f$
*> @param[out] K  Output vector \f$ K \f$
*> @param[out] C  Output value \f$ C \f$
************************************************************************
      Subroutine GetKandC(U1,U2,V1,V2,K,C)
      Implicit none
#include "debug.fh"
#include "real.fh"
      Real*8 U1(3),V1(3),U2(3),V2(3),K(3)
      Real*8 T1(3),T2(3)
      Real*8 C
      real*8 ddot_
      external ddot_

      call dcopy_(3,V1,1,T1,1)   ! T1=V1
      call dcopy_(3,V2,1,T2,1)   ! T2=V2

      Call daxpy_(3,-One,U1,1,T1,1) ! T1 = T1-U1 = V1-U1
      Call daxpy_(3,-One,U2,1,T2,1) ! T2 = T2-U2 = V2-U2

      Call cross(T1,T2,K)          ! K = T1 x T2 = (V1-U1) x (V2-U2)
      Call cross(U1,U2,T1)         ! T1 = U1 x U2
      C=ddot_(3,K,1,T1,1)           ! C = K . (U1 x U2)

      if (debug) then
        Call RecPrt("K",' ',K,3,1)
        Write(6,*) "C",C
      end if

      Return ! K and C

      End
