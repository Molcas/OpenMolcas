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
*  quatersolve
*
*> @brief
*>   Computes the quaternion which corresponds to the best rotation that transforms \p U1 in \p U2 and \p V1 in \p V2
*> @author Y. Carissan
*>
*> @details
*> Input vectors will be normalized and \p V2 changed.
*> Computes the quaternion which corresponds to the
*> rotation that transforms \p U1 in \p U2 and \p V1 in \p V2.
*> If such a rotation does not exists, i.e. if
*> \f[ U_1 \cdot V_1 \ne U_2 \cdot V_2 \f]
*> the quaternion contains the rotation
*> which transforms \p U1 into \p U2 and \p V1 into \p V2' such that
*> \f[ U_1 \cdot V_1 = U_2 \cdot V_2' \f]
*> and
*> \f[ U_1 \times V_1 = b (U_2 \times V_2') \f] with \f$ b \f$ real.
*>
*> @param[in,out] U1 vector \f$ U_1 \f$
*> @param[in,out] U2 vector \f$ U_2 \f$
*> @param[in,out] V1 vector \f$ V_1 \f$
*> @param[in,out] V2 vector \f$ V_2 \f$
*> @param[out]    Q  quaternion \f$ Q \f$
************************************************************************
      subroutine quatersolve(U1,U2,V1,V2,Q)
      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "real.fh"
      Real*8 U1(3),V1(3),U2(3),V2(3)
      Real*8 Q(0:3)
      Real*8 U3(3),V3(3)
      Real*8 Uref(3),Vref(3),K(3)
      Real*8 Vtmp(3)
      Real*8 c, thrs, ddot_
      external ddot_

      thrs=1d-3

      if (debug) then
        Call RecPrt("IN SOLVE U1",' ',U1,3,1)
        Call RecPrt("IN SOLVE V1",' ',V1,3,1)
        Call RecPrt("IN SOLVE U2",' ',U2,3,1)
        Call RecPrt("IN SOLVE V2",' ',V2,3,1)
      end if
      call quatersetup(U1,U2,V1,V2)
      if (debug) Call RecPrt("new V2",' ',V2,3,1)

      call dcopy_(3,U1,1,Uref,1)
      call dcopy_(3,V1,1,Vref,1)

      call GetKandC(U1,U2,V1,V2,K,C)

      if ( C.lt.thrs ) then
        Call cross(U1,U2,U3)
        Call cross(V1,V2,V3)
        Call getKandC(U1,V1,U3,V3,K,C)
        if ( C.lt.thrs ) then
          Call getKandC(U2,V2,U3,V3,K,C)
          if ( C.lt.thrs ) then
            call dcopy_(4,[Zero],0,Q,1)
            Q(0)=One
            Go To 999
          end if
          call dcopy_(3,U2,1,Uref,1)
          call dcopy_(3,V2,1,Vref,1)
        end if
      end if

      Q(1) = Half * K(1)/sqrt(C)
      Q(2) = Half * K(2)/sqrt(C)
      Q(3) = Half * K(3)/sqrt(C)

      Call Cross(Uref,Q(1),Vtmp) ! Vtmp=Uref x Q

      Q(0) = Half * ddot_(3,Vref,1,Vtmp,1) / ddot_(3,Vtmp,1,Vtmp,1)
c                     ! Q(0)=0.5 * Vref.(UrefxQ)/(UrefxQ)^2

999   continue ! normal endding
      Call CheckQuater(Q)
      Call setMatrix(Q)
      if(debug) Call RecPrt("Quaternion",' ',Q,4,1)

      return
      end
