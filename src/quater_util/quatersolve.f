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
      subroutine quatersolve(U1,U2,V1,V2,Q)
************************************************************
*
*   <DOC>
*     <Name>quatersolve</Name>
*     <Syntax>quatersolve(U1,U2,V1,V2,Q)</Syntax>
*     <Arguments>
*       \Argument{U1}{vector U1, Dimension(3)}{Real*8}{in-out}
*       \Argument{U2}{vector U2, Dimension(3)}{Real*8}{in-out}
*       \Argument{V1}{vector V1, Dimension(3)}{Real*8}{in-out}
*       \Argument{V2}{vector V2, Dimension(3)}{Real*8}{in-out}
*       \Argument{Q}{quaternion Q, Dimension(4)}{Real*8}{out}
*     </Arguments>
*     <Purpose>Computes the quaternion which corresponds to the
*          best rotation that transforms U1 in U2 and V1 in V2.
*          </Purpose>
*     <Dependencies>quater util blas and util</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Input vectors will be normalized and V2 changed.
*        Computes the quaternion which corresponds to the
*          rotation that transforms U1 in U2 and V1 in V2.
*          If such a rotation does not exists i.e. if
*          U1.V1 .NE. U2.V2
*          the quaternion contains the rotation
*          which transforms U1 into U2 and V1 into V2' such that
*             U1.V1 == U2.V2'
*          and
*             U1xV1=b * (U2xV2') with b real
*     </Description>
*    </DOC>
*
************************************************************

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
            call dcopy_(4,Zero,0,Q,1)
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
