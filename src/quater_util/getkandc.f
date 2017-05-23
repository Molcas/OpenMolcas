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
      Subroutine GetKandC(U1,U2,V1,V2,K,C)
************************************************************
*
*   <DOC>
*     <Name>GetKandC</Name>
*     <Syntax>GetKandC(U1,U2,V1,V2,K,C)</Syntax>
*     <Arguments>
*       \Argument{U1}{Input Vector U1, Dimension(3)}{Real*8}{in}
*       \Argument{U2}{Input Vector U2, Dimension(3)}{Real*8}{in}
*       \Argument{V1}{Input Vector V1, Dimension(3)}{Real*8}{in}
*       \Argument{V2}{Input Vector V2, Dimension(3)}{Real*8}{in}
*       \Argument{K}{Output Vector K, Dimension(3)}{Real*8}{out}
*       \Argument{C}{Output value C}{Real*8}{out}
*     </Arguments>
*     <Purpose>Utility routine for quaternion resolution</Purpose>
*     <Dependencies>blas routines and call cross</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Performs the following operation :
*  K = (V1-U1) x (V2-U2)
*  C = K . (U1 x U2)
*     </Description>
*    </DOC>
*
************************************************************
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
