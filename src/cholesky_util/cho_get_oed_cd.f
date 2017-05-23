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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE CHO_GET_OED_cd(incore,nOV,W,NVec,ID_bj,JVec,Y,R)
***********************************************************************
*
*   <DOC>
*     <Name>CHO\_GET\_OED\_cd</Name>
*     <Syntax>Call CHO\_GET\_ORD\_bj(incore,nOV,W,ID\_bj,NVec,JVec,Y,R)</Syntax>
*     <Arguments>
*       \Argument{incore}{Boolean for in or out of core decomposition}{Logical}{in}
*       \Argument{nOV}{Number of (occ,vir) pairs matching a given
*                      compound symmetry}{Integer}{in}
*       \Argument{W}{Array (nOV) of the orbital energy
*                    differences W(bj) = e(b)-e(j)}{Array Real*8}{in}
*       \Argument{NVec}{Number of Cholesky vectors requested}{Integer}{in}
*       \Argument{ID\_bj}{Index array (NVec) of the decomposition pattern}{Array Integer}{in}
*       \Argument{JVec}{Starting Cholesky vector}{Integer}{in}
*       \Argument{Y}{Matrix (nOV,NVec) of the Cholesky vectors
*          for the decomposition of the orbital energy denominators.}{Array Real*8}{inout}
*       \Argument{R}{Matrix (nOV) of the previous products, to be
*          used for out of core decomposition. Note: in case of out of
*          core, in the first call this array MUST be initialized to 1.0d0. For the in core
*          case, it MUST also coincide with the first column of Y and MUST
*          still be initialized in input.}{Array Real*8}{inout}
*     </Arguments>
*     <Purpose>
*        Computes the Cholesky vectors of the 2nd-rank orbital energy
*        denominators according to Eq. (7) of H. Koch, A. M. Sanchez,
*        JCP 113, 508 (2000).
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>F. Aquilante</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     </Description>
*    </DOC>
*
************************************************************

      Implicit Real*8 (a-h,o-z)

      Logical incore
      Integer nOV, NVec, JVec, ID_bj(*)
      Real*8  W(*), Y(*), R(*)

#include "WrkSpc.fh"
#include "warnings.fh"
*******************************************************************

      If (NVec.lt.1) Then
         Write(6,*)'Error in CHO_GET_OED_cd : in input NVec < 1 .'
         Call quit(_RC_CHO_LOG_)
      EndIf

      xtwo = sqrt(2.0d0)

* Compute  R(p,k) = R(p,k-1) * (W(p) - W(J[k-1]))/(W(p) + W(J[k-1]))
* ------------------------------------------------------------------
      If (incore) Then

         If (JVec.ne.1) Then
            Write(6,*)'CHO_GET_OED_cd : JVec must be 1 if incore .'
            Call quit(_RC_CHO_LOG_)
         EndIf
         Do Jk=2,NVec
            kp = NOV*(Jk-1)
            np = kp - NOV
            Jm = ID_bj(Jk-1)
            Do ip=1,NOV
               Y(kp+ip) = Y(np+ip)*(W(ip)-W(Jm))/(W(ip)+W(Jm))
            End Do
          End Do

      Else

         Jm = ID_bj(JVec-1)
         Do ip=1,NOV
            Y(ip) = R(ip)*(W(ip)-W(Jm))/(W(ip)+W(Jm))
         End Do
         Do Jk=2,NVec
            kp = NOV*(Jk-1)
            np = kp - NOV
            Jm = ID_bj(JVec+Jk-3)
            Do ip=1,NOV
               Y(kp+ip) = Y(np+ip)*(W(ip)-W(Jm))/(W(ip)+W(Jm))
            End Do
          End Do
          call dcopy_(NOV,Y(1+NOV*(NVec-1)),1,R,1)

      EndIf

* Compute  Y(p,k) = R(p,k) * sqrt(2*W(J[k]))/(W(p) + W(J[k]))
* -----------------------------------------------------------
      Do Jk=1,NVec
         kp = NOV*(Jk-1)
         Jm = ID_bj(JVec+Jk-1)
         Do ip=1,NOV
            Y(kp+ip) = Y(kp+ip)*xtwo*sqrt(W(Jm))/(W(ip)+W(Jm))
         End Do
      End Do

      Return
      End
