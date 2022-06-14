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
*  CHO_GET_OED_cd
*
*> @brief
*>   Computes the Cholesky vectors of the 2nd-rank orbital energy denominators according to Eq. (7) of \cite Koc2000-JCP-113-508
*> @author F. Aquilante
*>
*> @details
*> Computes the Cholesky vectors of the 2nd-rank orbital energy denominators
*> according to Eq. (7) of \cite Koc2000-JCP-113-508
*>
*> @param[in]     incore Boolean for in or out of core decomposition
*> @param[in]     nOV    Number of (occ,vir) pairs matching a given compound symmetry
*> @param[in]     W      Array (\p nOV) of the orbital energy differences \p W(bj) = ``e(b)-e(j)``
*> @param[in]     NVec   Number of Cholesky vectors requested
*> @param[in]     ID_bj  Index array (\p NVec) of the decomposition pattern
*> @param[in]     JVec   Starting Cholesky vector
*> @param[in,out] Y      Matrix (\p nOV,\p NVec) of the Cholesky vectors for the decomposition of the orbital energy denominators
*> @param[in,out] R      Matrix (\p nOV) of the previous products, to be used for out of core decomposition.
*>                       Note: in case of out of core, in the first call this array *must* be initialized to ``1.0d0``.
*>                       For the in core case, it *must* also coincide with the first column of \p Y and *must* still be initialized in input.
************************************************************************
      SUBROUTINE CHO_GET_OED_cd(incore,nOV,W,NVec,ID_bj,JVec,Y,R)
      Implicit Real*8 (a-h,o-z)

      Logical incore
      Integer nOV, NVec, JVec, ID_bj(*)
      Real*8  W(*), Y(*), R(*)

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
