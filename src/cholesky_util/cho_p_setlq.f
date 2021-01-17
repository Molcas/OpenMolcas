************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_P_SetLQ()
C
C     Purpose: set local qualified indices from known global qualified.
C
C     The local indices are stored in cholq.fh:
C
C     nQual_L(iSym)   : #qualified, irrep iSym (=1,2,..,nSym)
C     iQuAB_L(iQ,iSym): address of qualified iQ of sym. iSym in current
C                       local reduced set (i.e. reduced set at location
C                       2). Not symmetry reduced, i.e. includes the
C                       offset iiBstR(iSym,2).
C     iQL2G(iQ,iSym)  : iWork(ip_iQL2G-1+MaxQual*(iSym-1)+iQ)
C                       returns index of the qualified in the global
C                       list.
C
      use ChoSwp, only: iQuAB, iQuAB_L
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "cholq.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"

      Integer iSym, nQL, kOff2, iQ, iQG, i2, i, j, k

      Integer IndRed_G, IndRed, iL2G
      IndRed_G(i,j)=iWork(ip_IndRed_G-1+mmBstRT_G*(j-1)+i)
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      iL2G(i)=iWork(ip_iL2G-1+i)

      If (.not.Cho_Real_Par) Return ! not truely parallel...

      Call Cho_iZero(iQuAB_L,SIZE(iQuAB_L))
      Call Cho_iZero(iWork(ip_iQL2G),l_iQL2G)
      Do iSym = 1,nSym
         nQL = 0
         kOff2 = ip_iQL2G   + MaxQual*(iSym-1) - 1
         Do iQ = 1,nQual(iSym)
            iQG = IndRed_G(iQuAB(iQ,iSym),2) ! addr of qual in glob. rs1
            i2 = iiBstR(iSym,2) + nnBstR(iSym,2)
            i = iiBstR(iSym,2)
            Do While (i .lt. i2)
               i = i + 1
               j = IndRed(i,2) ! addr in local rs1
               k = iL2G(j)     ! addr in global rs1
               If (k .eq. iQG) Then ! found qual in local set
                  nQL = nQL + 1
                  iQuAB_L(nQL,iSym) = i
                  iWork(kOff2+nQL) = iQ
                  i = i2 ! break while loop
               End If
            End Do
         End Do
         nQual_L(iSym) = nQL
      End Do

      End
