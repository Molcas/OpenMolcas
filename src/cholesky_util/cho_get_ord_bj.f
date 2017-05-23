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
*               2002, Thomas Bondo Pedersen                            *
************************************************************************
      SUBROUTINE CHO_GET_ORD_bj(nOV,MaxNVec,thr,W,ID_bj,nVec,Dmax)
************************************************************************
*
*   <DOC>
*     <Name>CHO\_GET\_ORD\_bj</Name>
*     <Syntax>Call CHO\_GET\_ORD\_bj(nOV,MaxNVec,thr,W,ID\_bj,nVec,Dmax)</Syntax>
*     <Arguments>
*       \Argument{nOV}{Number of (occ,vir) pairs matching a given
*                      compound symmetry}{Integer}{in}
*       \Argument{MaxNVec}{Max number of Cholesky vectors}{Integer}{in}
*       \Argument{thr}{Threshold for the Cholesky decomposition}{Real*8}{in}
*       \Argument{W}{Array (nOV) of the orbital energy
*                    differences W(bj) = e(b)-e(j)}{Array Real*8}{in}
*       \Argument{ID\_bj}{Index array (MaxNVec) of the decomposition pattern}{Array Integer}{out}
*       \Argument{nVec}{Number of resulting Cholesky vectors}{Integer}{out}
*       \Argument{Dmax}{Max value of the remainder after the
*                       decomposition}{Real*8}{out}
*     </Arguments>
*     <Purpose>
*        Computes the decomposition pattern of the 2nd-rank orbital energy
*        denominators according to Eq. (6) of H. Koch, A. M. Sanchez,
*        JCP 113, 508 (2000).
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>F. Aquilante</Author>
*     <Modified_by>
*     Thomas Bondo Pedersen, December 2012: code restructured, fix of
*     minor bug.
*     </Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*            The decomposition is considered converged when
*            either of the two criteria (MaxNVec or thr) is fulfilled.
*            Therefore these two arguments should have appropriate
*            values depending on what the user expects:
*
*            ex.    MaxNVec := nOV  if the user specifies
*                                      a meaningful thr
*                   thr := zero  if the user requires a given number
*                                of vectors (or percentage of nOV)
*     </Description>
*    </DOC>
*
************************************************************

      Implicit Real*8 (a-h,o-z)

      Integer nOV, MaxNVec, ID_bj(*), NVec
      Real*8  thr, W(*), Dmax

      Parameter (zero = 0.0d0, half = 0.5d0)
#include "WrkSpc.fh"
*******************************************************************
      D(i) = Work(ipD+i)
*******************************************************************

* Initialize
* ----------
      nVec=0
      If (nOV.lt.1) Then
         Dmax=-9.987654321d0 ! dummy value
         Return
      End If

* Allocate diagonal
* -----------------
      Call GetMem('Diag','Allo','Real',ipDiag,NOV)
      ipD=ipDiag-1

* Compute diagonals
* -----------------
      Do ip=1,nOV
         If (W(ip).gt.zero) Then
            Work(ipD+ip)=half/W(ip)
         Else ! tbp: perhaps we should stop it here (matrix not PSD)
            Work(ipD+ip)=zero
         End If
      End Do

* Find ID (Jm) of max diagonal
* ----------------------------
      Jm=1
      Do ip=2,NOV
         If (D(ip).gt.D(Jm)) Then
            Jm=ip
         End If
      End Do

* Find CD pattern using the diagonal update:
*   D(p)[k] = D(p)[k-1] * ((W(p) - W(J[k-1]))/(W(p) + W(J[k-1])))^2
* -----------------------------------------------------------------
      Do While (nVec.lt.MaxNVec .and. D(Jm).gt.thr)
        ! update vector counter
        nVec=nVec+1
        ! save ID of vector = ID of current max diagonal
        ID_bj(nVec)=Jm
        ! update diagonals
        Do ip=1,nOV
           Work(ipD+ip)=D(ip)*((W(ip)-W(Jm))/(W(ip)+W(Jm)))**2
        End Do
        ! find ID (Jm) of max updated diagonal
        Jm=1
        Do ip=2,nOV
           If (D(ip).gt.D(Jm)) Then
              Jm=ip
           End If
        End Do
      End Do

* Set max diagonal (Dmax) to return to caller
* -------------------------------------------
      Dmax=D(Jm)

* Deallocate diagonal
* -------------------
      Call GetMem('Diag','Free','Real',ipDiag,NOV)

      Return
      End
