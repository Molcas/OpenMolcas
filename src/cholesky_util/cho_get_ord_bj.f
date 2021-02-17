
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
*  CHO_GET_ORD_bj
*
*> @brief
*>   Computes the decomposition pattern of the 2nd-rank orbital energy denominators according to Eq. (6) of \cite Koc2000-JCP-113-508
*> @author F. Aquilante
*> @modified_by Thomas Bondo Pedersen, December 2012: code restructured, fix of minor bug
*>
*> @details
*> Computes the decomposition pattern of the 2nd-rank orbital energy denominators
*> according to Eq. (6) of \cite Koc2000-JCP-113-508
*>
*> The decomposition is considered converged when
*> either of the two criteria (\p MaxNVec or \p thr) is fulfilled.
*> Therefore these two arguments should have appropriate
*> values depending on what the user expects:
*>
*> - \p MaxNVec = \p nOV if the user specifies a meaningful \p thr
*> - \p thr = ``0.0d0``  if the user requires a given number of vectors (or percentage of \p nOV)
*>
*> @param[in]  nOV     Number of (occ,vir) pairs matching a given compound symmetry
*> @param[in]  MaxNVec Max number of Cholesky vectors
*> @param[in]  thr     Threshold for the Cholesky decomposition
*> @param[in]  W       Array (\p nOV) of the orbital energy differences \p W(bj) = ``e(b)-e(j)``
*> @param[out] ID_bj   Index array (\p MaxNVec) of the decomposition pattern
*> @param[out] nVec    Number of resulting Cholesky vectors
*> @param[out] Dmax    Max value of the remainder after the decomposition
************************************************************************
      SUBROUTINE CHO_GET_ORD_bj(nOV,MaxNVec,thr,W,ID_bj,nVec,Dmax)
      Implicit Real*8 (a-h,o-z)

      Integer nOV, MaxNVec, ID_bj(*), NVec
      Real*8  thr, W(*), Dmax

#include "real.fh"
#include "stdalloc.fh"

      Real*8, Allocatable:: Diag(:)

* Initialize
* ----------
      nVec=0
      If (nOV.lt.1) Then
         Dmax=-9.987654321d0 ! dummy value
         Return
      End If

* Allocate diagonal
* -----------------
      Call mma_allocate(Diag,NOV,Label='Diag')

* Compute diagonals
* -----------------
      Do ip=1,nOV
         If (W(ip).gt.zero) Then
            Diag(ip)=half/W(ip)
         Else ! tbp: perhaps we should stop it here (matrix not PSD)
            Diag(ip)=zero
         End If
      End Do

* Find ID (Jm) of max diagonal
* ----------------------------
      Jm=1
      Do ip=2,NOV
         If (Diag(ip).gt.Diag(Jm)) Then
            Jm=ip
         End If
      End Do

* Find CD pattern using the diagonal update:
*   Diag(p)[k] = Diag(p)[k-1]
*              * ((W(p) - W(J[k-1]))/(W(p) + W(J[k-1])))^2
* -----------------------------------------------------------------
      Do While (nVec.lt.MaxNVec .and. Diag(Jm).gt.thr)
        ! update vector counter
        nVec=nVec+1
        ! save ID of vector = ID of current max diagonal
        ID_bj(nVec)=Jm
        ! update diagonals
        Do ip=1,nOV
           Diag(ip)=Diag(ip)*((W(ip)-W(Jm))/(W(ip)+W(Jm)))**2
        End Do
        ! find ID (Jm) of max updated diagonal
        Jm=1
        Do ip=2,nOV
           If (Diag(ip).gt.Diag(Jm)) Then
              Jm=ip
           End If
        End Do
      End Do

* Set max diagonal (Dmax) to return to caller
* -------------------------------------------
      Dmax=Diag(Jm)

* Deallocate diagonal
* -------------------
      Call mma_deallocate(Diag)

      Return
      End
