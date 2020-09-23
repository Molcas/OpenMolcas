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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Integer Function MemSO2_P(iCmp,jCmp,kCmp,lCmp,iAO,jAO,kAO,lAO)
************************************************************************
*  Object: to compile the number of SO block which will be generated   *
*          by the current shell quadruplet.                            *
*                                                                      *
*          Observe that the indices are canonically ordered at the     *
*          time of calling this routine!                               *
*                                                                      *
* Called from: Drv2El                                                  *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             February '90                                             *
************************************************************************
      use SOAO_Info, only: iAOtSO
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
*
      MemSO2_P = 0
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      If (nIrrep.eq.1) Then
*
         MemSO2_P=iCmp*jCmp*kCmp*lCmp
*
      Else
*
         Do i1 = 1, iCmp
            Do i2 = 1, jCmp
               Do i3 = 1, kCmp
                  Do i4 = 1, lCmp
*
*         Loop over irreps which are spanned by the basis function.
*         Again, the loop structure is restricted to ensure unique
*         integrals.
*
          Do 110 j1 = 0, nIrrep-1
             If (iAOtSO(iAO+i1,j1)<0) Cycle
             Do 210 j2 = 0, nIrrep-1
                If (iAOtSO(jAO+i2,j2)<0) Cycle
                j12 = iEor(j1,j2)
                Do 310 j3 = 0, nIrrep-1
                   If (iAOtSO(kAO+i3,j3)<0) Cycle
                   j4 = iEor(j12,j3)
                   If (iAOtSO(lAO+i4,j4)<0) Cycle
                   MemSO2_P = MemSO2_P + 1
*
 310            Continue
 210         Continue
 110      Continue
*
                  End Do
               End Do
            End Do
         End Do
*
      End If
*
      Return
      End
