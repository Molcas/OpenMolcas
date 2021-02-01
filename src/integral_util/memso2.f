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
      Integer Function MemSO2(iAng,jAng,kAng,lAng,
     &                        iCmp,jCmp,kCmp,lCmp,
     &                        iShell,jShell,kShell,lShell,
     &                        iAO,jAO,kAO,lAO)
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
      Logical Shij, Shkl, Shik, Shjl
*
      MemSO2 = 0
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      Shij = iShell.eq.jShell
      Shkl = kShell.eq.lShell
      Shik = iShell.eq.kShell
      Shjl = jShell.eq.lShell
*
      If (nIrrep.eq.1) Then
*
         Do i1 = 1, iCmp
            jCmpMx = jCmp
            If (Shij) jCmpMx = i1
            Do i2 = 1, jCmpMx
               kCmpMx = kCmp
               If (Shik .and. Shjl) kCmpMx = i1
               Do i3 = 1, kCmpMx
                  lCmpMx = lCmp
                  If (Shkl) lCmpMx = i3
                  If (Shik .and. i1.eq.i3 .and. Shjl) lCmpMx = i2
                  MemSO2 = MemSO2 + lCmpMx
               End Do
            End Do
         End Do
*
      Else
*
         Do i1 = 1, iCmp
            jCmpMx = jCmp
            If (Shij) jCmpMx = i1
            Do i2 = 1, jCmpMx
               kCmpMx = kCmp
               If (Shik .and. Shjl) kCmpMx = i1
               Do i3 = 1, kCmpMx
                  lCmpMx = lCmp
                  If (Shkl) lCmpMx = i3
                  If (Shik .and. i1.eq.i3 .and. Shjl) lCmpMx = i2
                  Do i4 = 1, lCmpMx
*
*         Loop over irreps which are spanned by the basis function.
*         Again, the loop structure is restricted to ensure unique
*         integrals.
*
          Do 110 j1 = 0, nIrrep-1
             If (iAOtSO(iAO+i1,j1)<0) Cycle
             j2Max = nIrrep-1
             If (Shij .and. i1.eq.i2) j2Max = j1
             Do 210 j2 = 0, j2Max
                If (iAOtSO(jAO+i2,j2)<0) Cycle
                j12 = iEor(j1,j2)
                j3Max = nIrrep-1
                If (Shik .and. i1.eq.i3 .and.
     &              Shjl .and. i2.eq.i4) j3Max = j1
                Do 310 j3 = 0, j3Max
                   If (iAOtSO(kAO+i3,j3)<0) Cycle
                   j4 = iEor(j12,j3)
                   If (iAOtSO(lAO+i4,j4)<0) Cycle
                   If (Shkl .and. i3.eq.i4 .and. j4.gt.j3) Go To 310
                   If (Shik .and. i1.eq.i3 .and. Shjl .and. i2.eq.i4
     &                      .and. j1.eq.j3 .and. j4.gt.j2) Go To 310
                   MemSO2 = MemSO2 + 1
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
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iAng)
         Call Unused_integer(jAng)
         Call Unused_integer(kAng)
         Call Unused_integer(lAng)
      End If
      End
