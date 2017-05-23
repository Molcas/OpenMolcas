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
      Subroutine Stblz(iChxyz,iOper,nIrrep,nStab,jStab,MaxDCR,iCoSet)
      Implicit Real*8 (a-h,o-z)
      Integer iOper(0:nIrrep-1), jStab(0:7), iCoSet(0:7,0:7)
*                                                                      *
************************************************************************
*                                                                      *
*     Find the Stabilizers of this center
*                                                                      *
************************************************************************
*                                                                      *
      iStab = 0
      Do i = 0, nIrrep-1
         If (iAnd(iChxyz,iOper(i)).eq.0) Then
            iStab = iStab + 1
            jStab(iStab-1) = iOper(i)
         End If
      End Do
      nStab = iStab
      MaxDCR = Max(MaxDCR,iStab)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate all possible (left) CoSet
*                                                                      *
************************************************************************
*                                                                      *
      Do i = 0, nIrrep-1
         Do j = 0, iStab-1
            iCoSet(i,j) = iEor(iOper(i),jStab(j))
         End Do
      End Do
      If (iStab.eq.1) Go To 39 ! skip if all are unique
*
*     Order the Coset so we will have the unique ones first
*
      nMax = 1
      If (nMax.eq.nIrrep/iStab) Go To 39 ! skip if there is only one
      iTest=iStab-1 ! Test on the last element
      Do 35 j = 1, nIrrep-1 !
*        Check uniqueness
         Do i = 0, nMax-1
            Do ielem = 0, iStab-1
               If (iCoSet(i,iTest).eq.iCoSet(j,ielem)) Go To 35
            End Do
         End Do
*        Move unique CoSet to nMax+1
         nMax = nMax + 1
         Do ielem = 0, iStab-1
            iTmp = iCoSet(nMax-1,ielem)
            iCoSet(nMax-1,ielem) = iCoSet(j,ielem)
            iCoSet(j,ielem) = iTmp
         End Do
         If (nMax.eq.nIrrep/iStab) Go To 39
 35   Continue
 39   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do i = 0, nIrrep/iStab-1
         iOpMn=iCoSet(i,0)
         Do j = 1, iStab-1
            iOpMn=iAnd(iOpMn,iCoset(i,j))
         End Do
         ip=0
         Do j = 0, iStab-1
            If (iOpMn.eq.iCoSet(i,j)) ip = j
         End Do
         iTmp=iCoSet(i,0)
         iCoSet(i,0)=iCoSet(i,ip)
         iCoSet(i,ip)=iTmp
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
