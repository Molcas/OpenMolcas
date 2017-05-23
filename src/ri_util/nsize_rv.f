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
      Integer Function nSize_Rv(kS,lS,nShBf,nShell,nIrrep,iOff,
     &                          nVec)
************************************************************************
*                                                                      *
*     Compute the size of Rv(nu,mu,K) and the offsets to the           *
*     different symmetry blocks.                                       *
*                                                                      *
************************************************************************
      Integer nShBf(0:nIrrep-1,nShell), iOff(0:nIrrep-1),
     &        nVec(0:nIrrep-1)
*
      nSize_Rv=0
*
      If (nIrrep.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
         iOff(0)=0
         If (kS.ne.lS) Then
            nK = nShBf(0,kS)
            nL = nShBf(0,lS)
            nKL = nK*nL
         Else
            nK = nShBf(0,kS)
            nKL = nK*(nK+1)/2
         End If
*
         nJ = nVec(0)
         nSize_Rv = nJ*nkl
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         Call IZero(iOff,nIrrep)
         Do klIrrep = 0, nIrrep-1
            iOff(klIrrep) = nSize_Rv
*
            nKL = 0
            If (kS.ne.lS) Then
               Do kIrrep = 0, nIrrep-1
                  nK = nShBf(kIrrep,kS)
                  lIrrep = iEor(klIrrep,kIrrep)
                  nL = nShBf(lIrrep,lS)
                  nKL = nKL + nK*nL
               End Do
            Else
               Do kIrrep = 0, nIrrep-1
                  nK = nShBf(kIrrep,kS)
                  lIrrep = iEor(klIrrep,kIrrep)
                  nL = nShBf(lIrrep,lS)
*
                  If (kIrrep.gt.lIrrep) Then
                     nKL = nKL + nK*nL
                  Else If (kIrrep.eq.lIrrep) Then
                     nKL = nKL + nK*(nK+1)/2
                  Else
                     nKL = nKL + 0
                  End If
*
               End Do
            End If
*
            nJ = nVec(klIrrep)
            nSize_Rv = nSize_Rv + nJ*nKL
*
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
