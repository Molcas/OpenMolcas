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
      Integer Function nSize_3C(kS,lS,nShBf,nShell, nIrrep,iOff,
     &                          nBas_Aux)
************************************************************************
*                                                                      *
*     Compute the size of ({nu,mu}|K) and the offsets to the           *
*     different symmetry blocks.                                       *
*                                                                      *
************************************************************************
      Integer nShBf(0:nIrrep-1,nShell), iOff(3,0:nIrrep-1),
     &        nBas_Aux(0:nIrrep-1)
*
      nSize_3C=0
*
      If (nIrrep.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
         Call IZero(iOff,3)
         If (kS.ne.lS) Then
            nK = nShBf(0,kS)
            nL = nShBf(0,lS)
            nKL = nK*nL
         Else
            nK = nShBf(0,kS)
            nKL = nK*(nK+1)/2
         End If
         iOff(1,0)=nKL
*
         nJ = nBas_Aux(0)-1
         iOff(2,0) = nJ
         nSize_3C = nJ*nkl
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         Call IZero(iOff,3*nIrrep)
         Do klIrrep = 0, nIrrep-1
            iOff(3,klIrrep) = nSize_3C
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
            iOff(1,klIrrep)=nKL
*
            nJ = nBas_Aux(klIrrep)
            If (klIrrep.eq.0) nJ = nJ-1
            iOff(2,klIrrep)=nJ
            nSize_3C = nSize_3C + nJ*nKL
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
