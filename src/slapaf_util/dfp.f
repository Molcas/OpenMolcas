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
      Subroutine DFP(B,nDim,Bd,Delta,Gamma)
      Implicit Real*8 (a-h,o-z)
      Real*8 B(nDim,nDim), Bd(nDim),Gamma(nDim),Delta(nDim)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: B',' ',B,nDim,nDim)
      Call RecPrt('DFP: Bd',' ',Bd,1,nDim)
      Call RecPrt('DFP: Gamma',' ',Gamma,1,nDim)
      Call RecPrt('DFP: Delta',' ',Delta,1,nDim)
#endif
      Call DGEMM_('N','N',
     &            nDim,1,nDim,
     &            1.0d0,B,nDim,
     &                  Delta,nDim,
     &            0.0d0,Bd,nDim)
      dd=DDot_(nDim,Delta,1,Delta,1)
      gd=DDot_(nDim,Gamma,1,Delta,1)
      dBd=DDot_(nDim,Delta,1,Bd,1)
*
*     Try to avoid numerical instability
*
      Thr=1.0D-8
      If (dBd.gt.Thr .and. ABS(gd).gt.Thr) Then
         Do i = 1, nDim
            Do j = 1, nDim
               B(i,j) = B(i,j) + (Gamma(i)*Gamma(j))/gd
     &                         - (Bd(i)*Bd(j))/dBd
            End Do
         End Do
      Else If (ABS(gd).gt.Thr) Then
         Do i = 1, nDim
            Do j = 1, nDim
               B(i,j) = B(i,j) + (Gamma(i)*Gamma(j))/gd
            End Do
         End Do
      Else If (dBd.gt.Thr) Then
         Do i = 1, nDim
            Do j = 1, nDim
               B(i,j) = B(i,j) - (Bd(i)*Bd(j))/dBd
            End Do
         End Do
      End If
*
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: B',' ',B,nDim,nDim)
#endif
      Return
      End
