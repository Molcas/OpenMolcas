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
      Real*8, Parameter :: Thr=1.0D-8
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: B',' ',B,nDim,nDim)
*     Call RecPrt('DFP: Bd',' ',Bd,1,nDim)
      Call RecPrt('DFP: Gamma',' ',Gamma,1,nDim)
      Call RecPrt('DFP: Delta',' ',Delta,1,nDim)
#endif
      Call DGEMM_('N','N',
     &            nDim,1,nDim,
     &            1.0d0,B,nDim,
     &                  Delta,nDim,
     &            0.0d0,Bd,nDim)
      gd=DDot_(nDim,Gamma,1,Delta,1)
      dBd=DDot_(nDim,Delta,1,Bd,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: Bd',' ',Bd,1,nDim)
      Write (6,*) 'gd=',gd
      Write (6,*) 'dBd=',dBd
      Write (6,*) 'Thr=',Thr
#endif
      If (gd<0.0D0) Then
         Call MSP(B,Gamma,Delta,nDim)
      Else
*
         Do i = 1, nDim
            Do j = 1, nDim
               B(i,j) = B(i,j) + (Gamma(i)*Gamma(j))/gd
     &                         - (Bd(i)*Bd(j))/dBd
            End Do
         End Do
      End If
*
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: B',' ',B,nDim,nDim)
#endif
      Return
      End
