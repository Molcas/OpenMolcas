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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine AOAdd_NQ(AOInt,iBas,iBas_Eff,jBas,jBas_Eff,PrpInt,nPrp,
     &                 iCmp,jCmp,iShell,jShell,iAO,jAO)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
************************************************************************
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
      Real*8 AOInt(0:iBas_Eff-1,0:jBas_Eff-1,iCmp,jCmp), PrpInt(nPrp)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*
#ifdef _DEBUG_
      Call qEnter('AOAdd')
      Call RecPrt(' In AOAdd:AOInt',' ',AOInt,iBas_Eff*jBas_Eff,
     &            iCmp*jCmp)
      Write (6,*) 'iBas_Eff,jBas_Eff,iCmp,jCmp=',iBas_Eff,jBas_Eff,
     &            iCmp,jCmp
      Write (6,*) 'nPrp=',nPrp
#endif
*
      iAdd = iBas-iBas_Eff
      jAdd = jBas-jBas_Eff
      Do i1 = 1, iCmp
         iSO1=iAOtSO(iAO+i1,0)
*
         jjMx = jCmp
         If (iShell.eq.jShell) jjMx = i1
         Do i2 = 1, jjMx
            iSO2=iAOtSO(jAO+i2,0)
*
            Do indAO1_Eff = 0, iBas_Eff-1
               indAO1=indAO1_Eff+iAdd
               Indi = iSO1+indAO1
*
               jBsMax = jBas_Eff - 1
               If (iSO1.eq.iSO2) jBsMax=indAO1_Eff
               Do indAO2_Eff = 0, jBsMax
                  indAO2=indAO2_Eff+jAdd
                  Indj = iSO2+indAO2
#ifdef _DEBUG_
                  Write (6,*) 'iC,jC,iB,jB=',i1,i2,indAO1+1,
     &                                             indAO2+1
                  Write (6,*) 'Indi,Indj=',Indi,Indj
#endif
*
*                 Add one matrix element
*
                  PrpInt(iTri(Indi,Indj))=PrpInt(iTri(Indi,Indj))
     &                   + AOInt(indAO1_Eff,indAO2_Eff,i1,i2)
*
               End Do ! indAO2
            End Do   ! indAO1
         End Do      ! i2
      End Do         ! i1
*
#ifdef _DEBUG_
      Call GetMem(' Exit AOAdd','CHECK','REAL',iDum,iDum)
      Call qExit('AOAdd')
#endif
      Return
      End
