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
      Subroutine ISCD_MakenMat(n_max,nOsc,lNMAT,lnTabDim,
     &                         Graph2,nTabDim,nMat0)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Integer Graph2(n_max+1,n_max+1,nOsc)
      Integer nTabDim(0:lnTabDim),  nMat0(nOsc)
      Integer nvTabDim
#include "WrkSpc.fh"
#include "io_mula.fh"
C!
C!---- Initialize.
      Do i=0,maxMax_n
        nIndex(1,i)=0
        nIndex(2,i)=0
        nIndex(3,i)=0
      EndDo
      Do i=1,nOsc
        nMat0(i)=0
      EndDo
      iIndex = 0
      nTabDim(0) = iIndex
      Call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
      Call GetMem('iVec','Allo','INTE',ipiVec,nOsc)
C!
C! -- Macrocycle on iQuanta
      nD_0 = 0
      Do iQuanta=1,n_max
        do jv=1,nOsc
          iWork(ipiVec+jv-1)=0
        enddo
        iQ=-1
        iWork(ipiVec)=-1
        Call TabDim2_drv(iQuanta,nOsc,nd)
        Call TabDim2_drv(iQuanta-1,nOsc,nvTabDim)
        nd=nd-nvTabDim
C!
C! ---  Microcycle on iDet
        Do iDet=1,nD
          iWork(ipiVec)=iWork(ipiVec)+1
          iQ=iQ+1
          If (iQ.gt.iQuanta) Then
            Do i=1,nOsc-1
              if(iQ.le.iQuanta) goto 99
              iQ=iQ-iWork(ipiVec+i-1)+1
              iWork(ipiVec+i-1)=0
              iWork(ipiVec+i)=iWork(ipiVec+i)+1
           End Do
          End If
  99      Continue
          iWork(ipiVec+nOsc-1)=iQuanta-iq
          iDNR = iDetnr(iWork(ipiVec),Graph2,nosc,n_max)
          iDNR = iDNR-nD_0
          do iv=1,nOsc
            nMat0(iv)=iWork(ipiVec+iv-1)
          enddo
          nTabDim(iDNR+nD_0) = iIndex
          Call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
        End Do
        nD_0 = nD_0 + nD
      End Do
      Call GetMem('iVec','Free','INTE',ipiVec,nOsc)
      Return
      End
