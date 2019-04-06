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
      Subroutine RasRasTrans(nB,nStatePrim,iEig2,iPrint)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension iTocBig(MxStOT)

      Character*30 OutLine

*
*--- Guten Tag.
*
      Write(6,*)'     ----- Transform from non-orthogonal RASSCF states'
     &//' to orthogonal RASSI states.'

*
*--- Set zeros and decide if transformation is at all possible.
*
      LuIn=66
      kaunt=0
      iSkutt=0
      iDisk=0
      Call DaName(LuIn,RassiM)
      nTriSP=nStatePrim*(nStatePrim+1)/2
      Call iDaFile(LuIn,2,iTocBig,nTriSP,iDisk)
      nSize=nB*(nB+1)/2
      nTriS=nState*(nState+1)/2
      nSizeBig=nSize*nTriS
      nSizeBigPrim=nSize*nTriSP
      Call GetMem('HOWMUCH','Max','Real',ipMAX,MEMMAX)
*
*-- This means that we do not have memory enough for TDM in contracted
*   form. Then there is no use to proceed at all.
*
      If(MEMMAX.le.nSizeBig) then
        Write(6,*)
        Write(6,*)'The transition density matrix is too big to put in m'
     &//'emory!'
        Write(6,*)'Either,'
        Write(6,*)'       (1) increase MOLCAS_MEM,'
        Write(6,*)'       (2) contract number of states further.'
        Call Quit(_RC_GENERAL_ERROR_)
*
*-- Here we go if there is enough memory for an in core transformation.
*
      Elseif(MEMMAX.ge.(nSizeBig+nSizeBigPrim+nTriSP+nTriS+nStatePrim**2
     &                 +nState*nStatePrim+nState**2)) then
        Call GetMem('ALLES','Allo','Real',iBigT,nSizeBig)
        Call GetMem('ALLESin','Allo','Real',iBigV,nSizeBigPrim)
        Call GetMem('Int1','Allo','Real',iInt1,nTriSP)
        Call GetMem('Int2','Allo','Real',iInt2,nTriS)
        Call GetMem('Square1','Allo','Real',iSnt1,nStatePrim**2)
        Call GetMem('Square2','Allo','Real',iSnt2,nState*nStatePrim)
        Call GetMem('Square3','Allo','Real',iSnt3,nState**2)
        call dcopy_(nSizeBig,[ZERO],iZERO,Work(iBigT),iONE)
        kaunt=0
        Do 78, i=1,nStatePrim
          Do 79, j=1,i
            kaunt=kaunt+1
            iDisk=iTocBig(kaunt)
            Call dDaFile(LuIn,2,Work(iBigV+(kaunt-1)*nSize),nSize,iDisk)
79        Continue
78      Continue
*
*---- A lot of printing of TDM if requested.
*
        If(iPrint.ge.25) then
          kaunt=0
          Do 17, i=1,nStatePrim
            Do 18, j=1,i
              Write(OutLine,'(A,I3,I3)')'TDM, Piece ',i,j
              Call TriPrt(OutLine,' ',Work(iBigV+kaunt),nB)
              kaunt=kaunt+nSize
18          Continue
17        Continue
        Endif
*
*---- Proceed with transformation.
*
        kaunt=0
        Do 10001, iBas=1,nB
          Do 10002, jBas=1,iBas
            call dcopy_(nTriSP,Work(iBigV+kaunt),nSize,Work(iInt1),iONE)
            Call Square(Work(iInt1),Work(iSnt1),iONE,nStatePrim
     &                                              ,nStatePrim)
            Call Dgemm_('T','N',nState,nStatePrim,nStatePrim,ONE
     &                ,Work(iEig2),nStatePrim,Work(iSnt1),nStatePrim
     &                ,ZERO,Work(iSnt2),nState)
            Call Dgemm_('N','N',nState,nState,nStatePrim,ONE,Work(iSnt2)
     &                ,nState,Work(iEig2),nStatePrim,ZERO,Work(iSnt3)
     &                ,nState)
            Call SqToTri_Q(Work(iSnt3),Work(iInt2),nState)
            call dcopy_(nTriS,Work(iInt2),iONE,Work(iBigT+kaunt),nSize)
            kaunt=kaunt+1
10002     Continue
10001   Continue
        Call GetMem('ALLESin','Free','Real',iBigV,nSizeBigPrim)
        Call GetMem('Int1','Free','Real',iInt1,nTriSP)
        Call GetMem('Int2','Free','Real',iInt2,nTriS)
        Call GetMem('Square1','Free','Real',iSnt1,nStatePrim**2)
        Call GetMem('Square2','Free','Real',iSnt2,nState*nStatePrim)
        Call GetMem('Square3','Free','Real',iSnt3,nState**2)
*
*-- Here we go if both TDM's can not be put in memory. Might be a bit
*   slow due to its nested nature with repeated IO.
*
      Else
        Call GetMem('ALLES','Allo','Real',iBigT,nSizeBig)
        Call GetMem('AOGamma','Allo','Real',ipAOG,nSize)
        call dcopy_(nSizeBig,[ZERO],iZERO,Work(iBigT),iONE)
        Do 11001, iiS=1,nStatePrim
          Do 11002, jjS=1,nStatePrim
            If(iiS.le.jjS) then
              indypop=jjS*(jjS+1)/2-jjS+iiS
            Else
              indypop=iiS*(iiS+1)/2-iiS+jjS
            Endif
            iDisk=iTocBig(indypop)
            Call dDaFile(LuIn,2,Work(ipAOG),nSize,iDisk)
            kaunter=0
            Do 11003, iB=1,nB
              Do 11004, jB=1,iB
                Do 11005, iS=1,nState
                  Do 11006, jS=1,iS
                    index=(iS*(iS-1)/2+jS-1)*nSize
                    index=index+kaunter
                    Work(iBigT+index)=Work(iBigT+index)
     &                +Work(iEig2+iiS-1+(iS-1)*nStatePrim)
     &                *Work(iEig2+jjS-1+(jS-1)*nStatePrim)
     &                *Work(ipAOG+kaunter)
11006             Continue
11005           Continue
                kaunter=kaunter+1
11004         Continue
11003       Continue
11002     Continue
11001   Continue
        Call GetMem('AOGamma','Free','Real',ipAOG,nSize)
      Endif

*
*--- Deallocations and finish up.
*
      Call GetMem('RedEigV1','Free','Real',iEig2,nStatePrim**2)
      Call DaClos(LuIn)

      Return
      End
