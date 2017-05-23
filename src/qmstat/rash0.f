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
*
*-- In this routine H_0 in RASSI basis is constructed, possibly with
*   external perturbation added on.
*
      Subroutine RasH0(nB)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension DiagH0(MxState)

      nBTri=nB*(nB+1)/2
      nTri=nState*(nState+1)/2
      If(.not.AddExt) then
        kaunter=0
        Do 1, i=1,nState
          Do 2, j=1,i
            kaunter=kaunter+1
2         Continue
          DiagH0(i)=HmatState(kaunter)
1       Continue
        Write(6,*)'     -----RASSI H_0 eigenvalues:'
        Write(6,99)(DiagH0(k),k=1,nState)
      Else
*
*---- Collect one-electron perturbations.
*
        Lu_One=49
        Lu_One=IsFreeUnit(Lu_One)
        Call OpnOne(irc,0,'ONEINT',Lu_One)
        Call GetMem('AOExt','Allo','Real',ipAOx,nBTri+4)
        Do 101, iExt=1,nExtAddOns
          irc=-1
          iopt=0
          iSmLbl=0
          Call RdOne(irc,iopt,ExtLabel(iExt),iCompExt(iExt),Work(ipAOx)
     &              ,iSmLbl)
          call dscal_(nBTri,ScalExt(iExt),Work(ipAOx),iONE)
          If(irc.ne.0) then
            Write(6,*)
            Write(6,*)'ERROR when reading ',ExtLabel(iExt),'.'
            Write(6,*)'Have Seward computed this integral?'
            Call Quit(_RC_IO_ERROR_READ_)
          Endif
*
*---- We need to know in which basis the TDM is and then transform
*     the one-electron integrals to RASSI-basis.
*
          If(.not.MoAveRed) then
            Call GetMem('Transition','Allo','Real',ipAOG,nBTri)
            kaunter=0
            Do 102, iS1=1,nState
              Do 103, iS2=1,iS1
                call dcopy_(nBTri,Work(iBigT+nBTri*kaunter),iONE
     &                    ,Work(ipAOG),iONE)
                Element=Ddot_(nBTri,Work(ipAOG),iONE,Work(ipAOx),iONE)
                kaunter=kaunter+1
                HmatState(kaunter)=HmatState(kaunter)+Element
103           Continue
102         Continue
            Call GetMem('Transition','Free','Real',ipAOG,nBTri)
          Else
            nSize=nRedMO*(nRedMO+1)/2
            Call GetMem('Transition','Allo','Real',ipMOG,nSize)
            Call GetMem('AUX','Allo','Real',iAUX,nRedMO*nB)
            Call GetMem('SquareAO','Allo','Real',iSqAO,nB**2)
            Call GetMem('SquareMO','Allo','Real',iSqMO,nRedMO**2)
            Call GetMem('MOExt','Allo','Real',ipMOx,nSize)
            Call Square(Work(ipAOx),Work(iSqAO),iONE,nB,nB)
            Call Dgemm_('T','N',nRedMO,nB,nB,ONE,Work(ipAvRed)
     &                ,nB,Work(iSqAO),nB,ZERO,Work(iAUX),nRedMO)
            Call Dgemm_('N','N',nRedMO,nRedMO,nB,ONE,Work(iAUX)
     &                ,nRedMO,Work(ipAvRed),nB,ZERO,Work(iSqMO),nRedMO)
            Call SqToTri_Q(Work(iSqMO),Work(ipMOx),nRedMO)
            kaunter=0
            Do 104, iS1=1,nState
              Do 105, iS2=1,nState
                call dcopy_(nSize,Work(iBigT+nSize*kaunter),iONE
     &                    ,Work(ipMOG),iONE)
                Element=Ddot_(nSize,Work(ipMOG),iONE,Work(ipMOx),iONE)
                kaunter=kaunter+1
                HmatState(kaunter)=HmatState(kaunter)+Element
105           Continue
104         Continue
            Call GetMem('Transition','Free','Real',ipMOG,nSize)
            Call GetMem('AUX','Free','Real',iAUX,nRedMO*nB)
            Call GetMem('SquareAO','Free','Real',iSqAO,nB**2)
            Call GetMem('SquareMO','Free','Real',iSqMO,nRedMO**2)
            Call GetMem('MOExt','Free','Real',ipMOx,nSize)
          Endif
101     Continue
        Call GetMem('AOExt','Free','Real',ipAOx,nBTri+4)
        Call ClsOne(irc,Lu_One)
*
*-- If sufficient print level, print HmatState with perturbation added.
*
        If(iPrint.ge.5) then
          Write(6,*)
          Call TriPrt('H_0+External perturbation',' ',HmatState,nState)
        Endif
      Endif

99    Format('            ',9(F12.7,'  '))

      Return
      End
