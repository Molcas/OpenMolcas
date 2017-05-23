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
*      Subroutine ContRASBas(nBas,nStatePrim,iNonH,iNonS,Eigis)
      Subroutine ContRASBas(nBas,nStatePrim,iNonH,iNonS,iEig2)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension nBas(MxSym)!,Eigis(MxState,MxState)

*
*--- Hi y'all
*
      Write(6,*)'     ----- Constructing CASSI eigenstates.'

*
*--- Diagonalize overlap matrix.
*
      Call GetMem('EigV1','Allo','Real',iEig1,nStatePrim**2)
      kaunter=0
      Do 1, i=1,nStatePrim
        Do 2, j=1,nStatePrim
          If(i.eq.j) then
            Work(iEig1+kaunter)=ONE
          Else
            Work(iEig1+kaunter)=ZERO
          Endif
          kaunter=kaunter+1
2       Continue
1     Continue
      Call Jacob(Work(iNonS),Work(iEig1),nStatePrim,nStatePrim)
      If(iPrint.ge.15) then
        Call TriPrt('Diagonal RASSCF overlap matrix',' '
     &             ,Work(iNonS),nStatePrim)
      Endif

*
*--- Construct TS^(-1/2) for canonical orthogonalization.
*
      ii=0
      Do 4, i=1,nStatePrim
        ii=ii+i
        x=1.0D00/Sqrt(Max(1.0D-14,Work(iNonS-1+ii)))
        Do 5, k=1,nStatePrim
          ind=k+nStatePrim*(i-1)-1
          Work(iEig1+ind)=x*Work(iEig1+ind)
5       Continue
4     Continue

*
*--- Make reductions if requested.
*
      Call GetMem('RedEigV1','Allo','Real',iEig2,nStatePrim**2)
      iT=0
      If(ContrStateB) then
        Do 61, iS=1,nStatePrim
          kaunt=iS*(iS+1)/2-1
          sss=Work(iNonS+kaunt)
          If(sss.gt.ThrsCont) then
            iT=iT+1
            call dcopy_(nStatePrim,Work(iEig1+nStatePrim*(iS-1)),iONE
     &                ,Work(iEig2+nStatePrim*(iT-1)),iONE)
          Endif
61      Continue
        nStateRed=iT
        Write(6,6199)'  ----- Contraction:',nStatePrim,' ---> '
     &                                     ,nStateRed
6199    Format(A,I3,A,I3)
      Else
        call dcopy_(nStatePrim**2,Work(iEig1),iONE,Work(iEig2),iONE)
        nStateRed=nStatePrim
      Endif

*
*--- Transform H and diagonalize in the original basis.
*
      nTri=nStateRed*(nStateRed+1)/2
      Call GetMem('TEMP','Allo','Real',iTEMP,nStatePrim**2)
      Call GetMem('SqH','Allo','Real',iSqH,nStatePrim**2)
      Call GetMem('RedHSq','Allo','Real',iRedHSq,nStateRed**2)
      Call GetMem('RedHTr','Allo','Real',iRedHTr,nTri)
      Call Square(Work(iNonH),Work(iSqH),iONE,nStatePrim,nStatePrim)
      Call Dgemm_('N','N',nStatePrim,nStateRed,nStatePrim,ONE,Work(iSqH)
     &          ,nStatePrim,Work(iEig2),nStatePrim,ZERO,Work(iTEMP)
     &          ,nStatePrim)
      Call Dgemm_('T','N',nStateRed,nStateRed,nStatePrim,ONE,Work(iEig2)
     &          ,nStatePrim,Work(iTEMP),nStatePrim,ZERO,Work(iRedHSq)
     &          ,nStateRed)
      Call SqToTri_Q(Work(iRedHSq),Work(iRedHTr),nStateRed)
      Call Jacob(Work(iRedHTr),Work(iEig2),nStateRed,nStatePrim)
      Call JacOrd(Work(iRedHTr),Work(iEig2),nStateRed,nStatePrim)

*
*--- At this stage we have eigenvectors to the CASSI states and their
*    eigenenergies, hence time to construct the first H_0 and store the
*    eigenvectors for subsequent transformations.
*
      kaunter=0
      nLvlInd=1
      Do 55, iState=1,nStateRed
        Do 56, jState=1,iState
          kaunter=kaunter+1
          HMatState(kaunter)=ZERO
56      Continue
        HMatState(kaunter)=Work(iRedHTr+kaunter-1)
*---- If requested, introduce level-shift of states.
        If(nLvlShift.ne.0) then
          If(iState.eq.iLvlShift(nLvlInd)) then
            HMatState(kaunter)=HMatState(kaunter)+dLvlShift(nLvlInd)
            nLvlInd=nLvlInd+1
          Endif
        Endif
55    Continue

*
*--- Print.
*
      If(iPrint.ge.10) then
        Call TriPrt('RASSI Hamiltonian',' ',HMatState,nStateRed)
        Write(6,*)
        Call RecPrt('RASSI eigenvectors',' ',Work(iEig2),nStatePrim
     &             ,nStateRed)
      Endif

*
*--- Deallocate.
*
      Call GetMem('EigV1','Free','Real',iEig1,nStatePrim**2)
      Call GetMem('TEMP','Free','Real',iTEMP,nStatePrim**2)
      Call GetMem('SqH','Free','Real',iSqH,nStatePrim**2)
      Call GetMem('RedHSq','Free','Real',iRedHSq,nStateRed**2)
      Call GetMem('RedHTr','Free','Real',iRedHTr,nTri)

*
*--- OBSERVE! CAUTION! ATTENTION! The variable nState is defined.
*
      nState=nStateRed

*
*--- No parasan!
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nBas)
      End
