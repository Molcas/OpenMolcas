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
*-- This is just an interface for the state transformation. Either
*   we use usual AO-basis route, or we take the reduced MO-basis
*   route.
*
      Subroutine StateMME(MoOrNot,nAObas,nMObas,nState,nTyp,iCi,iBigT
     &                   ,iMME,iCent,ipAvRed,Cha,Dip,Qua)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6),iCent(MxBas**2)
      Dimension Cha(MxStOT,MxQCen),Dip(MxStOT,3,MxQCen)
      Dimension Qua(MxStOT,6,MxQCen)
      Logical MoOrNot

      If(.not.MoOrNot) then
        Call StateMMEao(nAObas,nState,nTyp,iBigT,iMME,iCent,Cha,Dip,Qua)
      Else
        Call StateMMEmo(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME,iCent
     &                 ,ipAvRed,Cha,Dip,Qua)
      Endif

      Return
      End


*
*-- AO-basis route.
*
      Subroutine StateMMEao(nAObas,nState,nTyp,iBigT,iMME,iCent,Cha,Dip
     &                     ,Qua)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6),iCent(MxBas**2)
      Dimension Cha(MxStOT,MxQCen),Dip(MxStOT,3,MxQCen)
      Dimension Qua(MxStOT,6,MxQCen)

      kaunter=0
      nSize=nAObas*(nAObas+1)/2
      Call GetMem('Transition','Allo','Real',ipAOG,nSize)
      Call GetMem('OnTheWay','Allo','Real',ipO,nTyp)
*--- Loop over state pairs.
      Do 101, iS1=1,nState
        Do 102, iS2=1,iS1
          kaunter=kaunter+1
*------Collect this piece of the TDM in AO-basis.
          Call dCopy_(nSize,Work(iBigT+nSize*(kaunter-1)),iONE
     &              ,Work(ipAOG),iONE)
          kaunta=0
*------- Loop over AO-basis pairs and transform them as well as
*        distribute their multipoles. Observe that the array iCent
*        keeps track on where a certain AO-basis pair belongs.
          Do 103, iB1=1,nAObas
            Do 104, iB2=1,iB1
              PerAake=Work(ipAOG+kaunta)
              Do 105, iTyp=1,nTyp
                Work(ipO+iTyp-1)=Work(iMME(iTyp)+kaunta)*PerAake
105           Continue
              kaunta=kaunta+1
              Cha(kaunter,iCent(kaunta))=
     &           Cha(kaunter,iCent(kaunta))+Work(ipO)
              Dip(kaunter,1,iCent(kaunta))=
     &           Dip(kaunter,1,iCent(kaunta))+Work(ipO+1)
              Dip(kaunter,2,iCent(kaunta))=
     &           Dip(kaunter,2,iCent(kaunta))+Work(ipO+2)
              Dip(kaunter,3,iCent(kaunta))=
     &           Dip(kaunter,3,iCent(kaunta))+Work(ipO+3)
              Qua(kaunter,1,iCent(kaunta))=
     &           Qua(kaunter,1,iCent(kaunta))+Work(ipO+4)
              Qua(kaunter,2,iCent(kaunta))=
     &           Qua(kaunter,2,iCent(kaunta))+Work(ipO+5)
*----------- The reason why 7 and 6 are interchanged is that
*            QMSTAT uses the ordering xx,xy,yy,xz,yz,zz while
*            Seward uses the ordering xx,xy,xz,yy,yz,zz.
              Qua(kaunter,3,iCent(kaunta))=
     &           Qua(kaunter,3,iCent(kaunta))+Work(ipO+7)
              Qua(kaunter,4,iCent(kaunta))=
     &           Qua(kaunter,4,iCent(kaunta))+Work(ipO+6)
              Qua(kaunter,5,iCent(kaunta))=
     &           Qua(kaunter,5,iCent(kaunta))+Work(ipO+8)
              Qua(kaunter,6,iCent(kaunta))=
     &           Qua(kaunter,6,iCent(kaunta))+Work(ipO+9)
104         Continue
103       Continue
102     Continue
101   Continue
      Call GetMem('OnTheWay','Free','Real',ipO,nTyp)
      Call GetMem('Transition','Free','Real',ipAOG,nSize)

      Return
      End


*
*-- MO-basis route.
*
      Subroutine StateMMEmo(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME
     &                     ,iCent,ipAvRed,Cha,Dip,Qua)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6),iCent(MxBas**2)
      Dimension Cha(MxStOT,MxQCen),Dip(MxStOT,3,MxQCen)
      Dimension Qua(MxStOT,6,MxQCen)

      kaunter=0
      nSizeA=nAObas*(nAObas+1)/2
      nSizeM=nMObas*(nMObas+1)/2
      Call GetMem('Transition','Allo','Real',ipMOG,nSizeM)
      Call GetMem('SqMO','Allo','Real',ipMOG_s,nMObas**2)
      Call GetMem('TEMP','Allo','Real',iTEMP,nAObas*nMObas)
      Call GetMem('SqAO','Allo','Real',ipAOG_s,nAObas**2)
      Call GetMem('TransitionA','Allo','Real',ipAOG,nSizeA)
      Call GetMem('OnTheWay','Allo','Real',ipO,nTyp)
*
*--- Loop over state pairs.
*
      Do 101, iS1=1,nState
        Do 102, iS2=1,iS1
          kaunter=kaunter+1
*
*------- Collect the proper piece of the TDM in MO-basis.
*
          Call dCopy_(nSizeM,Work(iBigT+nSizeM*(kaunter-1)),iONE
     &              ,Work(ipMOG),iONE)
*
*------- Additional transformation step from MO to AO.
*
          Call Square(Work(ipMOG),Work(ipMOG_s),iONE,nMObas,nMObas)
          kk=0
          Do 403, i=1,nMObas
            Do 404, j=1,nMObas
              If(i.ne.j)Work(ipMOG_s+kk)=0.5d0*Work(ipMOG_s+kk)
              kk=kk+1
404         Continue
403       Continue
          Call Dgemm_('N','N',nAObas,nMObas,nMObas,ONE,Work(ipAvRed)
     &              ,nAObas,Work(ipMOG_s),nMObas,ZERO,Work(iTEMP)
     &              ,nAObas)
          Call Dgemm_('N','T',nAObas,nAObas,nMObas,ONE,Work(iTEMP)
     &              ,nAObas,Work(ipAvRed),nAObas,ZERO,Work(ipAOG_s)
     &              ,nAObas)
          kk=0
          Do 405, i=1,nAObas
            Do 406, j=1,nAObas
              If(i.ne.j)Work(ipAOG_s+kk)=2.0d0*Work(ipAOG_s+kk)
              kk=kk+1
406         Continue
405       Continue
          Call SqToTri_Q(Work(ipAOG_s),Work(ipAOG),nAObas)
*
*------- Loop over AO-basis pairs.
*
          kaunta=0
          Do 103, iB1=1,nAObas
            Do 104, iB2=1,iB1
              PerAake=Work(ipAOG+kaunta)
              Do 105, iTyp=1,nTyp
                Work(ipO+iTyp-1)=Work(iMME(iTyp)+kaunta)*PerAake
105           Continue
              kaunta=kaunta+1
              Cha(kaunter,iCent(kaunta))=
     &           Cha(kaunter,iCent(kaunta))+Work(ipO)
              Dip(kaunter,1,iCent(kaunta))=
     &           Dip(kaunter,1,iCent(kaunta))+Work(ipO+1)
              Dip(kaunter,2,iCent(kaunta))=
     &           Dip(kaunter,2,iCent(kaunta))+Work(ipO+2)
              Dip(kaunter,3,iCent(kaunta))=
     &           Dip(kaunter,3,iCent(kaunta))+Work(ipO+3)
              Qua(kaunter,1,iCent(kaunta))=
     &           Qua(kaunter,1,iCent(kaunta))+Work(ipO+4)
              Qua(kaunter,2,iCent(kaunta))=
     &           Qua(kaunter,2,iCent(kaunta))+Work(ipO+5)
*----------- The reason why 7 and 6 are interchanged is that
*            QMSTAT uses the ordering xx,xy,yy,xz,yz,zz while
*            Seward uses the ordering xx,xy,xz,yy,yz,zz.
              Qua(kaunter,3,iCent(kaunta))=
     &           Qua(kaunter,3,iCent(kaunta))+Work(ipO+7)
              Qua(kaunter,4,iCent(kaunta))=
     &           Qua(kaunter,4,iCent(kaunta))+Work(ipO+6)
              Qua(kaunter,5,iCent(kaunta))=
     &           Qua(kaunter,5,iCent(kaunta))+Work(ipO+8)
              Qua(kaunter,6,iCent(kaunta))=
     &           Qua(kaunter,6,iCent(kaunta))+Work(ipO+9)
104         Continue
103       Continue
102     Continue
101   Continue
      Call GetMem('Transition','Free','Real',ipMOG,nSizeM)
      Call GetMem('SqMO','Free','Real',ipMOG_s,nMObas**2)
      Call GetMem('TEMP','Free','Real',iTEMP,nAObas*nMObas)
      Call GetMem('SqAO','Free','Real',ipAOG_s,nAObas**2)
      Call GetMem('TransitionA','Free','Real',ipAOG,nSizeA)
      Call GetMem('OnTheWay','Free','Real',ipO,nTyp)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iCi)
      End


*
*-- The reduced MO-basis route. OBSOLOTE!! WORKS BUT IS SLOW!!!
*
      Subroutine StateMMEmo_NO(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME
     &                     ,iCent,ipAvRed,Cha,Dip,Qua)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6),iCent(MxBas**2)
      Dimension iAcc(MxMltp*(MxMltp+1)*(MxMltp+2)/6)
      Dimension Cha(MxStOT,MxQCen),Dip(MxStOT,3,MxQCen)
      Dimension Qua(MxStOT,6,MxQCen)

      kaunter=0
      nSizeA=nAObas*(nAObas+1)/2
      nSizeM=nMObas*(nMObas+1)/2
      Call GetMem('Transition','Allo','Real',ipMOG,nSizeM)
      Do 200, i=1,nTyp
        Call GetMem('Accumulate','Allo','Real',iAcc(i),nSizeA)
200   Continue
*--- Loop over state pairs.
      Do 201, iS1=1,nState
        Do 202, iS2=1,iS1
          kaunter=kaunter+1
*------- Collect the proper piece of the TDM in MO-basis.
          Call dCopy_(nSizeM,Work(iBigT+nSizeM*(kaunter-1)),iONE
     &              ,Work(ipMOG),iONE)
*------- Loop over centres in molecule. This is now necessary since
*        MOs are contrary to AOs not localized, hence the simple
*        construction with iCent used above, can not be used here.
          Do 203, iCentre=1,iCi
            Do 2031, i=1,nTyp
              Call dCopy_(nSizeA,[ZERO],iZERO,Work(iAcc(i)),iONE)
2031        Continue
            kaunta=0
*--------- Loop over AO-basis pairs.
            Do 204, iB1=1,nAObas
              Do 205, iB2=1,iB1
                kaunta=kaunta+1
*------------- If this basis pair belongs to the given centre,
*              accumulate all multipoles.
                If(iCent(kaunta).eq.iCentre) then
                  Do 2051, i=1,nTyp
                    Work(iAcc(i)+kaunta-1)=Work(iAcc(i)+kaunta-1)
     &                                    +Work(iMME(i)+kaunta-1)
2051              Continue
                Endif
205           Continue
204         Continue
*--------- Transform the MME on this centre to MO-basis. We hence
*          get the contribution to the density distributed on a
*          specific centre.
            Call MMEtoRMO(nAObas,nMObas,ipAvRed,iAcc)
*--------- Ordinary evaluations of expectation values.
            Cha(kaunter,iCentre)=Ddot_(nSizeM,Work(iAcc(1)),iONE
     &                                  ,Work(ipMOG),iONE)
            Dip(kaunter,1,iCentre)=Ddot_(nSizeM,Work(iAcc(2)),iONE
     &                                  ,Work(ipMOG),iONE)
            Dip(kaunter,2,iCentre)=Ddot_(nSizeM,Work(iAcc(3)),iONE
     &                                  ,Work(ipMOG),iONE)
            Dip(kaunter,3,iCentre)=Ddot_(nSizeM,Work(iAcc(4)),iONE
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,1,iCentre)=Ddot_(nSizeM,Work(iAcc(5)),iONE
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,2,iCentre)=Ddot_(nSizeM,Work(iAcc(6)),iONE
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,3,iCentre)=Ddot_(nSizeM,Work(iAcc(8)),iONE
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,4,iCentre)=Ddot_(nSizeM,Work(iAcc(7)),iONE
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,5,iCentre)=Ddot_(nSizeM,Work(iAcc(9)),iONE
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,6,iCentre)=Ddot_(nSizeM,Work(iAcc(10)),iONE
     &                                  ,Work(ipMOG),iONE)
203       Continue
202     Continue
201   Continue
      Call GetMem('Transition','Free','Real',ipMOG,nSizeM)
      Do 206, i=1,nTyp
        Call GetMem('Accumulate','Free','Real',iAcc(i),nSizeA)
206   Continue

      Return
      End
