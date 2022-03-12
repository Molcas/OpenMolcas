!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
!-- The reduced MO-basis route. OBSOLOTE!! WORKS BUT IS SLOW!!!
!
      Subroutine StateMMEmo_NO(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME &
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
!--- Loop over state pairs.
      Do 201, iS1=1,nState
        Do 202, iS2=1,iS1
          kaunter=kaunter+1
!------- Collect the proper piece of the TDM in MO-basis.
          Call dCopy_(nSizeM,Work(iBigT+nSizeM*(kaunter-1)),iONE        &
     &              ,Work(ipMOG),iONE)
!------- Loop over centres in molecule. This is now necessary since
!        MOs are contrary to AOs not localized, hence the simple
!        construction with iCent used above, can not be used here.
          Do 203, iCentre=1,iCi
            Do 2031, i=1,nTyp
              Call dCopy_(nSizeA,[ZERO],iZERO,Work(iAcc(i)),iONE)
2031        Continue
            kaunta=0
!--------- Loop over AO-basis pairs.
            Do 204, iB1=1,nAObas
              Do 205, iB2=1,iB1
                kaunta=kaunta+1
!------------- If this basis pair belongs to the given centre,
!              accumulate all multipoles.
                If(iCent(kaunta).eq.iCentre) then
                  Do 2051, i=1,nTyp
                    Work(iAcc(i)+kaunta-1)=Work(iAcc(i)+kaunta-1)       &
     &                                    +Work(iMME(i)+kaunta-1)
2051              Continue
                Endif
205           Continue
204         Continue
!--------- Transform the MME on this centre to MO-basis. We hence
!          get the contribution to the density distributed on a
!          specific centre.
            Call MMEtoRMO(nAObas,nMObas,ipAvRed,iAcc)
!--------- Ordinary evaluations of expectation values.
            Cha(kaunter,iCentre)=Ddot_(nSizeM,Work(iAcc(1)),iONE        &
     &                                  ,Work(ipMOG),iONE)
            Dip(kaunter,1,iCentre)=Ddot_(nSizeM,Work(iAcc(2)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Dip(kaunter,2,iCentre)=Ddot_(nSizeM,Work(iAcc(3)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Dip(kaunter,3,iCentre)=Ddot_(nSizeM,Work(iAcc(4)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,1,iCentre)=Ddot_(nSizeM,Work(iAcc(5)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,2,iCentre)=Ddot_(nSizeM,Work(iAcc(6)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,3,iCentre)=Ddot_(nSizeM,Work(iAcc(8)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,4,iCentre)=Ddot_(nSizeM,Work(iAcc(7)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,5,iCentre)=Ddot_(nSizeM,Work(iAcc(9)),iONE      &
     &                                  ,Work(ipMOG),iONE)
            Qua(kaunter,6,iCentre)=Ddot_(nSizeM,Work(iAcc(10)),iONE     &
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
