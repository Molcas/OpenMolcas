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
* Copyright (C) Anders Ohrn                                            *
************************************************************************
      Subroutine ParaRoot(Ract,BetaBol,Etot,CalledBefore,SampleThis)
************************************************************
*
*   <DOC>
*     <Name>ParaRoot</Name>
*     <Syntax>Call ParaRoot(Ract,BetaBol,Etot,CalledBefore,SampleThis)</Syntax>
*     <Arguments>
*       \Argument{Ract}{}{}{in}
*       \Argument{BetaBol}{}{}{in}
*       \Argument{Etot}{}{}{in}
*       \Argument{CalledBefore}{}{}{inout}
*       \Argument{SampleThis}{}{}{out}
*     </Arguments>
*     <Purpose>Manage the parallel tempering routine</Purpose>
*     <Dependencies></Dependencies>
*     <Author>A.Ohrn</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      If our system is difficult and has small transition elements
*      in the Markov chain, we can use the parallel tempering to
*      boost sampling. This routine is the root for this; it
*      mainly handles the various configurations for the different
*      temperature ensembles; also, manages the ensemble switch.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      External Ranf
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"
#include<constants.fh>

*     Parameter (BoltzK=1.0d-3*CONST_BOLTZMANN_/CONV_AU_TO_KJ_)
      Dimension iPermutation(2,MxParT)
      Dimension CordstTEMP(MxCen*MxPut,3)
      Logical CalledBefore,WeiterBitte,Accept,SampleThis
      Save iTemp

      BoltzK=1.0d-3*CONST_BOLTZMANN_/CONV_AU_TO_KJ_
      Dum1=0.0d0
*
*-- If this is first time to call on this routine.
*
      If(.not.CalledBefore) then
        iTemp=0
        CalledBefore=.true.
      Endif

*
*-- See what to do.
*
999   Continue
      If(iTemp.lt.nTemp) then
        WeiterBitte=.true.
        iTemp=iTemp+1
        Write(6,*)
        Write(6,*)'    Run a new temperature ensemble...',iTemp

*----A logical variable to make parallel tempering sampling correct.
        If(iTemp.eq.1) then
          SampleThis=.true.
        Else
          SampleThis=.false.
        Endif

      Else
        WeiterBitte=.false.
        iTemp=0
        Write(6,*)
        Write(6,*)'    Evaluate temperature ensemble interchanges.'
      Endif

*
*-- If we are to run a new ensemble.
*
      If(WeiterBitte) then
        iLuStIn=8+nStFilT(iTemp)
        iLuStUt=16+nStFilT(iTemp)
        Write(StFilIn(6:6),'(i1.1)')nStFilT(iTemp)
        Write(StFilUt(6:6),'(i1.1)')nStFilT(iTemp)

*---- Collect coordinates from proper startfile.
        Call Get8(Ract,Etot)

*---- Set temperature.
        BetaBol=1.0d0/(ParaTemps(iTemp)*BoltzK)

*
*-- If we are to attempt interchanges.
*
      Else

        Do 10, iPa=1,MxParT
          iPermutation(1,iPa)=iPa
          iPermutation(2,iPa)=iPa
10      Continue

*-- Construct permutations, treat nTemp.eq.2 as special case, the others
*   are obtained with general algorithm.
        If(nTemp.eq.2) then
          iPermutation(2,1)=2
          iPermutation(2,2)=1
          Go To 101
        Endif

        PerType=Ranf(iseed)
        If(PerType.lt.0.5D+0) then

          If(Mod(nTemp,2).eq.1) then
            mTemp=nTemp-1
          Else
            mTemp=nTemp
          Endif

*------ Construct permutation for odd iMac
          Do 12, iPa=1,mTemp,2
            iPermutation(2,iPa)=iPermutation(1,iPa+1)
            iPermutation(2,iPa+1)=iPermutation(1,iPa)
12        Continue

        Else

          mTemp=2*((nTemp-1)/2)
*------ Contruct permutation for even iMac
          Do 13, iPa=2,mTemp,2
            iPermutation(2,iPa)=iPermutation(1,iPa+1)
            iPermutation(2,iPa+1)=iPermutation(1,iPa)
13        Continue
        Endif

101     Continue

*
*-- Now attempt interchange.
*
        iEnsemb=1
2001    Continue
          If(iPermutation(1,iEnsemb).eq.iPermutation(2,iEnsemb)) then
            iEnsemb=iEnsemb+1
            Go To 2001
          Endif

*------ Collect energies for the permutations.
          iLuStIn=8+nStFilT(iPermutation(1,iEnsemb))
          iLuStUt=16+nStFilT(iPermutation(1,iEnsemb))
          Write(StFilIn(6:6),'(i1.1)')nStFilT(iPermutation(1,iEnsemb))
          Write(StFilUt(6:6),'(i1.1)')nStFilT(iPermutation(1,iEnsemb))
          Call Get8(Dum,E1)
          iLuStIn=8+nStFilT(iPermutation(2,iEnsemb))
          iLuStUt=16+nStFilT(iPermutation(2,iEnsemb))
          Write(StFilIn(6:6),'(i1.1)')nStFilT(iPermutation(2,iEnsemb))
          Write(StFilUt(6:6),'(i1.1)')nStFilT(iPermutation(2,iEnsemb))
          Call Get8(Dum,E2)
          T1=ParaTemps(iPermutation(1,iEnsemb))
          T2=ParaTemps(iPermutation(2,iEnsemb))
          B1=1.0d0/(BoltzK*T1)
          B2=1.0d0/(BoltzK*T2)

*------ Make the Metropolis thing.
          BigDelta=(B2-B1)*(E2-E1)
          Expe=exp(BigDelta)
          Accept=.true.
          If(Expe.lt.1.0D+0) then
            Expran=ranf(iseed)
            If(Expe.lt.Expran) Accept=.false.
          Endif

          If(Accept) then
*-----Ct=C2
            iLuStIn=8+nStFilT(iPermutation(2,iEnsemb))
            iLuStUt=16+nStFilT(iPermutation(2,iEnsemb))
            Write(StFilIn(6:6),'(i1.1)')nStFilT(iPermutation(2,iEnsemb))
            Write(StFilUt(6:6),'(i1.1)')nStFilT(iPermutation(2,iEnsemb))
            Call Get8(R2,E2)
            Do 221, i=1,3
              Do 222, j=1,nCent*nPart
                CordstTEMP(j,i)=Cordst(j,i)
222           Continue
221         Continue
*-----C2=C1
            iLuStIn=8+nStFilT(iPermutation(1,iEnsemb))
            iLuStUt=16+nStFilT(iPermutation(1,iEnsemb))
            Write(StFilIn(6:6),'(i1.1)')nStFilT(iPermutation(1,iEnsemb))
            Write(StFilUt(6:6),'(i1.1)')nStFilT(iPermutation(1,iEnsemb))
            Call Get8(R1,E1)
            iLuStIn=8+nStFilT(iPermutation(2,iEnsemb))
            iLuStUt=16+nStFilT(iPermutation(2,iEnsemb))
            Write(StFilIn(6:6),'(i1.1)')nStFilT(iPermutation(2,iEnsemb))
            Write(StFilUt(6:6),'(i1.1)')nStFilT(iPermutation(2,iEnsemb))
            Call Put8(R1,E1,Dum1,Dum1,Dum1)
*-----C1=Ct
            Do 223, i=1,3
              Do 224, j=1,nCent*nPart
                Cordst(j,i)=CordstTEMP(j,i)
224           Continue
223         Continue
            iLuStIn=8+nStFilT(iPermutation(1,iEnsemb))
            iLuStUt=16+nStFilT(iPermutation(1,iEnsemb))
            Write(StFilIn(6:6),'(i1.1)')nStFilT(iPermutation(1,iEnsemb))
            Write(StFilUt(6:6),'(i1.1)')nStFilT(iPermutation(1,iEnsemb))
            Call Put8(R2,E2,Dum1,Dum1,Dum1)
          Endif

          iEnsemb=iEnsemb+2
          If(Accept)Write(6,*)'            accepted!'
          If(.not.Accept)Write(6,*)'            not accepted!'
        If(iEnsemb.lt.nTemp) Go To 2001
        Write(6,*)

*
*-- Do some stuff before exit. The reason we go back up is that this
*   way we will collect the right coordinates from first startfile.
*   Observe that iTemp has been reset hence we are back at the
*   square one again.
*
        Go To 999

      Endif

      Return
      End
