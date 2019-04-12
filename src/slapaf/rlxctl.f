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
      Subroutine RlxCtl(iStop)
      Use Chkpnt
      Implicit Real*8 (a-h,o-z)
************************************************************************
*     Program for determination of the new molecular geometry          *
************************************************************************
#include "info_slapaf.fh"
      Parameter(nLabels=10*MxAtom,nLbl=10*MxAtom)
#include "real.fh"
#include "WrkSpc.fh"
#include "nadc.fh"
#include "weighting.fh"
#include "db.fh"
#include "print.fh"
      Logical Numerical, GoOn, PrQ, TSReg,
     &        Do_ESPF, Just_Frequencies, Found
      Character*8 GrdLbl, StpLbl, Labels(nLabels), Lbl(nLbl)
      Character*1 Step_trunc
      Integer AixRm, iNeg(2)
*
      Lu=6
      Call QEnter('RlxCtl')
      iRout = 32
      iPrint=nPrint(iRout)
      StpLbl=' '
      GrdLbl=' '
      Just_Frequencies=.False.
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*
      Call RdCtl_Slapaf(iRow,iInt,nFix,LuSpool,.False.)
*
      Call Close_LuSpool(LuSpool)
*
      Call Chkpnt_open()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Request_Alaska.or.Request_RASSI) Then
*
*        Alaska/RASSI only
         iStop=3

         Go To 999
      End If

      If (isFalcon) Then
         iStop=1
         Goto 999
      End if
*                                                                      *
************************************************************************
*                                                                      *
      jPrint=nPrint(iRout)
*
      If (nLbl.lt.nBVec) Then
         Call WarningMessage(2,'Error in RlxCtl')
         Write (Lu,*)
         Write (Lu,*) '**********************'
         Write (Lu,*) ' ERROR: nLbl.lt.nBVec '
         Write (Lu,*) ' nLbl=',nLbl
         Write (Lu,*) ' nBVec=',nBVec
         Write (Lu,*) '**********************'
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      PrQ= .Not.Request_Alaska
*     PrQ=.True.
*                                                                      *
************************************************************************
*                                                                      *
      If (lCtoF .AND. PrQ) Call Def_CtoF(.False.,Work(ipCM),
     &                                   nsAtom,AtomLbl,
     &                                   Work(ipCoor),nSym,iOper,
     &                                   jStab,nStab)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the Wilson B-matrices, this describe the transformations
*     between internal and Cartesian coordinates. Values of the
*     Internal coordinates are computed too.
*
      HSet=.True.
      BSet=.True.
      kIter=iter
*
*---- Compute number of steps for numerical differentiation
*
      NmIter=1
      If (lRowH)  NmIter=nRowH+1     ! Numerical only for some rows
      If (lNmHss) NmIter=2*mInt+1    ! Full numerical
      If (Cubic)  NmIter=2*mInt**2+1 ! Full cubic
*
      If (lTherm .and. iter.EQ.1) then
         Call Put_dArray('Initial Coordinates',Work(ipCoor),3*nsAtom)
      EndIf
*
*---- Fix the definition of internal during numerical differentiation
      If (lNmHss.and.iter.lt.NmIter.and.iter.ne.1) nPrint(122)=5
*
*---- Do not overwrite numerical Hessian
      If ((lNmHss.or.lRowH)
     &    .and.(iter.gt.NmIter.or.iter.lt.NmIter)) HSet = .False.
*
*---- Set logical to indicate status during numerical differentiation
      Numerical = lNmHss .and.iter.le.NmIter .and.iter.ne.1
*
      If (Numerical) nWndw=NmIter
      Call BMtrx(iRow,nBVec,ipB,nsAtom,mInt,ipqInt,Lbl,
     &           Work(ipCoor),nDimBC,Work(ipCM),AtomLbl,nSym,
     &           iOper,Smmtrc,Degen,BSet,HSet,iter,ipdqInt,ipShf,
     &           Work(ipGx),Work(ipCx),mTtAtm,iWork(ipANr),iOptH,
     &           User_Def,nStab,jStab,Curvilinear,Numerical,
     &           DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &           iCoSet,lOld,rHidden,nFix,nQQ,iRef,Redundant,nqInt,
     &           MaxItr,nWndw)
*
      nPrint(30) = nPrint(30)-1
*
      Call Put_dArray('BMtrx',Work(ipB),3*nsAtom*nQQ)
      Call Put_iScalar('No of Internal coordinates',nQQ)
*
*     Too many constraints?
*
      If (nLambda.gt.nQQ) Then
         Call WarningMessage(2,'Error in RlxCtl')
         Write (Lu,*)
         Write (Lu,*) '********************************************'
         Write (Lu,*) ' ERROR: nLambda.gt.nQQ'
         Write (Lu,*) ' nLambda=',nLambda
         Write (Lu,*) ' nQQ=',nQQ
         Write (Lu,*) ' There are more constraints than coordinates'
         Write (Lu,*) '********************************************'
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Reset_ThrGrd(nsAtom,nDimBC,Work(ipCM),nSym,iOper,Smmtrc,
     &                  Degen,Iter,Work(ipCx),mTtAtm,iWork(ipANr),
     &                  DDV_Schlegel,iOptC,rHidden,ThrGrd)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Compute the norm of the Cartesian gradient.
*
      ipOff = (iter-1)*3*nsAtom + ipGx
      Call G_Nrm(Work(ipOff),nsAtom,nQQ,Work(ipGNrm),iter,
     &           Work(ipdqInt),Degen,mIntEff)
      If (nPrint(116).ge.6) Call ListU(Lu,Lbl,Work(ipdqInt),nQQ,iter)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Accumulate gradient for complete or partial numerical
*     differentiation of the Hessian.
*
      If (lRowH.and.iter.lt.NmIter) Then
*
*----------------------------------------------------------------------*
*        I) Update geometry for selected numerical differentiation.    *
*----------------------------------------------------------------------*
*
         Call Freq1(iter,nQQ,nRowH,mRowH,Delta/2.5d0,Work(ipShf),
     &              Work(ipqInt))
         UpMeth='RowH  '
      Else If (lNmHss.and.iter.lt.NmIter) Then
*
*----------------------------------------------------------------------*
*        II) Update geometry for full numerical differentiation.       *
*----------------------------------------------------------------------*
*
         Call Freq2(iter,Work(ipdqInt),Work(ipShf),nQQ,Delta,Stop,
     &              Work(ipqInt),.Not.User_Def,nsAtom,Work(ipCM))
         UpMeth='NumHss'
      Else
         Go To 777
      End If
*
      Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nQQ,
     &            Work(ipdqInt+(iter-1)*nQQ),
     &            Work(ipShf+(iter-1)*nQQ),Lbl)
      iNeg(1)=-99
      iNeg(2)=-99
      HUpMet='None  '
      Stop = .False.
      nPrint(116)=nPrint(116)-3
      nPrint( 52)=nPrint( 52)-1  ! Status
      nPrint( 53)=nPrint( 53)-1
      nPrint( 54)=nPrint( 54)-1
      Write (6,*) ' Accumulate the gradient for selected '//
     &            'numerical differentiation.'
      Ed=Zero
      Step_Trunc=' '
         Go To 666
*
 777  Continue
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute updated geometry in Internal coordinates
*
      Step_Trunc=' '
      ed=zero
      If (lRowH.or.lNmHss) kIter = iter - (NmIter-1)
*
*     Update geometry
*
      Call Update_sl(Iter,MaxItr,NmIter,iInt,nFix,nQQ,Work(ipqInt),
     &               Work(ipShf),Work(ipdqInt),iOptC,Beta,
     &               Lbl,Work(ipGNrm),Work(ipEner),UpMeth,
     &               ed,Line_Search,Step_Trunc,nLambda,iRow_c,nsAtom,
     &               AtomLbl,nSym,iOper,mxdc,jStab,nStab,Work(ipB),
     &               Smmtrc,nDimBC,Work(ipL),ipCx,GrdMax,
     &               StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &               Labels,nLabels,FindTS,TSConstraints,nRowH,
     &               nWndw,Mode,ipMF,
     &               iOptH,HUpMet,kIter,GNrm_Threshold,
     &               IRC,Work(ipCM),HrmFrq_Show,
     &               CnstWght,Curvilinear,Redundant,Degen)
*
 666  Continue
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Transform the new internal coordinates to Cartesians
*
      Call GetMem(' DCF  ', 'Allo','Real',ipDCF, 3*nsAtom)
      Call GetMem(' dss  ', 'Allo','Real',ipdss, nQQ)
      Call GetMem(' qTemp', 'Allo','Real',ipTmp, nQQ)
      PrQ=.False.
      Call NewCar(Iter,nBVec,iRow,nsAtom,nDimBC,nQQ,Work(ipCoor),
     &            ipB,Work(ipCM),Lbl,Work(ipShf),ipqInt,
     &            ipdqInt,Work(ipDCF),Work(ipdss),Work(ipTmp),
     &            Stop,AtomLbl,iOper,nSym,iSym,Smmtrc,Degen,
     &            Work(ipGx),Work(ipCx),mTtAtm,iWork(ipANr),iOptH,
     &            User_Def,nStab,jStab,Curvilinear,Numerical,
     &            DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &            iCoSet,rHidden,ipRef,Redundant,nqInt,MaxItr)
      Call GetMem(' qTemp', 'Free','Real',ipTmp, nQQ)
      Call GetMem(' dss  ', 'Free','Real',ipdss, nQQ)
      Call GetMem(' DCF  ', 'Free','Real',ipDCF, 3*nsAtom)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*-----If this is a ESPF QM/MM job, the link atom coordinates are updated
*
      Do_ESPF = .False.
      Call DecideOnESPF(Do_ESPF)
      If (Do_ESPF) Then
       Call LA_Morok(nsAtom,ipCoor,2)
       call dcopy_(3*nsAtom,Work(ipCoor),1,Work(ipCx+3*nsAtom*Iter),1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Adjust some print levels
*
      If ((lNmHss.or.lRowH).and.iter.eq.NmIter) Then
*
*        If only frequencies no more output
*
         nPrint(21) = 5 ! Hessian already printed calling Update_sl
         If (kIter.gt.MxItr) Then
            Just_Frequencies=.True.
            nPrint(116)=nPrint(116)-3
            nPrint( 52)=nPrint( 52)-1
            nPrint( 53)=nPrint( 53)-1
            nPrint( 54)=nPrint( 54)-1
         End If
      End If
*
*     Fix correct reference structure in case of Numerical Hessian
*     optimization.
*
      If ((lNmHss.or.lRowH).and.kIter.eq.1) Then
         ip_1=ipCx
         ip_x=ipCx+(iter-1)*3*nsAtom
         call dcopy_(3*nsAtom,Work(ip_1),1,Work(ip_x),1)
      End If
*
*---- Print statistics and check on convergence
*
      GoOn = (lNmHss.and.iter.lt.NmIter).OR.(lRowH.and.iter.lt.NmIter)
      TSReg = iAnd(iOptC,8192).eq.8192
      Call Convrg(iter,kIter,nQQ,Work(ipqInt),Work(ipShf),
     &            Work(ipdqInt),Lbl,Work(ipGNrm),
     &            Work(ipEner),Stat,MaxItr,Stop,iStop,ThrCons,
     &            ThrEne,ThrGrd,MxItr,UpMeth,HUpMet,mIntEff,Baker,
     &            Work(ipCx),Work(ipGx),nsAtom,mTtAtm,iOper,nSym,ed,
     &            iNeg,GoOn,Step_Trunc,GrdMax,StpMax,GrdLbl,StpLbl,
     &            Analytic_Hessian,rMEP,MEP,nMEP,
     &            (lNmHss.or.lRowH).and.iter.le.NmIter,
     &            Just_Frequencies,FindTS,ipCoor,eMEPTest,nLambda,
     &            TSReg)
      Call Free_Work(ipShf)
*
************************************************************************
*                                                                      *
*                           EPILOGUE                                   *
*                                                                      *
************************************************************************
*
*-----Write information to files
*
      Call DstInf(iStop,Just_Frequencies,
     &            (lNmHss.or.lRowH) .and.iter.le.NmIter)
      If (lCtoF) Call Def_CtoF(.True.,Work(ipCM),
     &         nsAtom,AtomLbl,Work(ipCoor),nSym,iOper,jStab,nStab)
      If (.Not.User_Def .and.
     &   ((lNmHss.and.iter.ge.NmIter).or..Not.lNmHss)) Call cp_SpcInt
*
*-----After a numerical frequencies calculation, restore the original
*     runfile, but save the useful data (gradient and Hessian)
*
      If (lNmHss.and.iter.ge.NmIter) Then
         Call f_Inquire('RUNBACK',Found)
         If (Found) Then
*           Read data
            Call Get_Grad(ipGB,nGB)
            Call Qpg_dArray('Hss_X',Found,nHX)
            Call GetMem('HssX','Allo','Real',ipHX,nHX)
            Call Get_dArray('Hss_X',Work(ipHX),nHX)
            Call Qpg_dArray('Hss_Q',Found,nHQ)
            Call GetMem('HssQ','Allo','Real',ipHQ,nHQ)
            Call Get_dArray('Hss_Q',Work(ipHQ),nHQ)
            Call Qpg_dArray('KtB',Found,nKtB)
            Call GetMem('KtB','Allo','Real',ipKtB,nKtB)
            Call Get_dArray('KtB',Work(ipKtB),nKtB)
            Call Get_iScalar('No of Internal coordinates',nIntCoor)
*           Write data in backup file
            Call NameRun('RUNBACK')
            Call Put_Grad(Work(ipGB),nGB)
            Call Put_dArray('Hss_X',Work(ipHX),nHX)
            Call Put_dArray('Hss_Q',Work(ipHQ),nHQ)
            Call Put_dArray('Hss_upd',Work(ip_Dummy),0)
            Call Put_dArray('Hess',Work(ipHQ),nHQ)
            Call Put_dArray('KtB',Work(ipKtB),nKtB)
            Call Put_iScalar('No of Internal coordinates',nIntCoor)
*           Pretend the Hessian is analytical
            nHX2=Int(Sqrt(Dble(nHX)))
            iOff=0
            Do i=1,nHX2
               Do j=1,i
                  Work(ipHx+iOff)=Work(ipHX+(i-1)*nHX2+(j-1))
                  iOff=iOff+1
               End Do
            End Do
            Call Put_AnalHess(Work(ipHX),iOff)
            Call NameRun('#Pop')
            Call GetMem('Grad','Free','Real',ipGB,nGB)
            Call GetMem('HssX','Free','Real',ipHX,nHX)
            Call GetMem('HssQ','Free','Real',ipHQ,nHQ)
            Call GetMem('KtB','Free','Real',ipKtB,nKtB)
*
*           Restore and remove the backup runfile
*
            Call fCopy('RUNBACK','RUNFILE',iErr)
            If (iErr.ne.0) Call Abend
            If (AixRm('RUNBACK').ne.0) Call Abend
         End If
      End If
*
*-----Remove the GRADS file
*
      Call f_Inquire('GRADS',Found)
      If (Found) Then
         If (AixRm('GRADS').ne.0) Call Abend()
      End If
*
      Call Chkpnt_update()
      Call Chkpnt_close()
*                                                                      *
************************************************************************
*                                                                      *
*-----Deallocate memory
*
      Call GetMem(' B ',    'Free','Real',ipB,   (nsAtom*3)**2)
 999  Continue
*
      If (ip_B.ne.ip_Dummy) Call Free_Work(ip_B)
      If (ip_dB.ne.ip_Dummy) Call Free_Work(ip_dB)
      If (ip_iB.ne.ip_iDummy) Call Free_iWork(ip_iB)
      If (ip_idB.ne.ip_iDummy) Call Free_iWork(ip_idB)
      If (ip_nqB.ne.ip_iDummy) Call Free_iWork(ip_nqB)
*
      If (ipNADC.ne.ip_Dummy) Call Free_Work(ipNADC)
      If (Ref_Geom) Call Free_Work(ipRef)
      If (Ref_Grad) Call Free_Work(ipGradRef)
      If (lRP)      Call Free_Work(ipR12)
      If (ipqInt.ne.ip_Dummy) Then
          Call GetMem('dqInt', 'Free','Real',ipdqInt, nqInt)
          Call GetMem('qInt', 'Free','Real',ipqInt, nqInt)
      End If
      Call GetMem('Relax', 'Free','Real',ipRlx, Lngth)
      Call GetMem('Grad',  'Free','Real',ipGrd, 3*nsAtom)
      Call GetMem('Coord', 'Free','Real',ipCoor,3*nsAtom)
      Call GetMem('Anr',   'Free','Inte',ipANr, nsAtom)
      Call GetMem('Charge','Free','Real',ipCM,  nsAtom)
*     The weights array length is actually the total number of atoms,
*     not just symmetry-unique, but that doesn't matter for deallocation
      Call GetMem('Weights','Free','Real',ipWeights,nsAtom)
*
*-----Terminate the calculations.
*
      Call QExit('RlxCtl')
*
      Return
      End
