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
      Use kriging_mod, only: Kriging, nspAI
      Use Slapaf_Info, only: Cx, Coor, Shift, GNrm,
     &                       Free_Slapaf, qInt, dqInt
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
#include "stdalloc.fh"
      Logical Numerical, GoOn, PrQ, TSReg,
     &        Do_ESPF, Just_Frequencies, Found, Error
      Character*8 GrdLbl, StpLbl, Labels(nLabels), Lbl(nLbl)
      Character*1 Step_trunc
      Integer AixRm, iNeg(2)
      Integer nGB
      Real*8 rDum(1)
      Real*8, Allocatable:: GB(:), HX(:), HQ(:), KtB(:)
*
      Lu=6
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
      If (lCtoF .AND. PrQ) Call Def_CtoF(.False.,nsAtom,AtomLbl,
     &                                   Coor,jStab,nStab)
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
         Call Put_dArray('Initial Coordinates',Coor,3*nsAtom)
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
      iRef=0
      Call BMtrx(iRow,nBVec,ipB,nsAtom,mInt,Lbl,
     &           Coor,nDimBC,AtomLbl,
     &           Smmtrc,Degen,BSet,HSet,iter,
     &           mTtAtm,iOptH,
     &           User_Def,nStab,jStab,Curvilinear,Numerical,
     &           DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &           iCoSet,lOld,rHidden,nFix,nQQ,iRef,Redundant,
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
      Call Reset_ThrGrd(nsAtom,nDimBC,Smmtrc,Degen,Iter,
     &                  mTtAtm,DDV_Schlegel,iOptC,rHidden,
     &                  ThrGrd)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Compute the norm of the Cartesian gradient.
*
      Call G_Nrm(nsAtom,nQQ,GNrm,iter,dqInt,Degen,mIntEff)
      If (nPrint(116).ge.6) Call ListU(Lu,Lbl,dqInt,nQQ,iter)
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
         Call Freq1(iter,nQQ,nRowH,mRowH,Delta/2.5d0,qInt)
         UpMeth='RowH  '
      Else If (lNmHss.and.iter.lt.NmIter) Then
*
*----------------------------------------------------------------------*
*        II) Update geometry for full numerical differentiation.       *
*----------------------------------------------------------------------*
*
         Call Freq2(iter,dqInt,nQQ,Delta,Stop,qInt)
         UpMeth='NumHss'
      Else
         Go To 777
      End If
*
      Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nQQ,dqInt(:,iter),
     &            Shift(:,iter),Lbl)
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
      Write (6,'(1x,i5,1x,a,1x,i5)') iter,'of',NmIter
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
*define UNIT_MM
#ifdef UNIT_MM
      Call Init_UpdMask(Curvilinear, Redundant, nsAtom, nInter)
#endif
*
*     Update geometry
*
      If (Kriging .and. Iter.ge.nspAI) Then
         Call Update_Kriging(
     &               Iter,iInt,nFix,nQQ,
     &               iOptC,Beta,Beta_Disp,
     &               Lbl,UpMeth,
     &               ed,Line_Search,Step_Trunc,nLambda,iRow_c,nsAtom,
     &               AtomLbl,mxdc,jStab,nStab,Work(ipB),
     &               Smmtrc,nDimBC,
     &               GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &               Labels,nLabels,FindTS,TSConstraints,nRowH,
     &               nWndw,Mode,
     &               iOptH,HUpMet,GNrm_Threshold,
     &               IRC,HrmFrq_Show,
     &               CnstWght,Curvilinear,Degen,ThrEne,ThrGrd,iRow)
      Else
         Call Update_sl(
     &               Iter,NmIter,iInt,nFix,nQQ,
     &               iOptC,Beta,Beta_Disp,
     &               Lbl,UpMeth,
     &               ed,Line_Search,Step_Trunc,nLambda,iRow_c,nsAtom,
     &               AtomLbl,mxdc,jStab,nStab,Work(ipB),
     &               Smmtrc,nDimBC,GrdMax,
     &               StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &               Labels,nLabels,FindTS,TSConstraints,nRowH,
     &               nWndw,Mode,
     &               iOptH,HUpMet,kIter,GNrm_Threshold,
     &               IRC,HrmFrq_Show,
     &               CnstWght,Curvilinear,Degen)
      End If
*
#ifdef UNIT_MM
      Call Free_UpdMask()
#endif
*
 666  Continue
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Transform the new internal coordinates to Cartesians
*     (if not already done by Kriging)
*
      If (Kriging .and. Iter.ge.nspAI) Then
         Call dCopy_(3*nsAtom,Cx(1,1,Iter+1),1,Coor,1)
      Else
         PrQ=.False.
         Error=.False.
         iRef=0
         Call NewCar(Iter,nBVec,iRow,nsAtom,nDimBC,nQQ,Coor,
     &               ipB,Lbl,AtomLbl,iSym,Smmtrc,Degen,
     &               mTtAtm,iOptH,
     &               User_Def,nStab,jStab,Curvilinear,Numerical,
     &               DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &               iCoSet,rHidden,Redundant,MaxItr,iRef,Error)
      End If
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
       Call LA_Morok(nsAtom,Coor,2)
       call dcopy_(3*nsAtom,Coor,1,Cx(1,1,Iter+1),1)
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
         call dcopy_(3*nsAtom,Cx(1,1,1),1,Cx(1,1,iter),1)
      End If
*
*---- Print statistics and check on convergence
*
      GoOn = (lNmHss.and.iter.lt.NmIter).OR.(lRowH.and.iter.lt.NmIter)
      TSReg = iAnd(iOptC,8192).eq.8192
      Call Convrg(iter,kIter,nQQ,Lbl,MaxItr,Stop,iStop,ThrCons,
     &            ThrEne,ThrGrd,MxItr,UpMeth,HUpMet,mIntEff,Baker,
     &            nsAtom,mTtAtm,ed,
     &            iNeg,GoOn,Step_Trunc,GrdMax,StpMax,GrdLbl,StpLbl,
     &            Analytic_Hessian,rMEP,MEP,nMEP,
     &            (lNmHss.or.lRowH).and.iter.le.NmIter,
     &            Just_Frequencies,FindTS,eMEPTest,nLambda,
     &            TSReg,ThrMEP)
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
      If (lCtoF) Call Def_CtoF(.True.,nsAtom,AtomLbl,Coor,jStab,nStab)
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
            nGB=3*nsAtom
            Call mma_allocate(GB,nGB,Label='GB')
            Call Get_Grad(GB,nGB)

            Call Qpg_dArray('Hss_X',Found,nHX)
            Call mma_allocate(HX,nHX,Label='HX')
            Call Get_dArray('Hss_X',HX,nHX)

            Call Qpg_dArray('Hss_Q',Found,nHQ)
            Call mma_allocate(HQ,nHQ,Label='HQ')
            Call Get_dArray('Hss_Q',HQ,nHQ)

            Call Qpg_dArray('KtB',Found,nKtB)
            Call mma_allocate(KtB,nKtB,Label='KtB')
            Call Get_dArray('KtB',KtB,nKtB)

            Call Get_iScalar('No of Internal coordinates',nIntCoor)
*           Write data in backup file
            Call NameRun('RUNBACK')
            Call Put_Grad(GB,nGB)
            Call Put_dArray('Hss_X',HX,nHX)
            Call Put_dArray('Hss_Q',HQ,nHQ)
            Call Put_dArray('Hss_upd',rdum,0)
            Call Put_dArray('Hess',HQ,nHQ)
            Call Put_dArray('KtB',KtB,nKtB)
            Call Put_iScalar('No of Internal coordinates',nIntCoor)
*           Pretend the Hessian is analytical
            nHX2=Int(Sqrt(Dble(nHX)))
            iOff=0
            Do i=1,nHX2
               Do j=1,i
                  iOff=iOff+1
                  HX(iOff)=HX((i-1)*nHX2+j)
               End Do
            End Do

            Call Put_AnalHess(HX,iOff)
            Call NameRun('#Pop')

            Call mma_deallocate(GB)
            Call mma_deallocate(HX)
            Call mma_deallocate(HQ)
            Call mma_deallocate(KtB)
*
*           Restore and remove the backup runfile
*
            Call fCopy('RUNBACK','RUNFILE',iErr)
            If (iErr.ne.0) Call Abend()
            If (AixRm('RUNBACK').ne.0) Call Abend()
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
      If (Ref_Grad) Call Free_Work(ipGradRef)
      If (lRP)      Call Free_Work(ipR12)
      Call GetMem('Relax', 'Free','Real',ipRlx, Lngth)
      Call Free_Slapaf()
*
*-----Terminate the calculations.
*
*
      Return
      End
