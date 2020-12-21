**********************************************************************
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
      Use Slapaf_Info, only: Cx, Coor, Shift, GNrm, BMx, mRowH,
     &                       Free_Slapaf, qInt, dqInt, Lbl
      use Slapaf_Parameters, only: HUpMet, User_Def, iOptC, UpMeth,
     &                             HSet, BSet, PrQ, Numerical, iNeg,
     &                             E_Delta
      Implicit Real*8 (a-h,o-z)
************************************************************************
*     Program for determination of the new molecular geometry          *
************************************************************************
#include "info_slapaf.fh"
#include "real.fh"
#include "nadc.fh"
#include "weighting.fh"
#include "db.fh"
#include "stdalloc.fh"
#include "print.fh"
      Logical GoOn, TSReg, Do_ESPF, Just_Frequencies, Found, Error
      Character(LEN=1) Step_trunc
      Integer, External:: AixRm
      Integer nGB
      Real*8 rDum(1)
      Real*8, Allocatable:: GB(:), HX(:), HQ(:), KtB(:)
*                                                                      *
************************************************************************
*                                                                      *
      Lu=6
      Just_Frequencies=.False.
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*
      Call RdCtl_Slapaf(LuSpool,.False.)
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
         Call Free_Slapaf()
         Return
       Else If (isFalcon) Then
         iStop=1
         Call Free_Slapaf()
         Return
      End if
*                                                                      *
************************************************************************
*                                                                      *
      PrQ= .Not.Request_Alaska
*                                                                      *
************************************************************************
*                                                                      *
      If (lCtoF .AND. PrQ) Call Def_CtoF(.False.)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the Wilson B-matrices, which describe the transformations
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
      ! Numerical only for some rows
      If (Allocated(mRowH))  NmIter=SIZE(mRowH)+1
      If (lNmHss) NmIter=2*mInt+1    ! Full numerical
      If (Cubic)  NmIter=2*mInt**2+1 ! Full cubic
*
      If (lTherm .and. iter.EQ.1) then
         Call Put_dArray('Initial Coordinates',Coor,SIZE(Coor))
      EndIf
*
*---- Fix the definition of internal during numerical differentiation
      If (lNmHss.and.iter.lt.NmIter.and.iter.ne.1) nPrint(122)=5
*
*---- Do not overwrite numerical Hessian
      If ((lNmHss.or.Allocated(mRowH))
     &    .and.(iter.gt.NmIter.or.iter.lt.NmIter)) HSet = .False.
*
*---- Set logical to indicate status during numerical differentiation
      Numerical = lNmHss .and.iter.le.NmIter .and.iter.ne.1
*
      If (Numerical) nWndw=NmIter
      iRef=0
      Call BMtrx(SIZE(Coor,2),mInt,Coor,iter,mTtAtm,nQQ,iRef,nWndw)
*
      nPrint(30) = nPrint(30)-1
*
      Call Put_dArray('BMtrx',BMx,SIZE(Coor)*nQQ)
      Call Put_iScalar('No of Internal coordinates',nQQ)
*                                                                      *
************************************************************************
*                                                                      *
      Call Reset_ThrGrd(Iter,mTtAtm,ThrGrd)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Compute the norm of the Cartesian gradient.
*
      Call G_Nrm(nQQ,GNrm,iter,dqInt,mIntEff)
      If (nPrint(116).ge.6) Call ListU(Lu,Lbl,dqInt,nQQ,iter)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Accumulate gradient for complete or partial numerical
*     differentiation of the Hessian.
*
      If ((Allocated(mRowH).and.iter.lt.NmIter) .or.
     &        (lNmHss.and.iter.lt.NmIter)) Then

         If (Allocated(mRowH).and.iter.lt.NmIter) Then
*
*----------------------------------------------------------------------*
*        I) Update geometry for selected numerical differentiation.    *
*----------------------------------------------------------------------*
*
            Call Freq1(iter,nQQ,Delta/2.5d0,qInt)
            UpMeth='RowH  '
         Else
*
*----------------------------------------------------------------------*
*        II) Update geometry for full numerical differentiation.       *
*----------------------------------------------------------------------*
*
            Call Freq2(iter,dqInt,nQQ,Delta,Stop,qInt)
            UpMeth='NumHss'
         End If
*
         Call MxLbls(nQQ,dqInt(:,iter),Shift(:,iter),Lbl)
         iNeg(:)=-99
         HUpMet='None  '
         Stop = .False.
         nPrint(116)=nPrint(116)-3
         nPrint( 52)=nPrint( 52)-1  ! Status
         nPrint( 53)=nPrint( 53)-1
         nPrint( 54)=nPrint( 54)-1
         Write (6,*) ' Accumulate the gradient for selected '//
     &               'numerical differentiation.'
         Write (6,'(1x,i5,1x,a,1x,i5)') iter,'of',NmIter
         E_Delta=Zero
         Step_Trunc=' '
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*--------Compute updated geometry in Internal coordinates
*
         Step_Trunc=' '
         E_Delta=zero
         If (Allocated(mRowH).or.lNmHss) kIter = iter - (NmIter-1)
*define UNIT_MM
#ifdef UNIT_MM
         Call Init_UpdMask(nInter)
#endif
*
*        Update geometry
*
         If (Kriging .and. Iter.ge.nspAI) Then
            Call Update_Kriging(Iter,nQQ,Step_Trunc,nWndw)
         Else
            Call Update_sl(Iter,NmIter,nQQ,Step_Trunc,nWndw,kIter)
         End If
*
#ifdef UNIT_MM
         Call Free_UpdMask()
#endif
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Transform the new internal coordinates to Cartesians
*     (if not already done by Kriging)
*
      If (Kriging .and. Iter.ge.nspAI) Then
         Coor(:,:) = Cx(:,:,Iter+1)
      Else
         PrQ=.False.
         Error=.False.
         iRef=0
         Call NewCar(Iter,Size(Coor,2),nQQ,Coor,iSym,mTtAtm,iRef,Error)
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
       Call LA_Morok(SIZE(Coor,2),Coor,2)
       Cx(:,:,Iter+1) = Coor(:,:)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Adjust some print levels
*
      If ((lNmHss.or.Allocated(mRowH)).and.iter.eq.NmIter) Then
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
      If ((lNmHss.or.Allocated(mRowH)).and.kIter.eq.1) Then
         Cx(:,:,iter) = Cx(:,:,1)
      End If
*
*---- Print statistics and check on convergence
*
      GoOn = (lNmHss.and.iter.lt.NmIter) .OR.
     &       (Allocated(mRowH).and.iter.lt.NmIter)
      TSReg = iAnd(iOptC,8192).eq.8192
      Numerical=(lNmHss.or.Allocated(mRowH)).and.iter.le.NmIter
      Call Convrg(iter,kIter,nQQ,Stop,iStop,ThrCons,MxItr,mIntEff,Baker,
     &            mTtAtm,GoOn,Step_Trunc,rMEP,MEP,nMEP,Just_Frequencies,
     &            eMEPTest,TSReg,ThrMEP)
*
************************************************************************
*                                                                      *
*                           EPILOGUE                                   *
*                                                                      *
************************************************************************
*
*-----Write information to files
*
      Numerical = (lNmHss.or.Allocated(mRowH)) .and. iter.le.NmIter
      Call DstInf(iStop,Just_Frequencies)
      If (lCtoF) Call Def_CtoF(.True.)
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
            nGB=SIZE(Coor)
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
#ifdef _DEBUGPRINT_
            Call TriPrt('AnalHess',' ',HX,nHX2)
#endif

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
      Call Free_Slapaf()
*
*-----Terminate the calculations.
*
      Return
      End
