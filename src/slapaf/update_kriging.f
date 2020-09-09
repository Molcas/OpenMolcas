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
* Copyright (C) 2019,2020, Roland Lindh                                *
*               2020, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Update_kriging(
     &                     iter,MaxItr,NmIter,iInt,nFix,nInter,qInt,
     &                     Shift,
     &                     Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,
     &                     Energy,UpMeth,ed,Line_Search,Step_Trunc,
     &                     nLambda,iRow_c,nsAtom,AtomLbl,
     &                     mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                     rLambda,Cx,Gx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,MF,
     &                     iOptH,HUpMet,kIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Degen,ThrEne,ThrGrd)
************************************************************************
*                                                                      *
*     Object: to update coordinates                                    *
*                                                                      *
*    (see update_sl)                                                   *
************************************************************************
      Use kriging_mod, only: miAI, meAI
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 qInt(nInter,MaxItr), Shift(nInter,MaxItr),
     &       Grad(nInter,MaxItr), GNrm(MaxItr), Energy(MaxItr),
     &       BMx(3*nsAtom,3*nsAtom), rLambda(nLambda,MaxItr),
     &       dMass(nsAtom), Degen(3*nsAtom), dEner,
     &       Cx(3*nsAtom,Iter+1), Gx(3*nsAtom,Iter), MF(3*nsAtom)
      Integer jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),
     &        FindTS, TSC, HrmFrq_Show, Curvilinear,
     &        First_MicroIteration, Error
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6
      Character GrdLbl_Save*8
      Real*8 Dummy(1)
      Real*8, Allocatable:: Hessian(:,:), Temp(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Different hardwired kriging options
*
*#define _DEBUG_
*#define _OVERSHOOT_
*                                                                      *
************************************************************************
*                                                                      *
      Logical Kriging_Hessian, Not_Converged, Force_RS
#ifdef _OVERSHOOT_
      Real*8, Allocatable:: Step_k(:,:)
#endif
*
      iRout=153
      iPrint=nPrint(iRout)
      Call QEnter('Update')
*
#ifdef _DEBUG_
      Call RecPrt('Update_Kriging: qInt',' ',qInt,nInter,Iter)
      Call RecPrt('Update_Kriging: Shift',' ',Shift,nInter,Iter-1)
      Call RecPrt('Update_Kriging: GNrm',' ',GNrm,Iter,1)
#endif
*
      Kriging_Hessian =.TRUE.
      Force_RS=.FALSE.
      iOpt_RS=1   ! Activate restricted variance.
      iterAI=iter
      dEner=meAI
      nRaw=Min(iter,nWndw/2)
*     nRaw=1 ! Force HMF Hessian
      iFirst = iter - nRaw + 1
      iterK=0
      dqdq=Zero
      qBeta=Beta
      qBeta_Disp=Beta_Disp
      GrdMax_Save=GrdMax
      GrdLbl_Save=GrdLbl
#ifdef _DEBUG_
      Call RecPrt('qInt(0)',  ' ',qInt(1,iFirst),nInter,nRaw)
      Call RecPrt('Energy(0)',' ',Energy(iFirst),1,nRaw)
      Call RecPrt('Grad(0)',  ' ',Grad(1,iFirst),nInter,nRaw)
      Call RecPrt('Shift',  ' ',Shift(1,iFirst),nInter,nRaw)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up the HMF Hessian
*
      Call mma_Allocate(Hessian,nInter,nInter,Label='Hessian')
      Call Mk_Hss_Q()
      Call Get_dArray('Hss_Q',Hessian,nInter**2)
*     Call RecPrt('HMF Hessian',' ',Hessian,nInter,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*     Select the data points to pass to the GEK routine. Remember to
*     change the sign of the gradients.
*                                                                      *
************************************************************************
*                                                                      *
*     Note that we could have some kind of sorting here if we like!
*
************************************************************************
*                                                                      *
*     Pass the sample points to the GEK procedure and deallocate the
*     memory -- to be reused for another purpose later.
*
      Call DScal_(nInter*nRaw,-One,Grad(1,iFirst),1)
      Call Setup_Kriging(nRaw,nInter,qInt(1,iFirst),
     &                               Grad(1,iFirst),
     &                               Energy(iFirst),Hessian)
      Call DScal_(nInter*nRaw,-One,Grad(1,iFirst),1)
*                                                                      *
************************************************************************
*                                                                      *
*     Save initial gradient
*     This is for the case the gradient is already converged, we want
*     the micro iterations to still reduce the gradient
*
      FAbs_ini=Sqrt(DDot_(nInter,Grad(1,iter),1,
     &                           Grad(1,iter),1)/DBLE(nInter))
      GrdMx_ini=Zero
      Do iInter = 1, nInter
         GrdMx_ini=Max(GrdMx_ini,Abs(Grad(iInter,iterAI)))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Find the fraction value to be used to define the restricted
*     variance threshold.
*
*
      Beta_Disp_Min=1.0D-10
*     Let the accepted variance be set as a fraction of the
*     largest component in the gradient.
      tmp=0.0D0
      Do i = 1, 3*nsAtom
         tmp = Max(tmp,Abs(Gx(i,iter)))
      End Do
*
*     Temporary code until we have figured out this for constrained
*     optimizations.
*
      Beta_Disp_=Max(Beta_Disp_Min,tmp*Beta_Disp)
      Beta_=Min(1.0D3*GNrm(iter),Beta)
*
#ifdef _RS_RFO_
*     Switch over to RS-RFO once the gradient is low.
*
      tmp=99.0D0
      Do j = 1, iter
         tmp0=0.0D0
         Do i = 1, 3*nsAtom
            tmp0 = Max(tmp0,Abs(Gx(i,j)))
         End Do
         tmp=Min(tmp,tmp0)
      End Do

      If (tmp.lt.4.0D-4) Then
         iOpt_RS=0
         Beta_=0.03D0
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Start the Kriging loop.
*
      Not_Converged = .True.
      Step_Trunc='N'  ! not defined
      Do While (Not_Converged)
*                                                                      *
************************************************************************
*                                                                      *
         kIter_=iterAI
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) 'Do iterAI: ',iterAI
#endif
*
*        Compute the Kriging Hessian
*
         If (Kriging_Hessian) Then
            Call Hessian_Kriging_Layer(qInt(1,iterAI),Hessian,nInter)
            Call Put_dArray('Hss_Q',Hessian,nInter**2)
*           Make fixhess.f treat the Hessian as if it was analytic.
            Call Put_iScalar('HessIter',iterAI)
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the updated structure.
*
         First_MicroIteration=iterAI.eq.iter
         nWndw_=nWndw/2 + (iterAI-iter)
         Call Update_sl_(iterAI,iInt,nFix,nInter,qInt,Shift,Grad,iOptC,
     &                   Beta_,Beta_Disp_,Lbl,GNrm,Energy,
     &                   UpMeth,ed,Line_Search,Step_Trunc,nLambda,
     &                   iRow_c,nsAtom,AtomLbl,mxdc,jStab,
     &                   nStab,BMx,Smmtrc,nDimBC,rLambda,Cx,
     &                   GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &                   Labels,nLabels,FindTS,TSC,nRowH,nWndw_,Mode,
     &                   MF,iOptH,HUpMet,kIter_,GNrm_Threshold,IRC,
     &                   dMass,HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                   Kriging_Hessian,qBeta,iOpt_RS,
     &                   First_MicroIteration,iter,qBeta_Disp)
#ifdef _DEBUG_
         Write (6,*) 'Update_kriging: Step_Trunc',Step_Trunc
         Call RecPrt('New Coord',' ',qInt,nInter,iterAI+1)
         Call RecPrt('Grad',' ',Grad,nInter,iterAI)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        In case of a constrained optimization GrdMax and GNrm has been
*        updated to reflect the contribution from the constraint.
*        GrdMax, however, is a scalar and we need to save it such that
*        it has the correct value on exit. That is, the value of itera-
*        tion iter.
*
         If (iterAI.eq.iter) Then
            GrdMax_Save=GrdMax
            GrdLbl_Save=GrdLbl
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        During the micro iterations we need to set GNrm to some
*        reasonable value. Normally it is set by the force routine or
*        the con_opt routine.
*
         If (nLambda.eq.0 .and. iterAI.gt.iter) Then
            GNrm(iterAI)=Sqrt(DDot_(nInter,Grad(1,iterAI),1,
     &                                    Grad(1,iterAI),1))
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the step length from the last ab initio point
*        to the most recent kriging point.
*
         dqdq=Zero
         Do iInter=1,nInter
            dqdq=(qInt(iInter,iterAI+1)
     &           -qInt(iInter,iFirst+nRaw-1))**2
         End Do
         If (iterK.eq.0.and.Step_trunc.eq.'*') dqdq=qBeta**2
*                                                                      *
************************************************************************
*                                                                      *
*        Transform the new internal coordinates to Cartesians
*
         Error=(iterK.ge.1)
         Call NewCar_Kriging(iterAI,iRow,nsAtom,nDimBC,nInter,BMx,dMass,
     &                       Lbl,Shift,qInt,Grad,AtomLbl,
     &                       Cx,.True.,iter,Error)
         If (Error) Then
            Step_Trunc='@'
            Exit
         End If
*
*        Compute the energy and gradient according to the
*        surrogate model for the new coordinates.
*
         Call Energy_Kriging_layer(qInt(1,iterAI+1),Energy(iterAI+1),
     &                             nInter)
         Call Dispersion_Kriging_Layer(qInt(1,iterAI+1),E_Disp,nInter)
         Call Gradient_Kriging_layer(qInt(1,iterAI+1),
     &                               Grad(1,iterAI+1),nInter)
         Call DScal_(nInter,-One,Grad(1,iterAI+1),1)
*                                                                      *
************************************************************************
*                                                                      *
         iterK  = iterK  + 1
         iterAI = iterAI + 1
         dEner = Energy(iterAI) - Energy(iterAI-1)
#ifdef _DEBUG_
         Call RecPrt('qInt(x):',' ',qInt,nInter,iterAI)
         Call RecPrt('Ener(x):',' ',Energy,1,iterAI)
         Call RecPrt('Grad(x):',' ',Grad,nInter,iterAI)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Check on convergence criteria.
*
         If (ThrEne.gt.Zero) Then
*           Convergence on energy criteria.
            Not_Converged = Abs(dEner).ge.ThrEne
            Not_Converged = Not_Converged .and. dqdq.lt.qBeta**2
         Else
*           Use standard convergence criteria
            FAbs=GNrm(iterAI-1)/SQRT(DBLE(nInter-nLambda))
            RMS =Sqrt(DDot_(nInter,Shift(1,iterAI-1),1,
     &                             Shift(1,iterAI-1),1)/DBLE(nInter))
            RMSMx=Zero
            Do iInter = 1, nInter
               RMSMx=Max(RMSMx,Abs(Shift(iInter,iterAI-1)))
            End Do
            GrdMx=GrdMax
#ifdef _DEBUG_
            Write (6,*)
            Write (6,*) 'iter=',iterAI-1
            Write (6,*) 'FAbs=',FAbs
            Write (6,*) 'GrdMx=',GrdMx
            Write (6,*) 'RMS=',RMS
            Write (6,*) 'RMSMx=',RMSMx
            Write (6,*) FAbs.gt.Min(ThrGrd,FAbs_ini)
            Write (6,*) GrdMx.gt.Min(ThrGrd*OneHalf,GrdMx_ini)
            Write (6,*) RMS.gt.ThrGrd*Four
            Write (6,*) RMSMx.gt.ThrGrd*Six
#endif
*           Ensure that the initial gradient is reduced,
*           except in the last micro iteration
            If (iterK.lt.miAI) Then
               Not_Converged = FAbs.gt.Min(ThrGrd,FAbs_ini)
               Not_Converged = Not_Converged .or.
     &                         GrdMx.gt.Min(ThrGrd*OneHalf,GrdMx_ini)
            Else
               Not_Converged = FAbs.gt.ThrGrd
               Not_Converged = Not_Converged .or.
     &                         GrdMx.gt.ThrGrd*OneHalf
            End If
            Not_Converged = Not_Converged .or.
     &                      RMS.gt.ThrGrd*Four
            Not_Converged = Not_Converged .or.
     &                      RMSMx.gt.ThrGrd*Six
            Not_Converged = Not_Converged .or.
     &                      Step_Trunc.ne.' '
         End If
         If (Step_Trunc.eq.'.') Step_Trunc=' '
*        Check total displacement from initial structure
         If (Force_RS) Then
            Call mma_allocate(Temp,3*nsAtom,1)
            Temp(:,1)=Cx(:,iterAI)-Cx(:,iter)
            RMS =Sqrt(DDot_(3*nsAtom,Temp,1,Temp,1)/DBLE(3*nsAtom))
            If (RMS.gt.(Three*Beta_)) Step_trunc='*'
            Call mma_deAllocate(Temp)
         End If
         If (Not_Converged.and.(Step_trunc.eq.' ').and.(iterK.ge.miAI))
     &      Step_trunc='#'
#ifdef _DEBUG_
         Write (6,*) 'Not_Converged=',Not_Converged
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        If the step restriction is invoked or there is some other
*        reason signaled by Step_Trunc, terminate anyhow.
*
         If (Step_trunc.ne.' ') Not_Converged=.False.
*
*        If RS rather than RV do not micro iterate
*
         If (iOpt_RS.eq.0) Not_Converged=.False.
*        Not_Converged=.False. ! Force single iteration scheme.
#ifdef _DEBUG_
         Write (6,*) 'Not_Converged=',Not_Converged
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     End of the micro iteration loop
*
      End Do  ! Do While (Not_Converged)
*
*     Change label of updating method
*
      UpMeth='RVO   '
      Write(UpMeth(4:6),'(I3)') iterK
*
#ifdef _OVERSHOOT_
*                                                                      *
************************************************************************
*                                                                      *
*     Attempt overshooting
*
      Call mma_allocate(Step_k,nInter,2)
      i=iFirst+nRaw-1
      Call dCopy_(nInter,qInt(1,iterAI),1,Step_k(1,2),1)
      Call dAXpY_(nInter,-One,qInt(1,i),1,Step_k(1,2),1)
*     Call RecPrt('q(i+1)-q(f)','',Step_k(1,2),nInter,1)
      If (i.gt.1) Then
         Call dCopy_(nInter,qInt(1,i),1,Step_k(1,1),1)
         Call dAXpY_(nInter,-One,qInt(1,i-1),1,Step_k(1,1),1)
*        Call RecPrt('q(f)-q(f-1)','',Step_k(1,1),nInter,1)
         dsds=dDot_(nInter,Step_k(1,1),1,Step_k(1,2),1)
         dsds=dsds/Sqrt(ddot_(nInter,Step_k(1,1),1,Step_k(1,1),1))
         dsds=dsds/Sqrt(ddot_(nInter,Step_k(1,2),1,Step_k(1,2),1))
         Write(6,*) 'dsds = ',dsds
      Else
         dsds=Zero
      End If
      dsds_min=0.9D0
      If ((dsds.gt.dsds_min).And.(Step_Trunc.eq.' ')) Then
        Do Max_OS=9,0,-1
         OS_Factor=Max_OS*((dsds-dsds_min)/(One-dsds_min))**4
         Call dCopy_(nInter,qInt(1,i),1,qInt(1,iterAI),1)
         Call dAXpY_(nInter,One+OS_Factor,Step_k(1,2),1,
     &                                    qInt(1,iterAI),1)
         Call Energy_Kriging_Layer(qInt(1,iterAI),OS_Energy,nInter)
         Call Dispersion_Kriging_Layer(qInt(1,iterAI),OS_Disp,nInter)
         Write(6,*) 'Max_OS=',Max_OS
         Write(6,*) OS_Disp,E_Disp,Beta_Disp_
         If ((OS_Disp.gt.E_Disp).And.(OS_Disp.lt.Beta_Disp_)) Then
            Call dAXpY_(nInter,OS_Factor,Step_k(1,2),1,
     &                                   Shift(1,iterAI-1),1)
            Call NewCar_Kriging(iterAI-1,iRow,nsAtom,nDimBC,nInter,BMx,
     &                          dMass,Lbl,Shift,qInt,Grad,AtomLbl,
     &                          Cx,.True.,iter)
            Energy(iterAI)=OS_Energy
            If (Max_OS.gt.0) Then
               If (UpMeth(4:4).ne.' ') UpMeth(5:6)='**'
               UpMeth(4:4)='+'
            End If
            Exit
         End If
        End Do
      End If
      Call mma_deallocate(Step_k)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Save the optimized kriging coordinates as the coordinates
*     for the next macro iteration.
*
      qInt(:,iter+1)=qInt(:,iterAI)
      Cx(:,iter+1)=Cx(:,iterAI)
*
*     Update the shift vector
*
      Shift(:,iter)=qInt(:,iter+1)-qInt(:,iter)
*
*     Update the predicted energy change
*
      ed = Energy(iterAI)-Energy(iter)
*
      Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter,
     &            Grad(1,iter),Shift(1,iter),Lbl)
*
*     Stick in the correct value for GrdMax, which might contain a
*     contribution due to constraints.
*
      GrdMax=GrdMax_Save
      GrdLbl=GrdLbl_Save
#ifdef _DEBUG_
      Call RecPrt('qInt(3):',' ',qInt,nInter,iter+1)
      Call RecPrt('Shift:',' ',Shift,nInter,iter)
      Write(6,*) 'UpMeth=',UpMeth
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocating memory used by Kriging
*
      Call mma_deallocate(Hessian)
      Call Close_Kriging()
*                                                                      *
************************************************************************
*                                                                      *
*-----Write out the shift in internal coordinate basis.
*
      If (iPrint.ge.99) Then
         Call RecPrt(
     &      'Shifts in internal coordinate basis / au or rad',
     &      ' ',Shift,nInter,Iter)
         Call RecPrt(
     &      'qInt in internal coordinate basis / au or rad',
     &      ' ',qInt,nInter,Iter+1)
      End If
*
*---- Remove unneeded fields from the runfile
      Dummy(1)=-Zero
      Call Put_dArray('BMxOld',Dummy(1),0)
      Call Put_dArray('TROld',Dummy(1),0)
*
      Call QExit('Update')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(kIter)
      If (.False.) Call Unused_integer(NmIter)
*
      Return
      End Subroutine Update_Kriging
