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
* Copyright (C) 2019, Roland Lindh                                     *
************************************************************************
      Subroutine Update_kriging(
     &                     iter,MaxItr,NmIter,iInt,nFix,nInter,qInt,
     &                     Shift,
     &                     Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,
     &                     Energy,UpMeth,ed,Line_Search,Step_Trunc,
     &                     nLambda,iRow_c,nsAtom,AtomLbl,nSym,iOper,
     &                     mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                     rLambda,ipCx,Gx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,ipMF,
     &                     iOptH,HUpMet,kIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Degen,ThrEne,ThrGrd,iRow)
************************************************************************
*                                                                      *
*     Object: to update coordinates                                    *
*                                                                      *
*    (see update_sl)                                                   *
************************************************************************
      Use kriging_mod, only: miAI, meAI, blavAI, set_l
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 qInt(nInter,MaxItr), Shift(nInter,MaxItr),
     &       Grad(nInter,MaxItr), GNrm(MaxItr), Energy(MaxItr),
     &       BMx(3*nsAtom,3*nsAtom), rLambda(nLambda,MaxItr),
     &       dMass(nsAtom), Degen(3*nsAtom), dEner,
     &       Gx(3*nsatom,Iter)
      Integer iOper(0:nSym-1), jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),
     &        FindTS, TSC, HrmFrq_Show, Curvilinear
*    &        FindTS, TSC, HrmFrq_Show, Curvilinear, Print_it
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6
      Real*8, Allocatable:: Hessian(:,:), U(:,:), HTri(:), Temp(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Different hardwired kriging options
*
*     Note that turning off the sorting will result in a poorer kriging!
*
*#define _UNSORTED_
#define _DIAG_HESS_
*                                                                      *
************************************************************************
*                                                                      *
      Logical Kriging_Hessian, Not_Converged
      Real*8, Allocatable:: Array_l(:)
      Real*8, Allocatable:: Energy_s(:)
      Real*8, Allocatable:: qInt_s(:,:), Grad_s(:,:), Shift_s(:,:)
*
      iRout=153
      iPrint=nPrint(iRout)
      Call QEnter('Update')
*
      If (iPrint.ge.99) Then
         Call RecPrt('Update_K: qInt',' ',qInt,nInter,Iter)
         Call RecPrt('Update_K: Shift',' ',Shift,nInter,Iter-1)
         Call RecPrt('Update_K: GNrm',' ',GNrm,Iter,1)
         Call RecPrt('Update_K: GNrm',' ',GNrm,Iter,1)
      End If
*
      Kriging_Hessian =.TRUE.
      iOpt_RS=1   ! Activate restricted variance.
      iterAI=iter
      dEner=meAI
      nRaw=Min(iter,nWndw/2)
      iFirst = iter - nRaw + 1
      iterK=0
      dqdq=Zero
      qBeta=Beta
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
*
      Call mma_allocate(U,nInter,nInter,Label='U')
      U(:,:)=Zero
      Forall (i=1:nInter) U(i,i)=One
#ifdef _DIAG_HESS_
*
      Call mma_allocate(HTri,nInter*(nInter+1)/2,Label='HTri')
      Do iInter = 1, nInter
         Do jInter = 1, iInter
            ij = iInter*(iInter-1)/2 + jInter
            HTri(ij)=Hessian(iInter,jInter)
         End Do
      End Do
*     Call TriPrt('HTri(raw)',' ',HTri,nInter)
      Call NIDiag_new(HTri,U,nInter,nInter,0)
*     U(:,:)=Zero
*     Forall (i=1:nInter) U(i,i)=One
*     Call TriPrt('HTri',' ',HTri,nInter)
      Hessian(:,:) = Zero
      Forall (i=1:nInter) Hessian(i,i)=HTri(i*(i+1)/2)
      Call mma_deallocate(HTri)
*     Call RecPrt('Hessian(D)',' ',Hessian,nInter,nInter)
#endif
#ifdef _DEBUG_
      Call RecPrt('U','',U,nInter,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Select between setting all l's to a single value or go in
*     multiple l-value mode in which the l-value is set such that
*     the kriging hessian reproduce the diagonal value of the HMF
*     Hessian of the current structure.
*
      If (Set_l) Then
         Call Get_dScalar('Value_l',Value_l)
      Else
         Call mma_Allocate(Array_l,nInter,Label='Array_l')
         Call Set_l_Array(Array_l,nInter,blavAI,Hessian)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Pass the data points to the GEK routine. Remember to
*     change the sign of the gradients.
*                                                                      *
************************************************************************
*                                                                      *
*     Note that we could have some kind of sorting here if we like!
*
      Call mma_Allocate(qInt_s,nInter,nRaw,Label="qInt_s")
      Call mma_Allocate(Grad_s,nInter,nRaw,Label="Grad_s")
      Call mma_Allocate(Energy_s,nRaw,Label="Energy_s")
#ifdef _UNSORTED_
*
*     Transform to the basis which diagonalizes the HMF Hessian.
*
      Call Trans_K(U,qInt(1,iFirst),qInt_s,nInter,nRaw)
      Call Trans_K(U,Grad(1,iFirst),Grad_s,nInter,nRaw)
      Call DScal_(nInter*nRaw,-One,Grad_s,1)
      Call DCopy_(nRaw,Energy(iFirst),1,Energy_s,1)
*
#else
*                                                                      *
************************************************************************
*                                                                      *
*     Sort the data so that the points are ordered by distance
*     to the last point. Make sure that the reference point is
*     the last point. This is important when the bias is defined
*     relative to the last point energy.
*
*     This code will have to be cleaned up up later.
*
      ipCx_Ref=ipCx + (iter-1)*(3*nsAtom)
*
      Call DCopy_(nInter,qInt(1,iter),1,qInt_s(1,nRaw),1)
      Call DCopy_(nInter,Grad(1,iter),1,Grad_s(1,nRaw),1)
      Energy_s(nRaw)=Energy(iter)
*
*     Pick up the coordinates in descending order starting with the ones
*     that are the closest to the current structure.
*
*     Print_it=.False.
*666  Continue
      iSt=Max(1,iter-nWndw+1)
      Thr_low = Zero
      Thr_high= 1.0D99
      Do iRaw = nRaw-1, 1, -1
*
         kter=-1
         Do jter = iSt, iter-1
*
*           Compute the distance in Cartesian coordinates.
*
#ifdef _CARTESIAN_
            Distance=Zero
            Do ix = 1, 3*nsAtom
               Distance= Distance +
     &                   Degen(ix) *
     &                  (Work(ipCx_ref+ix-1) -
     &                   Work(ipCx + (jter-1)*3*nsAtom+ix-1))**2
            End Do
#else
            Distance=0.0d0
            If (set_l) Then
               Do inter = 1, nInter
                  Distance = Distance +
     &                       (
     &                        ( qInt(inter,iter) - qInt(inter,jter) )
     &                       / Value_l
     &                       )**2
               End Do
            Else
               Do inter = 1, nInter
                  Distance = Distance +
     &                       (
     &                        ( qInt(inter,iter) - qInt(inter,jter) )
     &                       / Array_l(inter)
     &                       )**2
               End Do
            End If
#endif
            Distance = sqrt(Distance)
            If (Print_it) Then
               Write (6,*) 'iRaw,jter,kter=',iRaw,jter,kter
               Write (6,*) 'Distance=',Distance
               Write (6,*) 'Thr_low=',Thr_low
               Write (6,*) 'Thr_high=',Thr_high
            End If
*
            If (Distance.gt.Thr_low .and.
     &          Distance.lt.Thr_high) Then
               kter=jter
               Thr_high=Distance
            End If
*           If (Print_it) Then
*              Write (6,*) 'iRaw,jter,kter=',iRaw,jter,kter
*              Write (6,*) 'Thr_high=',Thr_high
*           End If
         End Do
         If (kter.eq.-1) Then
            Write (6,*) 'kter not set!'
*           If (Print_it) Then
               Call Abend()
*           Else
*              Print_it=.True.
*              Go To 666
*           End if
         Else
*           Write (6,*) 'Use iteration: kter=',kter
            Call DCopy_(nInter,qInt(1,kter),1,qInt_s(1,iRaw),1)
            Call DCopy_(nInter,Grad(1,kter),1,Grad_s(1,iRaw),1)
            Energy_s(iRaw)=Energy(kter)
            Thr_low=Thr_high
            Thr_high= 1.0D99
         End If
      End Do
#ifdef _DEBUG_
      Call RecPrt('qInt_s(s)',  ' ',qInt_s,nInter,nRaw)
      Call RecPrt('Energy_s(s)',' ',Energy_s,1,nRaw)
      Call RecPrt('Grad_s(s)',  ' ',Grad_s,nInter,nRaw)
#endif
*
*     Transform to the basis which diagonalizes the HMF Hessian.
*
      Call mma_allocate(Temp,nInter,nRaw,Label='Temp')
      Call Trans_K(U,qInt_s,Temp,nInter,nRaw)
      qInt_s(:,:) = Temp
      Call Trans_K(U,Grad_s,Temp,nInter,nRaw)
      Grad_s(:,:) = -Temp
      Call mma_deallocate(Temp)
*
#endif
      Call Start_Kriging(nRaw,nInter,qInt_s,Grad_s,Energy_s)
*
      Call mma_deAllocate(Energy_s)
      Call mma_deAllocate(qInt_s)
      Call mma_deAllocate(Grad_s)
*                                                                      *
************************************************************************
*                                                                      *
*     There is a small tweak here. We need to send the transformed
*     coordinates, shifts, and gradients to update_sl_
*
      Call mma_allocate(qInt_s,nInter,MaxItr,Label='qInt_s')
      Call mma_allocate(Grad_s,nInter,MaxItr,Label='Grad_s')
      Call mma_allocate(Shift_s,nInter,MaxItr,Label='Shift_s')
      qInt_s(:,:)=Zero
      Grad_s(:,:)=Zero
      Shift_s(:,:)=Zero
      Call Trans_K(U,qInt,qInt_s,nInter,iter)
      Call Trans_K(U,Grad,Grad_s,nInter,iter)
      Call Trans_K(U,Shift,Shift_s,nInter,iter)
*                                                                      *
************************************************************************
*                                                                      *
*     Select between setting all l's to a single value or go in
*     multiple l-value mode in which the l-value is set such that
*     the kriging hessian reproduce the diagonal value of the HMF
*     Hessian of the current structure.
*
      If (Set_l) Then
         Call set_l_kriging([Value_l],1)
         Call Put_dScalar('Value_l',Value_l)
      Else
         Call Set_l_Kriging(Array_l,nInter)
         Call mma_deAllocate(Array_l)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Let the accepted variance be set as a linear function of the
*     largest component in the gradient.
*     Make sure that the variance restriction never is too tight, hence
*     the limit to 0.001 au = 0.63 kcal/mol
*
      tmp=0.0D0
      Do i = 1, 3*nsAtom
         tmp = Max(tmp,Abs(Gx(i,iter)))
      End Do
      Beta_Disp_=Max(0.001D0,tmp*Beta_Disp)
*
      Beta_=Beta
*
#ifdef _RS_RFO_
*     Switch over to RS-RFO once the gradient is low.
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
      do while (Not_Converged)
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
            Call Hessian_Kriging(qInt_s(1,iterAI),Hessian,nInter)
            Call Put_dArray('Hss_Q',Hessian,nInter**2)
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the updated structure.
*
         Call Update_sl_(iterAI,iInt,nFix,nInter,
     &                qInt_s,Shift_s,Grad_s,iOptC,Beta_,
     &                Beta_Disp_,
     &                Lbl,GNrm,Energy,
     &                UpMeth,ed,Line_Search,Step_Trunc,nLambda,
     &                iRow_c,nsAtom,AtomLbl,nSym,iOper,mxdc,jStab,
     &                nStab,BMx,Smmtrc,nDimBC,rLambda,ipCx,
     &                GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &                Labels,nLabels,FindTS,TSC,nRowH,
     &                nWndw/2,Mode,ipMF,
     &                iOptH,HUpMet,kIter_,GNrm_Threshold,IRC,dMass,
     &                HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                Kriging_Hessian,qBeta,iOpt_RS)
*
*        Change label of updating method if kriging points have
*        been used.
*
         If (iterK.gt.0) Then
            UpMeth='GEK   '
            Write (UpMeth(4:6),'(I3)') iterK
         Else If(iOpt_RS.eq.1) Then
            UpMeth='RV-RFO'
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the step length from the last ab initio point
*        to the most recent kriging point.
*
         dqdq=Zero
         Do iInter=1,nInter
            dqdq=(qInt_s(iInter,iterAI+1)
     &           -qInt_s(iInter,iFirst+nRaw-1))**2
         End Do
         If (iterK.eq.0.and.Step_trunc.eq.'*') dqdq=qBeta**2
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the energy and gradient according to the
*        surrogate model for the new coordinates.
*
         Call Energy_Kriging(qInt_s(1,iterAI+1),Energy(iterAI+1),nInter)
         Call Dispersion_Kriging(qInt_s(1,iterAI+1),E_Disp,nInter)
         Call Gradient_Kriging(qInt_s(1,iterAI+1),
     &                         Grad_s(1,iterAI+1),nInter)
         Call DScal_(nInter,-One,Grad_s(1,iterAI+1),1)
*                                                                      *
************************************************************************
*                                                                      *
         iterK  = iterK  + 1
         iterAI = iterAI + 1
         dEner = Energy(iterAI) - Energy(iterAI-1)
#ifdef _DEBUG_
         Call RecPrt('qInt(x):',' ',qInt_s,nInter,iterAI)
         Call RecPrt('Ener(x):',' ',Energy,1,iterAI)
         Call RecPrt('Grad(x):',' ',Grad_s,nInter,iterAI)
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
            FAbs=Sqrt(DDot_(nInter,Grad_s(1,iterAI),1,
     &                            Grad_s(1,iterAI),1)/DBLE(nInter))
            RMS =Sqrt(DDot_(nInter,Shift_s(1,iterAI-1),1,
     &                          Shift_s(1,iterAI-1),1)/DBLE(nInter))
            GrdMx=Zero
            RMSMx=Zero
            Do iInter = 1, nInter
               GrdMx=Max(GrdMx,Abs(Grad_s(iInter,iterAI)))
               RMSMx=Max(RMSMx,Abs(Shift_s(iInter,iterAI-1)))
            End Do
*
            Not_Converged = FAbs.gt.ThrGrd
            Not_Converged = Not_Converged .or.
     &                      GrdMx.gt.ThrGrd*OneHalf
            Not_Converged = Not_Converged .or.
     &                      RMS.gt.ThrGrd*Four
            Not_Converged = Not_Converged .or.
     &                      RMSMx.gt.ThrGrd*Six
         End If
         Not_Converged = Not_Converged .and. iterK.lt.miAI
*                                                                      *
************************************************************************
*                                                                      *
*        If the step restriction is invoked, terminate anyhow.
*
         If (Step_trunc.eq.'*') Not_Converged=.False.
*
*        If RS rather than RV don not micro iterate
*
         If (iOpt_RS.eq.0) Not_Converged=.False.
*                                                                      *
************************************************************************
*                                                                      *
*     End of the micro iteration loop
*
      End Do  ! Do While
*                                                                      *
************************************************************************
*                                                                      *
*     Save the optimized kriging coordinates as the coordinates
*     for the next macro iteration.
*
      Call mma_allocate(Temp,nInter,1,Label='Temp')
      Call BackTrans_K(U,qInt_s(1,iterAI),Temp,nInter,1)
      qInt(:,iter+1)=Temp(:,1)
      Call mma_deallocate(Temp)
*
*     Update the shift vector
*
      Call DCopy_(nInter,qInt(1,iter+1),1,Shift(1,iter),1)
      Call DaXpY_(nInter,-One,qInt(1,iter),1,Shift(1,iter),1)
*
      Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter,
     &            Grad(1,iter),Shift(1,iter),Lbl)
#ifdef _DEBUG_
      Call RecPrt('qInt(3):',' ',qInt,nInter,iter+1)
      Call RecPrt('Shift:',' ',Shift,nInter,iter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Reset Step_trunc flag if gradient is below the convergence
*     theshold. If not we might end up in a loop of no displacements
*     and eventually numerical problems.
*
      FAbs=Sqrt(
     &          DDot_(nInter,Grad(1,iter),1,Grad(1,iter),1)
     &         / DBLE(nInter)
     &         )
      GrdMax=Abs(GrdMax)
      If (GrdMax.le.ThrGrd*OneHalf .and. Fabs.lt.ThrGrd) Step_trunc=' '
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocating memory used by Kriging
*
      Call mma_deallocate(qInt_s)
      Call mma_deallocate(Grad_s)
      Call mma_deallocate(Shift_s)
      Call mma_deallocate(U)
      Call mma_deallocate(Hessian)
      Call Finish_Kriging()
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
      Call Put_dArray('BMxOld',Work(ip_Dummy),0)
      Call Put_dArray('TROld',Work(ip_Dummy),0)
*
      Call QExit('Update')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(kIter)
      If (.False.) Call Unused_integer(NmIter)
*
      Return
      End Subroutine Update_Kriging
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Trans_K(U,X,Y,nInter,nIter)
      Implicit None
      Integer nInter, nIter
      Real*8 U(nInter,nInter), X(nInter,nIter), Y(nInter,nIter)
*
*     Call RecPrt('U',' ',U,nInter,nInter)
*     Call RecPrt('X',' ',X,nInter,nIter)
      Call DGEMM_('T','N',nInter,nIter,nInter,
     &            1.0D0,U,nInter,
     &                  X,nInter,
     &            0.0D0,Y,nInter)
*     Call RecPrt('Y',' ',Y,nInter,nIter)
*
      Return
      End Subroutine Trans_K
*
      Subroutine BackTrans_K(U,X,Y,nInter,nIter)
      Implicit None
      Integer nInter, nIter
      Real*8 U(nInter,nInter), X(nInter,nIter), Y(nInter,nIter)
*
*     Call RecPrt('U',' ',U,nInter,nInter)
*     Call RecPrt('X',' ',X,nInter,nIter)
      Call DGEMM_('N','N',nInter,nIter,nInter,
     &            1.0D0,U,nInter,
     &                  X,nInter,
     &            0.0D0,Y,nInter)
*     Call RecPrt('Y',' ',Y,nInter,nIter)
*
      End Subroutine BackTrans_K
