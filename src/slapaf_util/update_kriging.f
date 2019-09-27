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
     &                     rLambda,ipCx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,ipMF,
     &                     iOptH,HUpMet,kIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Degen,ThrEne,ThrGrd)
************************************************************************
*                                                                      *
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      iter           : iteration counter                              *
*      MaxItr         : max number of iteration                        *
*      NmIter         : number of iteration in numerical approach      *
*      iInt           : number of internal coordinates to vary         *
*      nFix           : number of frozen internal coordinates          *
*      nInter         : total number of internal coordinates           *
*      qInt(*,iter )  : the internal coordinates                       *
*      Grad(*,iter )  : the gradient in the internal coordinates       *
*      iOptC          : option flag for update methods                 *
*      Beta           : damping factor                                 *
*      Lbl            : character labels for internal coordinates      *
*      nLbl           : length of Lbl                                  *
*      GNrm           : the norm of the gradient in each iteration     *
*      Energy         : the energy of each iteration                   *
*      Line_Search    : logical flag for line search                   *
*      nLambda        : number of contraints                           *
*      iRow_c         : number of lines on the UDC file                *
*      nsAtom         : number of symmetry unique atoms                *
*      AtomLbl        : character string with atom labels              *
*      nSym           : number of irreps                               *
*      iOper          : integer representations of symmetry operators  *
*      mxdc           : max number of nsAtom                           *
*      jStab          : integer list of stabilizers                    *
*      nStab          : number of stabilizers                          *
*      BMx            : the so-called Wilson B matrix                  *
*      Smmtrc         : logical flag for symmetry properties           *
*      nDimBC         : dimension of redundant coordinates(?)          *
*      rLambda        : vector for Lagrange multipliers                *
*      ipCx           : pointer to cartesian coordinates               *
*      iNeg           : Hessian index                                  *
*      Labels         : character string of primitive int. coord.      *
*      nLabels        : length of Labels                               *
*      CnstWght       : constraints weight                             *
*                                                                      *
*    OutPut:                                                           *
*      Shift(*,iter ) : the shift of the internal coordinates          *
*      qInt(*,iter +1): the internal coordinates to be used in the     *
*                       next iteration                                 *
*      UpMeth         : character label with update method abrivation  *
*      ed             : estimated energy change to the next point      *
*      Step_Trunc     : character label to denote truncation type      *
*      GrdMax         : largest gradient component                     *
*      StpMax         : largest step component                         *
*      GrdLbl         : character label of component with GrdMax       *
*      StpLbl         : character label of component with StpLbl       *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             2000                                                     *
************************************************************************
      Use AI, only: miAI, meAI, nspAI
      Implicit Real*8 (a-h,o-z)
      External Restriction_Step, Restriction_Dispersion
      Real*8 Restriction_Step, Restriction_Dispersion
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 qInt(nInter,MaxItr), Shift(nInter,MaxItr),
     &       Grad(nInter,MaxItr), GNrm(MaxItr), Energy(MaxItr),
     &       BMx(3*nsAtom,3*nsAtom), rLambda(nLambda,MaxItr),
     &       dMass(nsAtom), Degen(3*nsAtom), dEner
      Integer iOper(0:nSym-1), jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),
     &        FindTS, TSC, HrmFrq_Show, Curvilinear
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6
*
      Logical Kriging_Hessian, Not_Converged, Single_l_value
      Real*8, Allocatable:: Array_l(:)
*define _UNSORTED_
#ifndef _UNSORTED_
      Real*8, Allocatable:: Energy_s(:), qInt_s(:,:), Grad_s(:,:)
#endif
*
      iRout=153
      iPrint=nPrint(iRout)
      Call QEnter('Update')
*
      If (iPrint.ge.99) Then
         Call RecPrt('Update_K: qInt',' ',qInt,nInter,Iter)
         Call RecPrt('Update_K: Shift',' ',Shift,nInter,Iter-1)
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
*     Pass the data point to the GEK routine. Remember to
*     change the sign of the gradients.
*                                                                      *
************************************************************************
*                                                                      *
*     Note that we could have some kind of sorting here if we
*     like!
*
#ifdef _UNSORTED_
      Call DScal_(nInter*nRaw,-One,Grad(1,iFirst),1)
      Call Start_Kriging(nRaw,nInter,
     &                      qInt(1,iFirst),
     &                      Grad(1,iFirst),
     &                      Energy(iFirst))
      Call DScal_(nInter*nRaw,-One,Grad(1,iFirst),1)
#else
*                                                                      *
************************************************************************
*                                                                      *
*     Sort the data so that it some in an order of the points
*     closest to the last point. Make sure that the reference
*     point is the last point such that is always will define
*     the reference point. This is improtant when the bias is
*     defined relative to the last points energy.
*
*     This code will have to be cleanup up later.
*
      ipCx_Ref=ipCx + (iter-1)*(3*nsAtom)
      Call mma_Allocate(Energy_s,nRaw,Label="Energy_s")
      Call mma_Allocate(qInt_s,nInter,nRaw,Label="qInt_s")
      Call mma_Allocate(Grad_s,nInter,nRaw,Label="Grad_s")
*
      Call DCopy_(nInter,qInt(1,iter),1,qInt_s(1,nRaw),1)
      Call DCopy_(nInter,Grad(1,iter),1,Grad_s(1,nRaw),1)
      Energy_s(nRaw)=Energy(iter)
*
*     Pick up the coordinates in descending order with the ones
*     that are the closest to the current structure.
*
      iSt=Max(1,iter-nWndw+1)
      Thr_low = 0.0D0
      Thr_high= 99.0D0
      Do iRaw = nRaw-1, 1, -1
*
         kter=-1
         Do jter = iSt, iter-1
*
*           Compute the distance in Cartesian coordinates.
*
            Distance=Zero
            Do ix = 1, 3*nsAtom
               Distance= Distance +
     &                   Degen(ix) *
     &                  (Work(ipCx_ref+ix-1) -
     &                   Work(ipCx + (jter-1)*3*nsAtom+ix-1))**2
            End Do
            Distance = sqrt(Distance)
*
            If (Distance.gt.Thr_low .and.
     &          Distance.lt.Thr_high) Then
               kter=jter
               Thr_high=Distance
            End If
         End Do
         If (kter.eq.-1) Then
            Write (6,*) 'kter not set!'
         Else
            Call DCopy_(nInter,qInt(1,kter),1,qInt_s(1,iRaw),1)
            Call DCopy_(nInter,Grad(1,kter),1,Grad_s(1,iRaw),1)
            Energy_s(iRaw)=Energy(kter)
            Thr_low=Thr_high
            Thr_high= 99.0D0
         End If
      End Do
#ifdef _DEBUG_
      Call RecPrt('qInt_s(s)',  ' ',qInt_s,nInter,nRaw)
      Call RecPrt('Energy_s(s)',' ',Energy_s,1,nRaw)
      Call RecPrt('Grad_s(s)',  ' ',Grad_s,nInter,nRaw)
#endif
*
      Call DScal_(nInter*nRaw,-One,Grad_s,1)
      Call Start_Kriging(nRaw,nInter,
     &                      qInt_s,
     &                      Grad_s,
     &                      Energy_s)
*
      Call mma_deAllocate(Energy_s)
      Call mma_deAllocate(qInt_s)
      Call mma_deAllocate(Grad_s)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Update the l value dynamically. Here we compare the actual,
*     ab inito value with the GEK prediction of the gradient.
*
      Call Get_dScalar('Value_l',Value_l)
*
      If (iter.gt.nspAI) Then
         iOld=iFirst+nRaw-2
         xxx=DDot_(nInter,Grad(1,iOld),1,Grad(1,iOld),1)
         iNew=iFirst+nRaw-1
         yyy=DDot_(nInter,Grad(1,iNew),1,Grad(1,iNew),1)
#define _UPDATE_L_
#ifdef _UPDATE_L_
         If (yyy.gt.xxx) Then
            Value_l=Value_l * 0.95D0
         Else
            Value_l=Value_l * 1.05D0
         End If
#endif
*
*        Update the restricted variance threshold.
*
*        Beta_Disp=Max(Abs(Energy(iNew)-Energy(iOld)),
*    &                 1.0D-6)
      End If
*
*     Select between single or multiple l values.
*     Multiple option still experimental.
*
      Single_l_value=.True.
*     Single_l_value=.False.
      If (Single_l_value) Then
         Call set_l_kriging([Value_l],1)
      Else
         Call mma_Allocate(Array_l,nInter,Label='Array_l')
*        Call DCopy_(nInter,[1.0D0],0,Array_l,1)
         Call Set_l_Array(Array_l,nInter)
*        Call RecPrt('l values',' ',Array_l,nInter,1)
*        Write (6,*) 'Value_l=',Value_l
         Call DScal_(nInter,Value_l,Array_l,1)
*        Call RecPrt('l values',' ',Array_l,nInter,1)
         Call set_l_kriging(Array_l,nInter)
         Call mma_DeAllocate(Array_l)
      End If
      Call Put_dScalar('Value_l',Value_l)
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
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the updated structure.
*
         Call Update_sl_(iterAI,iInt,nFix,nInter,
     &                qInt,Shift,Grad,iOptC,Beta,Beta_Disp,
     &                Lbl,GNrm,Energy,
     &                UpMeth,ed,Line_Search,Step_Trunc,nLambda,
     &                iRow_c,nsAtom,AtomLbl,nSym,iOper,mxdc,jStab,
     &                nStab,BMx,Smmtrc,nDimBC,rLambda,ipCx,
     &                GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &                Labels,nLabels,FindTS,TSC,nRowH,
     &                nWndw/2,Mode,ipMF,
     &                iOptH,HUpMet,kIter_,GNrm_Threshold,IRC,dMass,
     &                HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                Kriging_Hessian,qBeta,Restriction_dispersion,
     &                iOpt_RS)
*
*        Change lable of updating method if kriging points has
*        been used.
*
         If (iterK.gt.0) Then
            UpMeth='GPR   '
            Write (UpMeth(4:6),'(I3)') iterK
         Else
            UpMeth='RV-RFO'
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the step length from the last ab inito point
*        to the most recent krining point.
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
*        Compute the energy and gradient according to the
*        surrogate model for the new coordinates.
*
         Call Energy_Kriging(qInt(1,iterAI+1),Energy(iterAI+1),nInter)
         Call Dispersion_Kriging(qInt(1,iterAI+1),Dummy,nInter)
         E_Disp = Dummy
         Call Gradient_Kriging(qInt(1,iterAI+1),Grad(1,iterAI+1),nInter)
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
*        Check on convergence criterions.
*
         If (ThrEne.gt.Zero) Then
*           Convergence on energy criterions.
            Not_Converged = Abs(dEner).ge.ThrEne
            Not_Converged = Not_Converged .and. iterK.lt.miAI
            Not_Converged = Not_Converged .and. dqdq.lt.qBeta**2
         Else
*           Use standard convergence criterions
            FAbs=Sqrt(DDot_(nInter,Grad(1,iterAI),1,
     &                            Grad(1,iterAI),1)/DBLE(nInter))
            RMS =Sqrt(DDot_(nInter,Shift(1,iterAI-1),1,
     &                          Shift(1,iterAI-1),1)/DBLE(nInter))
            GrdMx=Zero
            RMSMx=Zero
            Do iInter = 1, nInter
               GrdMx=Max(GrdMx,Abs(Grad(iInter,iterAI)))
               RMSMx=Max(RMSMx,Abs(Shift(iInter,iterAI-1)))
            End Do
*
            Not_Converged = FAbs.gt.ThrGrd
            Not_Converged = Not_Converged .or.
     &                      GrdMx.gt.ThrGrd*OneHalf
            Not_Converged = Not_Converged .or.
     &                      RMS.gt.ThrGrd*Four
            Not_Converged = Not_Converged .or.
     &                      RMSMx.gt.ThrGrd*Six
            Not_Converged = Not_Converged .and. iterK.le.miAI
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        If the step restriction is envoked terminate anyhow.
*
         If (Step_trunc.eq.'*') Not_Converged=.False.
*                                                                      *
************************************************************************
*                                                                      *
*     End of the micro iteration loop
*
      End Do  ! Do While
*                                                                      *
************************************************************************
*                                                                      *
*     Reduce the l-value(s) if the step was restricted.
*
#ifdef _UPDATE_L_
      If (Step_trunc.eq.'*') Then
         Call Get_dScalar('Value_l',Value_l)
         Value_l=Value_l * 0.95D0
*        Write (6,*) ' Set l value to:',Value_l
         Call Put_dScalar('Value_l',Value_l)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Save the optimized kriging coordinates as the coordinates
*     for the next macro iteration.
*
      Call DCopy_(nInter,qInt(1,iterAI),1,qInt(1,iter+1),1)
      Call DCopy_(nInter,qInt(1,iter+1),1,Shift(1,iter),1)
      Call DaXpY_(nInter,-One,qInt(1,iter),1,Shift(1,iter),1)
      Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter,
     &            Grad(1,iter),Shift(1,iter),Lbl)
#ifdef _DEBUG_
      Call RecPrt('qInt(3):',' ',qInt,nInter,iter+1)
      Call RecPrt('Shift(3):',' ',Shift,nInter,iter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocating memory used by Kriging
*
      Call Finish_Kriging()
*                                                                      *
************************************************************************
*                                                                      *
*-----Write out the shift in internal coordinates basis.
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
      End
