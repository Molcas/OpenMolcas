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
      Use AI, only: miAI, meAI, nspAI, blavAI, set_l
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
      Real*8, Allocatable:: Hessian(:,:), U(:,:), HTri(:), Temp(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Different hardwired kriging options
*
*     Note that turning of the sorting will result in a poorer
*     kriging!
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
*     Call RecPrt('Hessian',' ',Hessian,nInter,nInter)
*
      Call mma_allocate(U,nInter,nInter,Label='U')
      U(:,:)=0.0D0
      For all (i=1:nInter) U(i,i)=1.0D0
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
*     U(:,:)=0.0D0
*     For all (i=1:nInter) U(i,i)=1.0D0
*     Call TriPrt('HTri',' ',HTri,nInter)
      Hessian(:,:) = 0.0D0
      For all (i=1:nInter) Hessian(i,i)=HTri(i*(i+1)/2)
      Call mma_deallocate(HTri)
*     Call RecPrt('Hessian(D)',' ',Hessian,nInter,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Pass the data points to the GEK routine. Remember to
*     change the sign of the gradients.
*                                                                      *
************************************************************************
*                                                                      *
*     Note that we could have some kind of sorting here if we
*     like!
*
      Call mma_Allocate(qInt_s,nInter,nRaw,Label="qInt_s")
      Call mma_Allocate(Grad_s,nInter,nRaw,Label="Grad_s")
      Call mma_Allocate(Energy_s,nRaw,Label="Energy_s")
#ifdef _UNSORTED_
*
*     Tranform to the basis which diagonalize the HMF Hessian.
*
      Call Trans(U,qInt(1,iFirst),qInt_s,nInter,nRaw)
      Call Trans(U,Grad(1,iFirst),Grad_s,nInter,nRaw)
      Call DScal_(nInter*nRaw,-One,Grad_s,1)
      Call DCopy_(nRaw,Energy(iFirst),1,Energy_s,1)
*
#else
*                                                                      *
************************************************************************
*                                                                      *
*     Sort the data so that it some in an order of the points
*     closest to the last point. Make sure that the reference
*     point is the last point such that is always will define
*     the reference point. This is improtant when the bias is
*     defined Backrelative to the last points energy.
*
*     This code will have to be cleanup up later.
*
      ipCx_Ref=ipCx + (iter-1)*(3*nsAtom)
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
*
*
*     Tranform to the basis which diagonalize the HMF Hessian.
*
      Call mma_allocate(Temp,nInter,nRaw,Label='Temp')
      Call Trans(U,qInt_s,Temp,nInter,nRaw)
      qInt_s(:,:) = Temp
      Call Trans(U,Grad_s,Temp,nInter,nRaw)
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
*     coordinates, shifts, and gradients to update_sl_.f
*
      Call mma_allocate(qInt_s,nInter,MaxItr,Label='qInt_s')
      Call mma_allocate(Grad_s,nInter,MaxItr,Label='Grad_s')
      Call mma_allocate(Shift_s,nInter,MaxItr,Label='Shift_s')
      qInt_s(:,:)=0.0D0
      Grad_s(:,:)=0.0D0
      Shift_s(:,:)=0.0D0
      Call Trans(U,qInt,qInt_s,nInter,iter)
      Call Trans(U,Grad,Grad_s,nInter,iter)
      Call Trans(U,Shift,Shift_s,nInter,iter)
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
         Call set_l_kriging([Value_l],1)
         Call Put_dScalar('Value_l',Value_l)
      Else
         Call mma_Allocate(Array_l,nInter,Label='Array_l')
         Call Set_l_Array(Array_l,nInter,blavAI,Hessian)
         Call Set_l_Kriging(Array_l,nInter)
         Call mma_deAllocate(Array_l)
      End If
      Call mma_deallocate(Hessian)
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
         Beta_Disp_=MaxVal(GNrm,iterAI)*Beta_Disp
         Call Update_sl_(iterAI,iInt,nFix,nInter,
     &                qInt_s,Shift_s,Grad_s,iOptC,Beta,
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
         Call Dispersion_Kriging(qInt_s(1,iterAI+1),Dummy,nInter)
         E_Disp = Dummy
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
*        Check on convergence criterions.
*
         If (ThrEne.gt.Zero) Then
*           Convergence on energy criterions.
            Not_Converged = Abs(dEner).ge.ThrEne
            Not_Converged = Not_Converged .and. iterK.lt.miAI
            Not_Converged = Not_Converged .and. dqdq.lt.qBeta**2
         Else
*           Use standard convergence criterions
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
*     Save the optimized kriging coordinates as the coordinates
*     for the next macro iteration.
*
      Call mma_allocate(Temp,nInter,1,Label='Temp')
      Call BackTrans(U,qInt_s(1,iterAI),Temp,nInter,1)
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
*     Deallocating memory used by Kriging
*
      Call mma_deallocate(qInt_s)
      Call mma_deallocate(Grad_s)
      Call mma_deallocate(Shift_s)
      Call mma_deallocate(U)
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
*                                                                      *
************************************************************************
*                                                                      *
      Contains
      Subroutine Trans(U,X,Y,nInter,nIter)
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
      End Subroutine Trans
*
      Subroutine BackTrans(U,X,Y,nInter,nIter)
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
      Return
      End Subroutine BackTrans
*
      End Subroutine Update_Kriging
