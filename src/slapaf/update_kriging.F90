!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019,2020, Roland Lindh                                *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************
      Subroutine Update_kriging(Step_Trunc,nWndw)
!***********************************************************************
!                                                                      *
!     Object: to update coordinates                                    *
!                                                                      *
!    (see update_sl)                                                   *
!***********************************************************************
      Use kriging_mod, only: Max_Microiterations,                       &
     &                       Thr_microiterations, nSet
      Use Slapaf_Info, only: Cx, Gx, Shift, GNrm, Energy, qInt, dqInt,  &
     &                       Lbl
      use Slapaf_Parameters, only: UpMeth, Beta_Seed => Beta,           &
     &                             Beta_Disp_Seed => Beta_Disp, GrdLbl, &
     &                             GrdMax, E_Delta, ThrEne, ThrGrd,     &
     &                             ThrCons, nLambda, iter, NADC
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 dEner, E_Disp
      Integer HessIter
      Logical First_MicroIteration, Error, Found, CheckCons
      Character Step_Trunc
      Character GrdLbl_Save*8
      Real*8 Dummy(1), MaxErr, MaxErr_Ini
      Real*8, Allocatable:: ETemp(:,:), Hessian(:,:,:), Temp(:,:,:)
      Real*8, Parameter :: Beta_Disp_Min=1.0D-10
!                                                                      *
!***********************************************************************
!                                                                      *
!     Different hardwired kriging options
!
!#define _DEBUGPRINT_
!#define _OVERSHOOT_
!                                                                      *
!***********************************************************************
!                                                                      *
      Logical Kriging_Hessian, Not_Converged, Force_RS
#ifdef _OVERSHOOT_
      Real*8 :: OS_Disp(1), OS_Energy(1)
      Real*8, Allocatable:: Step_k(:,:)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface

         Subroutine SetUp_Kriging(nRaw,nInter,qInt,Grad,Energy,         &
     &                            Hessian_HMF,HDiag)
         Import :: nSet
         Integer, Intent(in) :: nRaw, nInter
         Real*8, Intent(in) :: qInt(nInter,nRaw),                       &
     &                         Grad(nInter,nRaw,nSet), Energy(nRaw,nSet)
         Real*8, Intent(inout), Optional:: Hessian_HMF(nInter,nInter),  &
     &                                     HDiag(nInter)
         End Subroutine SetUp_Kriging

      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *
      nQQ = SIZE(qInt,1)
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Update_Kriging: qInt',' ',qInt,nQQ,Iter)
      Call RecPrt('Update_Kriging: Shift',' ',Shift,nQQ,Iter-1)
      Call RecPrt('Update_Kriging: GNrm',' ',GNrm,Iter,1)
#endif
!
      Kriging_Hessian =.TRUE.
      Force_RS=.FALSE.
      iOpt_RS=1   ! Activate restricted variance.
      iterAI=iter
      dEner=Thr_microiterations
      iterK=0
      dqdq=Zero
      qBeta=Beta_Seed
      qBeta_Disp=Beta_Disp_Seed
      FAbs_Ini=Zero
      GrdMax_Save=GrdMax
      GrdLbl_Save=GrdLbl
      Call Qpg_iScalar('HessIter',Found)
      If (Found) Then
         Call Get_iScalar('HessIter',HessIter)
      Else
         HessIter = 0
      End If
      MaxErr_Ini=-One
      CheckCons=.False.
!                                                                      *
!***********************************************************************
!                                                                      *
!     Pick up the HMF Hessian
!
      Call mma_Allocate(Hessian,nQQ,nQQ,nSet,Label='Hessian')
      Call Mk_Hss_Q()
      Call Get_dArray('Hss_Q',Hessian(:,:,1),nQQ**2)
!     Call RecPrt('HMF Hessian',' ',Hessian(:,:,1),nQQ,nQQ)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Pass the sample points to the GEK procedure
!
      nRaw=Min(iter,nWndw/2)
      iFirst = iter - nRaw + 1
#ifdef _DEBUGPRINT_
      Call RecPrt('qInt(0)',  ' ',qInt(:,iFirst),nQQ,nRaw)
      Call RecPrt('Energy(0)',' ',Energy(iFirst),1,nRaw)
      Call RecPrt('dqInt(0)',  ' ',dqInt(1,iFirst),nQQ,nRaw)
      Call RecPrt('Shift',  ' ',Shift(1,iFirst),nQQ,nRaw)
#endif

      Call mma_allocate(ETemp,nRaw,nSet,Label='ETemp')
      Call mma_allocate(Temp,nQQ,nRaw,nSet,Label='Temp')
      Call Prepare_Kriging(ETemp,Temp,nRaw,nQQ,iFirst)
      Call Setup_Kriging(nRaw,nQQ,qInt(:,iFirst),Temp,ETemp,            &
     &                   Hessian_HMF=Hessian(:,:,1))
      Call mma_deallocate(ETemp)
      Call mma_deallocate(Temp)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Find the fraction value to be used to define the restricted
!     variance threshold.
!
!     Let the accepted variance be set as a fraction of the
!     largest component in the gradient.

      tmp=Zero
      Do i = 1, SIZE(Gx,2)
         Do j = 1, 3
            tmp = Max(tmp,Abs(Gx(j,i,iter)))
         End Do
      End Do
!
!     Temporary code until we have figured out this for constrained
!     optimizations.
!     (Note that for constrained optimizations, Beta is changed again
!      in con_opt, once the constrained gradient is available)
!
      Beta_Disp=Max(Beta_Disp_Min,tmp*Beta_Disp_Seed)
      Beta=Min(1.0D3*GNrm(iter),Beta_Seed)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     Start the Kriging loop.   The micro iterations
!
      Not_Converged = .True.
      Step_Trunc='N'  ! not defined
      Do While (Not_Converged)
!                                                                      *
!***********************************************************************
!                                                                      *
         kIter=iterAI
#ifdef _DEBUGPRINT_
         Write (6,*)
         Write (6,*) 'Do iterAI: ',iterAI
         Write (6,*) 'Beta=',Beta
#endif
!
!        Compute the Kriging Hessian
!
         First_MicroIteration=iterAI.eq.iter
         If (Kriging_Hessian) Then
            Call Hessian_Kriging_Layer(qInt(:,iterAI),Hessian,nQQ)
            If ((nSet > 3).And.NADC)                                    &
     &         Hessian(:,:,1) = Half*(Hessian(:,:,1)+Hessian(:,:,2))
            Call Put_dArray('Hss_Q',Hessian(:,:,1),nQQ**2)
!           Make fixhess treat the Hessian as if it was analytic.
            Call Put_iScalar('HessIter',iterAI)
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!        Compute the updated structure.
!
         nWndw_=nWndw/2 + (iterAI-iter)

         Call Update_inner(iterAI,Beta,Beta_Disp,Step_Trunc,nWndw_,     &
     &                     kIter,Kriging_Hessian,qBeta,iOpt_RS,         &
     &                     First_MicroIteration,iter,qBeta_Disp,.False.)

#ifdef _DEBUGPRINT_
         Write (6,*) 'After Update_inner: Step_Trunc=',Step_Trunc
         Call RecPrt('New Coord',' ',qInt,nQQ,iterAI+1)
         Call RecPrt('dqInt',' ',dqInt,nQQ,iterAI)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Save initial gradient
!     This is for the case the gradient is already converged, we want
!     the micro iterations to still reduce the gradient
!
!     Save the initial error in the constraints,
!     this is the value that should be used for global convergence
!
!     Try to enforce convergence in the constraints if the initial
!     gradient is close to convergence
!
         If (First_Microiteration) Then
            FAbs_Ini=GNrm(iterAI)/SQRT(DBLE(nQQ))
            If (nLambda.gt.0) Then
               If (FAbs_ini.le.Ten*ThrGrd) CheckCons=.True.
               Call Qpg_dScalar('Max error',Found)
               If (Found) Call Get_dScalar('Max error',MaxErr_Ini)
            End If
         End If
!
!        Attempt to break oscillations
!
         If (IterK.gt.2) Then
            dqdq1 = ddot_(nQQ,Shift(:,iterAI),1,Shift(:,iterAI-1),1)
!           If the overlap between three consecutive steps is negative,
!           cut down the step in half.
            If (dqdq1.lt.Zero) Then
               dqdq2 = ddot_(nQQ,Shift(:,iterAI-1),1,                   &
     &                       Shift(:,iterAI-2),1)
               dqdq3 = ddot_(nQQ,Shift(:,iterAI-2),1,                   &
     &                       Shift(:,iterAI-3),1)
               If ((dqdq2.lt.Zero).And.(dqdq3.lt.Zero)) Then
                  Shift(:,iterAI) = Half*Shift(:,iterAI)
                  qInt(:,iterAI+1) = qInt(:,iterAI)+Shift(:,iterAI)
#ifdef _DEBUGPRINT_
                  Write(6,*) 'Oscillation detected. '//                 &
     &                       'Step cut down in half!'
                  Call RecPrt('New Coord',' ',qInt,nQQ,iterAI+1)
#endif
               End If
            End If
         End If
!
!        Transform the new internal coordinates to Cartesians
!        (this updates the new internal coordinates, in case they were
!        not totally consistent)
!
         Error=(iterK.ge.1)
         iRef_Save=iRef
         iRef=iter ! Set the reference geometry
         Call NewCar_Kriging(iterAI,.True.,Error)
         iRef = iRef_Save
#ifdef _DEBUGPRINT_
         Call RecPrt('New Coord (after NewCar)','',qInt,nQQ,iterAI+1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!        In case of a constrained optimization GrdMax and GNrm have been
!        updated to reflect the contribution from the constraints.
!        GrdMax, however, is a scalar and we need to save it such that
!        it has the correct value on exit. That is, the value of itera-
!        tion iter.
!
         If (iterAI.eq.iter) Then
            GrdMax_Save=GrdMax
            GrdLbl_Save=GrdLbl
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!        During the micro iterations we need to set GNrm to some
!        reasonable value. Normally it is set by the force routine or
!        the con_opt routine.
!
         If (nLambda.eq.0 .and. iterAI.gt.iter) Then
            GNrm(iterAI)=Sqrt(DDot_(nQQ,dqInt(1,iterAI),1,              &
     &                                     dqInt(1,iterAI),1))
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!        Compute the step length from the last ab initio point
!        to the most recent kriging point.
!
         dqdq=Zero
         Do iInter=1,nQQ
            dqdq=dqdq+(qInt(iInter,iterAI+1)-qInt(iInter,iter))**2
         End Do
         If (iterK.eq.0.and.Step_trunc.eq.'*') dqdq=qBeta**2
!
!        In case of error during conversion to Cartesians
!
         If (Error) Then
            Step_Trunc='@'
            Exit
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!        Compute the energy and gradient according to the
!        surrogate model for the new coordinates.
!
         Call Kriging_Update(nQQ,iterAI+1,qInt(:,iterAI+1),E_Disp)
!                                                                      *
!***********************************************************************
!                                                                      *
         iterK  = iterK  + 1
         iterAI = iterAI + 1
         dEner = Energy(iterAI) - Energy(iterAI-1)
#ifdef _DEBUGPRINT_
         Call RecPrt('qInt(x):',' ',qInt,nQQ,iterAI)
         Call RecPrt('Ener(x):',' ',Energy,1,iterAI)
         Call RecPrt('dqInt(x):',' ',dqInt,nQQ,iterAI)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!        Check on convergence criteria.
!
         If (ThrEne.gt.Zero) Then
!           Convergence on energy criteria.
            Not_Converged = Abs(dEner).ge.ThrEne
            Not_Converged = Not_Converged .and. dqdq.lt.qBeta**2
         Else
!           Use standard convergence criteria
            FAbs=GNrm(iterAI-1)/SQRT(DBLE(nQQ))
            RMS =Sqrt(DDot_(nQQ,Shift(1,iterAI-1),1,                    &
     &                             Shift(1,iterAI-1),1)/DBLE(nQQ))
            RMSMx=Zero
            Do iInter = 1, nQQ
               RMSMx=Max(RMSMx,Abs(Shift(iInter,iterAI-1)))
            End Do
            GrdMx=GrdMax
#ifdef _DEBUGPRINT_
            Write (6,*)
            Write (6,*) 'iter=',iterAI-1
            Write (6,*) 'FAbs=',FAbs
            Write (6,*) 'GrdMx=',GrdMx
            Write (6,*) 'RMS=',RMS
            Write (6,*) 'RMSMx=',RMSMx
            Write (6,*) FAbs.gt.Min(ThrGrd,FAbs_ini)
            Write (6,*) GrdMx.gt.ThrGrd*OneHalf
            Write (6,*) RMS.gt.ThrGrd*Four
            Write (6,*) RMSMx.gt.ThrGrd*Six
            Write (6,*) 'Step_Trunc=',Step_Trunc
#endif
!           Ensure that the initial gradient is reduced,
!           except in the last micro iteration
            If (iterK.lt.Max_MicroIterations) Then
               Not_Converged = FAbs.gt.Min(ThrGrd,FAbs_ini)
            Else
               Not_Converged = FAbs.gt.ThrGrd
            End If
            Not_Converged = Not_Converged .or.                          &
     &                      GrdMx.gt.ThrGrd*OneHalf
            Not_Converged = Not_Converged .or.                          &
     &                      RMS.gt.ThrGrd*Four
            Not_Converged = Not_Converged .or.                          &
     &                      RMSMx.gt.ThrGrd*Six
            Not_Converged = Not_Converged .or.                          &
     &                      Step_Trunc.ne.' '
            If (CheckCons.and.(.not.Not_Converged).and.(nLambda.gt.0)   &
     &          .and.(iterK.lt.Max_Microiterations)) Then
               Call Qpg_dScalar('Max error',Found)
               If (Found) Then
                  Call Get_dScalar('Max error',MaxErr)
                  If (MaxErr.gt.ThrCons) Not_Converged=.True.
               End If
            End If
         End If
         If (Step_Trunc.eq.'.') Step_Trunc=' '
!        Check total displacement from initial structure
         If (Force_RS) Then
            Call mma_allocate(Temp,3,SIZE(Cx,2),1,Label='Temp')
            Temp(:,:,1)=Cx(:,:,iterAI)-Cx(:,:,iter)
            RMS = Sqrt(DDot_(3*SIZE(Cx,2),Temp,1,Temp,1)                &
     &          / DBLE(3*SIZE(Cx,2)) )
            If (RMS.gt.(Three*Beta)) Step_trunc='*'
            Call mma_deAllocate(Temp)
         End If
         If (Not_Converged.and.(Step_trunc.eq.' ') .and.                &
     &                   (iterK.ge.Max_Microiterations))  Step_trunc='#'
#ifdef _DEBUGPRINT_
         Write (6,*) 'Not_Converged=',Not_Converged
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!        If the step restriction is invoked or there is some other
!        reason signaled by Step_Trunc, terminate anyhow.
!
         If (Step_trunc.ne.' ') Not_Converged=.False.
!
!        If RS rather than RV do not micro iterate
!
         If (iOpt_RS.eq.0) Not_Converged=.False.
!        Not_Converged=.False. ! Force single iteration scheme.
#ifdef _DEBUGPRINT_
         Write (6,*) 'Not_Converged=',Not_Converged
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     End of the micro iteration loop
!
      End Do  ! Do While (Not_Converged)
!
!     Change label of updating method
!
      UpMeth='RVO   '
      Write(UpMeth(4:6),'(I3)') iterK
!
#ifdef _OVERSHOOT_
!                                                                      *
!***********************************************************************
!                                                                      *
!     Attempt overshooting
!
      Call mma_allocate(Step_k,nQQ,2)
      Call dCopy_(nQQ,qInt(1,iterAI),1,Step_k(1,2),1)
      Call dAXpY_(nQQ,-One,qInt(1,iter),1,Step_k(1,2),1)
!     Call RecPrt('q(iter+1)-q(f)','',Step_k(1,2),nQQ,1)
      If (iter.gt.1) Then
         Call dCopy_(nQQ,qInt(1,iter),1,Step_k(1,1),1)
         Call dAXpY_(nQQ,-One,qInt(1,iter-1),1,Step_k(1,1),1)
!        Call RecPrt('q(f)-q(f-1)','',Step_k(1,1),nQQ,1)
         dsds=dDot_(nQQ,Step_k(1,1),1,Step_k(1,2),1)
         dsds=dsds/Sqrt(ddot_(nQQ,Step_k(1,1),1,Step_k(1,1),1))
         dsds=dsds/Sqrt(ddot_(nQQ,Step_k(1,2),1,Step_k(1,2),1))
         Write(6,*) 'dsds = ',dsds
      Else
         dsds=Zero
      End If
      dsds_min=0.9D0
      If ((dsds.gt.dsds_min).And.(Step_Trunc.eq.' ')) Then
        Do Max_OS=9,0,-1
         OS_Factor=Max_OS*((dsds-dsds_min)/(One-dsds_min))**4
         Call dCopy_(nQQ,qInt(1,iter),1,qInt(1,iterAI),1)
         Call dAXpY_(nQQ,One+OS_Factor,Step_k(1,2),1,                   &
     &                                    qInt(1,iterAI),1)
         Call Energy_Kriging_Layer(qInt(1,iterAI),OS_Energy,nQQ)
         Call Dispersion_Kriging_Layer(qInt(1,iterAI),OS_Disp,nQQ)
         Write(6,*) 'Max_OS=',Max_OS
         Write(6,*) OS_Disp(1),E_Disp,Beta_Disp
         If ((OS_Disp(1).gt.E_Disp).And.(OS_Disp(1).lt.Beta_Disp))      &
     &      Then
            Call dAXpY_(nQQ,OS_Factor,Step_k(1,2),1,                    &
     &                                   Shift(1,iterAI-1),1)
            iRef_Save=iRef
            iRef=iter ! Set the reference geometry
            Call NewCar_Kriging(iterAI-1,.True.,Error)
            iRef = iRef_Save
            Energy(iterAI)=OS_Energy(1)
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
!                                                                      *
!***********************************************************************
!                                                                      *
!     Save the optimized kriging coordinates as the coordinates
!     for the next macro iteration.
!
      qInt(:,iter+1)=qInt(:,iterAI)
      Cx(:,:,iter+1)=Cx(:,:,iterAI)
!
!     Update the shift vector
!
      Shift(:,iter)=qInt(:,iter+1)-qInt(:,iter)
!
!     Update the predicted energy change
!
      E_Delta = Energy(iterAI)-Energy(iter)
!
      Call MxLbls(nQQ,dqInt(1,iter),Shift(1,iter),Lbl)
!
!     Stick in the correct value for GrdMax, which might contain a
!     contribution due to constraints.
!
      GrdMax=GrdMax_Save
      GrdLbl=GrdLbl_Save
#ifdef _DEBUGPRINT_
      Call RecPrt('qInt(3):',' ',qInt,nQQ,iter+1)
      Call RecPrt('Shift:',' ',Shift,nQQ,iter)
      Write(6,*) 'UpMeth=',UpMeth
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Deallocating memory used by Kriging
!
      Call mma_deallocate(Hessian)
      Call Finish_Kriging()
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Write out the shift in internal coordinate basis.
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Shifts in internal coordinate basis / au or rad',    &
     &            ' ',Shift,nQQ,Iter)
      Call RecPrt('qInt in internal coordinate basis / au or rad',      &
     &            ' ',qInt,nQQ,Iter+1)
#endif
!
!---- Remove unneeded fields from the runfile
      Dummy(1)=-Zero
      Call Put_dArray('BMxOld',Dummy(1),0)
      Call Put_dArray('TROld',Dummy(1),0)
!---- Restore previous values
      Call Put_iScalar('HessIter',HessIter)
      If (MaxErr_Ini.ge.Zero) Call Put_dScalar('Max error',MaxErr_Ini)
!
      Return
      End Subroutine Update_Kriging
