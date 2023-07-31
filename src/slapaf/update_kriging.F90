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

subroutine Update_kriging(Step_Trunc,nWndw)
!***********************************************************************
!                                                                      *
!     Object: to update coordinates                                    *
!                                                                      *
!    (see update_sl)                                                   *
!***********************************************************************

use Slapaf_Info, only: Beta_Disp_Seed => Beta_Disp, Beta_Seed => Beta, Cx, dqInt, E_Delta, Energy, GNrm, GrdLbl, GrdMax, Gx, iter, &
                       Lbl, NADC, nLambda, qInt, Shift, ThrCons, ThrEne, ThrGrd, UpMeth
use Kriging_mod, only: Max_Microiterations, nSet, Thr_microiterations
use Kriging_procedures, only: Setup_Kriging
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Four, Six, Ten, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nWndw
character, intent(out) :: Step_Trunc
integer(kind=iwp) :: HessIter, i, iFirst, iInter, iOpt_RS, iRef, iRef_Save, iterAI, iterK, j, kIter, nQQ, nRaw, nWndw_
real(kind=wp) :: Beta, Beta_Disp, dEner, dqdq, dqdq1, dqdq2, dqdq3, Dummy(1), E_Disp, FAbs, FAbs_Ini, GrdMax_Save, GrdMx, MaxErr, &
                 MaxErr_Ini, qBeta, qBeta_Disp, RMS, RMSMx, tmp
character(len=8) :: GrdLbl_Save
logical(kind=iwp) :: CheckCons, Error, First_MicroIteration, Force_RS, Found, Kriging_Hessian, Not_Converged
real(kind=wp), allocatable :: ETemp(:,:), Hessian(:,:,:), Temp(:,:,:)
real(kind=wp), parameter :: Beta_Disp_Min = 1.0e-10_wp
real(kind=wp), external :: DDot_
!                                                                      *
!***********************************************************************
!                                                                      *
!     Different hardwired kriging options
!
!#define _DEBUGPRINT_
!#define _OVERSHOOT_
#ifdef _OVERSHOOT_
real(kind=wp) :: OS_Disp(1), OS_Energy(1)
real(kind=wp), allocatable :: Step_k(:,:)
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
nQQ = size(qInt,1)

#ifdef _DEBUGPRINT_
call RecPrt('Update_Kriging: qInt',' ',qInt,nQQ,Iter)
call RecPrt('Update_Kriging: Shift',' ',Shift,nQQ,Iter-1)
call RecPrt('Update_Kriging: GNrm',' ',GNrm,Iter,1)
#endif

Kriging_Hessian = .true.
Force_RS = .false.
iOpt_RS = 1   ! Activate restricted variance.
iterAI = iter
dEner = Thr_microiterations
iterK = 0
dqdq = Zero
qBeta = Beta_Seed
qBeta_Disp = Beta_Disp_Seed
FAbs_Ini = Zero
GrdMax_Save = GrdMax
GrdLbl_Save = GrdLbl
call Qpg_iScalar('HessIter',Found)
if (Found) then
  call Get_iScalar('HessIter',HessIter)
else
  HessIter = 0
end if
MaxErr_Ini = -One
CheckCons = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up the HMF Hessian

call mma_Allocate(Hessian,nQQ,nQQ,nSet,Label='Hessian')
call Mk_Hss_Q()
call Get_dArray('Hss_Q',Hessian(:,:,1),nQQ**2)
!call RecPrt('HMF Hessian',' ',Hessian(:,:,1),nQQ,nQQ)
!                                                                      *
!***********************************************************************
!                                                                      *
! Pass the sample points to the GEK procedure

nRaw = min(iter,nWndw/2)
iFirst = iter-nRaw+1
#ifdef _DEBUGPRINT_
call RecPrt('qInt(0)',' ',qInt(:,iFirst),nQQ,nRaw)
call RecPrt('Energy(0)',' ',Energy(iFirst),1,nRaw)
call RecPrt('dqInt(0)',' ',dqInt(:,iFirst),nQQ,nRaw)
call RecPrt('Shift',' ',Shift(:,iFirst),nQQ,nRaw)
#endif

call mma_allocate(ETemp,nRaw,nSet,Label='ETemp')
call mma_allocate(Temp,nQQ,nRaw,nSet,Label='Temp')
call Prepare_Kriging(ETemp,Temp,nRaw,nQQ,iFirst)
call Setup_Kriging(nRaw,nQQ,qInt(:,iFirst),Temp,ETemp,Hessian_HMF=Hessian(:,:,1))
call mma_deallocate(ETemp)
call mma_deallocate(Temp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the fraction value to be used to define the restricted
! variance threshold.

! Let the accepted variance be set as a fraction of the
! largest component in the gradient.

tmp = Zero
do i=1,size(Gx,2)
  do j=1,3
    tmp = max(tmp,abs(Gx(j,i,iter)))
  end do
end do

! Temporary code until we have figured out this for constrained
! optimizations.
! (Note that for constrained optimizations, Beta is changed again
!  in con_opt, once the constrained gradient is available)

Beta_Disp = max(Beta_Disp_Min,tmp*Beta_Disp_Seed)
Beta = min(1.0e3_wp*GNrm(iter),Beta_Seed)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Start the Kriging loop.   The micro iterations

Not_Converged = .true.
Step_Trunc = 'N'  ! not defined
do while (Not_Converged)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  kIter = iterAI
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Do iterAI: ',iterAI
  write(u6,*) 'Beta=',Beta
# endif

  ! Compute the Kriging Hessian

  First_MicroIteration = iterAI == iter
  if (Kriging_Hessian) then
    call Hessian_Kriging_Layer(qInt(:,iterAI),Hessian,nQQ)
    if ((nSet > 3) .and. NADC) Hessian(:,:,1) = Half*(Hessian(:,:,1)+Hessian(:,:,2))
    call Put_dArray('Hss_Q',Hessian(:,:,1),nQQ**2)
    ! Make fixhess treat the Hessian as if it was analytic.
    call Put_iScalar('HessIter',iterAI)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the updated structure.

  nWndw_ = nWndw/2+(iterAI-iter)

  call Update_inner(iterAI,Beta,Beta_Disp,Step_Trunc,nWndw_,kIter,Kriging_Hessian,qBeta,iOpt_RS,First_MicroIteration,iter, &
                    qBeta_Disp,.false.)

# ifdef _DEBUGPRINT_
  write(u6,*) 'After Update_inner: Step_Trunc=',Step_Trunc
  call RecPrt('New Coord',' ',qInt,nQQ,iterAI+1)
  call RecPrt('dqInt',' ',dqInt,nQQ,iterAI)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !   Save initial gradient
  !   This is for the case the gradient is already converged, we want
  !   the micro iterations to still reduce the gradient
  !
  !   Save the initial error in the constraints,
  !   this is the value that should be used for global convergence
  !
  !   Try to enforce convergence in the constraints if the initial
  !   gradient is close to convergence

  if (First_Microiteration) then
    FAbs_Ini = GNrm(iterAI)/sqrt(real(nQQ,kind=wp))
    if (nLambda > 0) then
      if (FAbs_ini <= Ten*ThrGrd) CheckCons = .true.
      call Qpg_dScalar('Max error',Found)
      if (Found) call Get_dScalar('Max error',MaxErr_Ini)
    end if
  end if

  ! Attempt to break oscillations

  if (IterK > 2) then
    dqdq1 = ddot_(nQQ,Shift(:,iterAI),1,Shift(:,iterAI-1),1)
    ! If the overlap between three consecutive steps is negative,
    ! cut down the step in half.
    if (dqdq1 < Zero) then
      dqdq2 = ddot_(nQQ,Shift(:,iterAI-1),1,Shift(:,iterAI-2),1)
      dqdq3 = ddot_(nQQ,Shift(:,iterAI-2),1,Shift(:,iterAI-3),1)
      if ((dqdq2 < Zero) .and. (dqdq3 < Zero)) then
        Shift(:,iterAI) = Half*Shift(:,iterAI)
        qInt(:,iterAI+1) = qInt(:,iterAI)+Shift(:,iterAI)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Oscillation detected. Step cut down in half!'
        call RecPrt('New Coord',' ',qInt,nQQ,iterAI+1)
#       endif
      end if
    end if
  end if

  ! Transform the new internal coordinates to Cartesians
  ! (this updates the new internal coordinates, in case they were
  ! not totally consistent)

  Error = (iterK >= 1)
  iRef_Save = iRef
  iRef = iter ! Set the reference geometry
  call NewCar_Kriging(iterAI,.true.,Error)
  iRef = iRef_Save
# ifdef _DEBUGPRINT_
  call RecPrt('New Coord (after NewCar)','',qInt,nQQ,iterAI+1)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! In case of a constrained optimization GrdMax and GNrm have been
  ! updated to reflect the contribution from the constraints.
  ! GrdMax, however, is a scalar and we need to save it such that
  ! it has the correct value on exit. That is, the value of iteration
  ! iter.

  if (iterAI == iter) then
    GrdMax_Save = GrdMax
    GrdLbl_Save = GrdLbl
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! During the micro iterations we need to set GNrm to some
  ! reasonable value. Normally it is set by the force routine or
  ! the con_opt routine.

  if ((nLambda == 0) .and. (iterAI > iter)) GNrm(iterAI) = sqrt(DDot_(nQQ,dqInt(:,iterAI),1,dqInt(:,iterAI),1))
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the step length from the last ab initio point
  ! to the most recent kriging point.

  dqdq = Zero
  do iInter=1,nQQ
    dqdq = dqdq+(qInt(iInter,iterAI+1)-qInt(iInter,iter))**2
  end do
  if ((iterK == 0) .and. (Step_trunc == '*')) dqdq = qBeta**2

  ! In case of error during conversion to Cartesians

  if (Error) then
    Step_Trunc = '@'
    exit
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the energy and gradient according to the
  ! surrogate model for the new coordinates.

  call Kriging_Update(nQQ,iterAI+1,qInt(:,iterAI+1),E_Disp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  iterK = iterK+1
  iterAI = iterAI+1
  dEner = Energy(iterAI)-Energy(iterAI-1)
# ifdef _DEBUGPRINT_
  call RecPrt('qInt(x):',' ',qInt,nQQ,iterAI)
  call RecPrt('Ener(x):',' ',Energy,1,iterAI)
  call RecPrt('dqInt(x):',' ',dqInt,nQQ,iterAI)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Check on convergence criteria.

  if (ThrEne > Zero) then
    ! Convergence on energy criteria.
    Not_Converged = abs(dEner) >= ThrEne
    Not_Converged = Not_Converged .and. (dqdq < qBeta**2)
  else
    ! Use standard convergence criteria
    FAbs = GNrm(iterAI-1)/sqrt(real(nQQ,kind=wp))
    RMS = sqrt(DDot_(nQQ,Shift(:,iterAI-1),1,Shift(:,iterAI-1),1)/real(nQQ,kind=wp))
    RMSMx = Zero
    do iInter=1,nQQ
      RMSMx = max(RMSMx,abs(Shift(iInter,iterAI-1)))
    end do
    GrdMx = GrdMax
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'iter=',iterAI-1
    write(u6,*) 'FAbs=',FAbs
    write(u6,*) 'GrdMx=',GrdMx
    write(u6,*) 'RMS=',RMS
    write(u6,*) 'RMSMx=',RMSMx
    write(u6,*) FAbs > min(ThrGrd,FAbs_ini)
    write(u6,*) GrdMx > ThrGrd*OneHalf
    write(u6,*) RMS > ThrGrd*Four
    write(u6,*) RMSMx > ThrGrd*Six
    write(u6,*) 'Step_Trunc=',Step_Trunc
#   endif
    ! Ensure that the initial gradient is reduced,
    ! except in the last micro iteration
    if (iterK < Max_MicroIterations) then
      Not_Converged = FAbs > min(ThrGrd,FAbs_ini)
    else
      Not_Converged = FAbs > ThrGrd
    end if
    Not_Converged = Not_Converged .or. (GrdMx > ThrGrd*OneHalf)
    Not_Converged = Not_Converged .or. (RMS > ThrGrd*Four)
    Not_Converged = Not_Converged .or. (RMSMx > ThrGrd*Six)
    Not_Converged = Not_Converged .or. (Step_Trunc /= ' ')
    if (CheckCons .and. (.not. Not_Converged) .and. (nLambda > 0) .and. (iterK < Max_Microiterations)) then
      call Qpg_dScalar('Max error',Found)
      if (Found) then
        call Get_dScalar('Max error',MaxErr)
        if (MaxErr > ThrCons) Not_Converged = .true.
      end if
    end if
  end if
  if (Step_Trunc == '.') Step_Trunc = ' '
  ! Check total displacement from initial structure
  if (Force_RS) then
    call mma_allocate(Temp,3,size(Cx,2),1,Label='Temp')
    Temp(:,:,1) = Cx(:,:,iterAI)-Cx(:,:,iter)
    RMS = sqrt(DDot_(3*size(Cx,2),Temp,1,Temp,1)/real(3*size(Cx,2),kind=wp))
    if (RMS > (Three*Beta)) Step_trunc = '*'
    call mma_deAllocate(Temp)
  end if
  if (Not_Converged .and. (Step_trunc == ' ') .and. (iterK >= Max_Microiterations)) Step_trunc = '#'
# ifdef _DEBUGPRINT_
  write(u6,*) 'Not_Converged=',Not_Converged
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If the step restriction is invoked or there is some other
  ! reason signaled by Step_Trunc, terminate anyhow.

  if (Step_trunc /= ' ') Not_Converged = .false.

  ! If RS rather than RV do not micro iterate

  if (iOpt_RS == 0) Not_Converged = .false.
  ! Not_Converged=.False. ! Force single iteration scheme.
# ifdef _DEBUGPRINT_
  write(u6,*) 'Not_Converged=',Not_Converged
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !   End of the micro iteration loop

end do  ! Do While (Not_Converged)

! Change label of updating method

UpMeth = 'RVO   '
write(UpMeth(4:6),'(I3)') iterK

#ifdef _OVERSHOOT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Attempt overshooting

call mma_allocate(Step_k,nQQ,2)
Step_k(:,2) = qInt(:,iterAI)-qInt(:,iter)
!call RecPrt('q(iter+1)-q(f)','',Step_k(:,2),nQQ,1)
if (iter > 1) then
  Step_k(:,1) = qInt(:,iter)-qInt(:,iter-1)
  !call RecPrt('q(f)-q(f-1)','',Step_k(:,1),nQQ,1)
  dsds = dDot_(nQQ,Step_k(:,1),1,Step_k(:,2),1)
  dsds = dsds/sqrt(ddot_(nQQ,Step_k(:,1),1,Step_k(:,1),1))
  dsds = dsds/sqrt(ddot_(nQQ,Step_k(:,2),1,Step_k(:,2),1))
  write(u6,*) 'dsds = ',dsds
else
  dsds = Zero
end if
dsds_min = 0.9_wp
if ((dsds > dsds_min) .and. (Step_Trunc == ' ')) then
  do Max_OS=9,0,-1
    OS_Factor = Max_OS*((dsds-dsds_min)/(One-dsds_min))**4
    qInt(:,iterAI) = qInt(:,iter)+(One+OS_Factor)*Step_k(:,2)
    call Energy_Kriging_Layer(qInt(:,iterAI),OS_Energy,nQQ)
    call Dispersion_Kriging_Layer(qInt(:,iterAI),OS_Disp,nQQ)
    write(u6,*) 'Max_OS=',Max_OS
    write(u6,*) OS_Disp(1),E_Disp,Beta_Disp
    if ((OS_Disp(1) > E_Disp) .and. (OS_Disp(1) < Beta_Disp)) then
      Shift(:,iterAI-1) = Shift(:,iterAI-1)+OS_Factor*Step_k(:,2)
      iRef_Save = iRef
      iRef = iter ! Set the reference geometry
      call NewCar_Kriging(iterAI-1,.true.,Error)
      iRef = iRef_Save
      Energy(iterAI) = OS_Energy(1)
      if (Max_OS > 0) then
        if (UpMeth(4:4) /= ' ') UpMeth(5:6) = '**'
        UpMeth(4:4) = '+'
      end if
      exit
    end if
  end do
end if
call mma_deallocate(Step_k)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Save the optimized kriging coordinates as the coordinates
! for the next macro iteration.

qInt(:,iter+1) = qInt(:,iterAI)
Cx(:,:,iter+1) = Cx(:,:,iterAI)

! Update the shift vector

Shift(:,iter) = qInt(:,iter+1)-qInt(:,iter)

! Update the predicted energy change

E_Delta = Energy(iterAI)-Energy(iter)

call MxLbls(nQQ,dqInt(:,iter),Shift(:,iter),Lbl)

! Stick in the correct value for GrdMax, which might contain a
! contribution due to constraints.

GrdMax = GrdMax_Save
GrdLbl = GrdLbl_Save
#ifdef _DEBUGPRINT_
call RecPrt('qInt(3):',' ',qInt,nQQ,iter+1)
call RecPrt('Shift:',' ',Shift,nQQ,iter)
write(u6,*) 'UpMeth=',UpMeth
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocating memory used by Kriging

call mma_deallocate(Hessian)
call Finish_Kriging()
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the shift in internal coordinate basis.

#ifdef _DEBUGPRINT_
call RecPrt('Shifts in internal coordinate basis / au or rad',' ',Shift,nQQ,Iter)
call RecPrt('qInt in internal coordinate basis / au or rad',' ',qInt,nQQ,Iter+1)
#endif

! Remove unneeded fields from the runfile
Dummy(1) = -Zero
call Put_dArray('BMxOld',Dummy(1),0)
call Put_dArray('TROld',Dummy(1),0)
! Restore previous values
call Put_iScalar('HessIter',HessIter)
if (MaxErr_Ini >= Zero) call Put_dScalar('Max error',MaxErr_Ini)

return

end subroutine Update_Kriging
