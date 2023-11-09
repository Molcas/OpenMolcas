!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Convrg(iter,kIter,nInter,iStop,MxItr,mIntEff,mTtAtm,GoOn,Step_Trunc,Just_Frequencies)

use Slapaf_Info, only: Analytic_Hessian, ApproxNADC, Baker, Coor, Cx, dqInt, E_Delta, EDiffZero, eMEPTest, Energy, FindTS, GNrm, &
                       GrdMax, Gx, HUpMet, iNeg, iOptC, Lbl, MaxItr, MEP, NADC, nLambda, nMEP, Numerical, qInt, rMEP, Shift, &
                       SlStop, ThrCons, ThrEne, ThrGrd, ThrMEP
use Chkpnt, only: Chkpnt_update_MEP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Four, Six, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iter, kIter, nInter, MxItr, mIntEff, mTtAtm
integer(kind=iwp), intent(out) :: iStop
logical(kind=iwp), intent(in) :: GoOn, Just_Frequencies
character, intent(in) :: Step_Trunc
#include "print.fh"
#include "warnings.h"
integer(kind=iwp) :: i, iFile, iMEP, iOff_Iter, iPrint, IRC, iRout, iSaddle, iSaddle_p, iSaddle_r, iter_S, j, jSaddle, kkIter, Lu, &
                     LuInput, nAtom, nBackward, nConst, nForward, nIRC, nSaddle, nSaddle_Max
real(kind=wp) :: CumLen, E, E0, E1, E2, E_Prod, E_Reac, echng, eDiffMEP, Fabs, Maxed, MaxErr, prevDist, rDum(1,1,1,1), refDist, &
                 ResGrad, RMS, RMSMax, Thr1, Thr2, Thr3, Thr4, Thr5, Val1, Val2, Val3, Val4, Val5
logical(kind=iwp) :: BadConstraint, Conv1, Conv2, ConvTmp, eTest, Found, IRCRestart, Last_Energy, Saddle, Terminate, TSReg, TurnBack
character(len=80) :: Point_Desc
character(len=16) :: MEP_Text, StdIn
character(len=5) :: ConLbl(5)
character(len=8) :: Temp
integer(kind=iwp), allocatable :: Information(:)
real(kind=wp), allocatable :: C_IRC(:,:,:), C_P(:,:), C_R(:,:), C_S(:,:,:), Coor1(:,:), Coor2(:,:), Cu_MEP(:), E_IRC(:), E_MEP(:), &
                              E_P(:), E_R(:), E_S(:), G_IRC(:,:,:), G_MEP(:,:,:), G_P(:,:), G_R(:,:), G_S(:,:,:), L_MEP(:), Tmp(:)
real(kind=wp), allocatable, target :: C_MEP(:,:,:), OfRef(:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
nAtom = size(Cx,2)
TSReg = btest(iOptC,13)
!                                                                      *
!***********************************************************************
!                                                                      *
Lu = u6
nSaddle_Max = 100
iRout = 116
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt('Convrg: Energy',' ',Energy,1,iter)
  call RecPrt('Convrg: Shift',' ',Shift,nInter,iter)
  call RecPrt('Convrg: qInt',' ',qInt,nInter,iter+1)
  call RecPrt('Convrg: dqInt',' ',dqInt,nInter,iter)
  call RecPrt('Convrg: Cx',' ',Cx,3*nAtom,iter+1)
  call RecPrt('Convrg: Gx',' ',Gx,3*nAtom,iter+1)
end if

call Get_iScalar('Saddle Iter',iter_S)
if (iter_S == 0) then
  iter_S = 1
  call Put_iScalar('Saddle Iter',iter_S)
  call f_Inquire('RUNFILE2',Found)
  if (Found) then
    call NameRun('RUNFILE2')
    call Put_iScalar('Saddle Iter',iter_S)
    call NameRun('#Pop')
  end if
end if
Temp = ' '
if (Analytic_Hessian) then
  if ((HUPMET == '  No  ') .or. (HUPMET == ' None ')) then
    !Temp = 'Analytic'
    Temp = 'Computed'
  else
    Temp(1:6) = HUPMET(1:6)
  end if
else
  Temp(1:6) = HUPMET(1:6)
end if

call mma_allocate(Coor1,3,mTtAtm,Label='Coor1')
call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')
call AtmLst(Cx(:,:,iter),nAtom,Coor1,mTtAtm)
call AtmLst(Cx(:,:,iter+1),nAtom,Coor2,mTtAtm)
call OptRMS_Slapaf(Coor1,Coor2,mTtAtm,RMS,RMSMax)
call mma_deallocate(Coor1)
call mma_deallocate(Coor2)

if ((kIter /= iter) .and. (kIter == 1)) then
  Fabs = GNrm(kiter)
  E = Energy(kiter)
else
  Fabs = GNrm(iter)
  E = Energy(iter)
end if
Fabs = max(Zero,Fabs)
E0 = E+E_delta
Energy(iter+1) = E0
Gx(:,:,iter+1) = Zero
if (kiter == 1) then
  eChng = Zero
else
  if ((kIter /= iter) .and. (kIter == 2)) then
    eChng = Energy(iter)-Energy(1)
  else
    eChng = Energy(iter)-Energy(iter-1)
  end if
end if

eDiffMEP = Zero
if (MEP .or. rMEP) then
  Saddle = .false.
  iMEP = 0
  call Qpg_iScalar('nMEP',Found)
  if (Found) call Get_iScalar('nMEP',iMEP)
  if (iMEP == 0) then
    iOff_Iter = 0
    call Put_iScalar('iOff_Iter',iOff_Iter)
  else
    call Get_iScalar('iOff_Iter',iOff_Iter)
  end if
else
  iOff_Iter = 0
  iSaddle = 0
  call qpg_dArray('Saddle',Saddle,nSaddle)
  Saddle = Saddle .and. (.not. Just_Frequencies)
  if (Saddle) then
    call Qpg_iScalar('nMEP',Found)
    if (Found) call Get_iScalar('nMEP',iSaddle)
    call Get_iScalar('iOff_Iter',iOff_Iter)
  end if
end if

! Convergence criteria
!
! Too many iterations
!
! or
!
! 1) a la Baker
!      Abs(GrdMax) < ThrGrd
!      and
!      Abs(eChng) < ThrEne or RMSMax < ThrGrd
!
! 2) a la Gaussian
!      Abs(Fabs/Sqrt(mIntEff)) < ThrGrd
!      and
!      Abs(GrdMax) < ThrGrd*1.5
!      and
!      ((RMS < ThrGrd*4
!        and
!        RMSMax < ThrGrd*6)
!       or
!       Abs(eChng) < ThrEne)

if (Baker) then
  Val1 = abs(eChng)
  Thr1 = ThrEne
  if (kIter <= 1) then
    ConLbl(1) = ' --- '
  else if (Val1 < Thr1) then
    ConLbl(1) = ' Yes '
  else
    ConLbl(1) = ' No  '
  end if
  Val2 = RMSMax
  Thr2 = ThrGrd
  if (Val2 < Thr2) then
    if (Step_Trunc /= ' ') then
      ConLbl(2) = ' No *'
    else
      ConLbl(2) = ' Yes '
    end if
  else
    ConLbl(2) = ' No  '
  end if
  Val3 = abs(GrdMax)
  Thr3 = ThrGrd
  if (Val3 < Thr3) then
    ConLbl(3) = ' Yes '
  else
    ConLbl(3) = ' No  '
  end if
  Conv1 = (Val1 < Thr1) .and. (kIter > 1)
  Conv1 = Conv1 .or. ((Val2 < Thr2) .and. (Step_Trunc == ' '))
  Conv1 = Conv1 .and. (Val3 < Thr3)
else
  Val2 = abs(Fabs/sqrt(real(mIntEff,kind=wp)))
  Thr2 = ThrGrd
  Conv1 = Val2 < Thr2
  if (Conv1) then
    ConLbl(2) = ' Yes '
  else
    ConLbl(2) = ' No  '
  end if
  Conv1 = Conv1 .and. (abs(GrdMax) < ThrGrd*1.5_wp)
  Val4 = abs(GrdMax)
  Thr4 = ThrGrd*1.5_wp
  ConvTmp = Val4 < Thr4
  Conv1 = Conv1 .and. ConvTmp
  if (ConvTmp) then
    ConLbl(4) = ' Yes '
  else
    ConLbl(4) = ' No  '
  end if
  Conv2 = (RMS < ThrGrd*Four) .and. (Step_Trunc == ' ')
  Val1 = RMS
  Thr1 = ThrGrd*Four
  ConvTmp = Val1 < Thr1
  Conv2 = ConvTmp .and. (Step_Trunc == ' ')
  if (ConvTmp) then
    if (Step_Trunc /= ' ') then
      ConLbl(1) = ' No *'
    else
      ConLbl(1) = ' Yes '
    end if
  else
    ConLbl(1) = ' No  '
  end if
  Val3 = RMSMax
  Thr3 = ThrGrd*Six
  ConvTmp = Val3 < Thr3
  Conv2 = Conv2 .and. ConvTmp
  if (ConvTmp) then
    if (Step_Trunc /= ' ') then
      ConLbl(3) = ' No *'
    else
      ConLbl(3) = ' Yes '
    end if
  else
    ConLbl(3) = ' No  '
  end if
  Val5 = abs(eChng)
  Thr5 = ThrEne
  ConvTmp = (Val5 < Thr5) .and. (kIter > 1)
  if (ConvTmp) then
    ConLbl(5) = ' Yes '
  else
    if (kIter > 1) then
      ConLbl(5) = ' No  '
    else
      ConLbl(5) = ' --- '
    end if
  end if
  Conv2 = Conv2 .or. ConvTmp
  Conv1 = Conv1 .and. Conv2
end if

SlStop = Conv1 .or. (kIter-iOff_Iter >= MxItr) ! CGG
iStop = 1
if (kIter-iOff_Iter >= MxItr) iStop = 16     ! CGG
if (Conv1 .or. Just_Frequencies) iStop = 0

if (GoOn) then
  SlStop = .false.
  iStop = 1
else
  if (Just_Frequencies) SlStop = .true.
  nPrint(52) = nPrint(52)+1
  nPrint(54) = nPrint(54)+1
  iPrint = iPrint+1
  nPrint(53) = nPrint(53)+1
end if
if (.not. Just_Frequencies) call SlStatus(kIter-iOff_Iter,E,Fabs,E0,MaxItr-1,eChng,Temp,Step_Trunc,.not. Numerical)

if (Baker) then
  if (iPrint >= 5) then
    write(Lu,'(A)') '                +----------------------------------+'
    write(Lu,'(A)') '                +  Value      Threshold Converged? +'
    write(Lu,'(A)') '+---------------+----------------------------------+'
    write(Lu,4) '+ Max. gradient +',Val3,' ',Thr3,'    ',ConLbl(3),'  +'
    write(Lu,'(A)') '+---------------+----------------------------------+'
    write(Lu,4) '+ Max. disp.    +',Val2,' ',Thr2,'    ',ConLbl(2),'  +'
    write(Lu,'(A)') '+---------------+----------------------------------+'
    write(Lu,4) '+ Energy diff.  +',Val1,' ',Thr1,'    ',ConLbl(1),'  +'
    write(Lu,'(A)') '+---------------+----------------------------------+'
  end if
else
  if (iPrint >= 5) then
    write(Lu,'(A)') '       +----------------------------------+----------------------------------+'
    write(Lu,'(A)') '       +    Cartesian Displacements       +    Gradient in internals         +'
    write(Lu,'(A)') '       +  Value      Threshold Converged? +  Value      Threshold Converged? +'
    write(Lu,'(A)') ' +-----+----------------------------------+----------------------------------+'
    write(Lu,5) ' + RMS +',Val1,' ',Thr1,'    ',ConLbl(1),'  +',Val2,' ',Thr2,'    ',ConLbl(2),'  +'
    write(Lu,'(A)') ' +-----+----------------------------------+----------------------------------+'
    write(Lu,5) ' + Max +',Val3,' ',Thr3,'    ',ConLbl(3),'  +',Val4,' ',Thr4,'    ',ConLbl(4),'  +'
    write(Lu,'(A)') ' +-----+----------------------------------+----------------------------------+'
    if (ThrEne > Zero) then
      write(Lu,5) ' + dE  +',Val5,' ',Thr5,'    ',ConLbl(5),'  +'
      write(Lu,'(A)') ' +-----+----------------------------------+'
    end if
    write(Lu,*)
  end if
end if

if (SlStop .and. Conv1) then
  call Qpg_dScalar('Max error',Found)
  if (Found) call Get_dScalar('Max error',MaxErr)
  if (MaxErr > ThrCons) then
    iStop = 1
    Conv1 = .false.
    SlStop = .false.
    write(Lu,'(A,ES11.4)') 'Maximum constraint error: ',MaxErr
    write(Lu,*)
  end if
end if

nConst = 0
if (iNeg(1) == 0) then
  if (EDiffZero) then
    if (NADC) then
      nConst = 2
    else
      nConst = 1
    end if
    Point_Desc = 'Minimum Energy Crossing Point Structure'
  else
    Point_Desc = 'Minimum Structure'
  end if
else if (iNeg(1) == 1) then
  Point_Desc = 'Transition State Structure'
else
  Point_Desc = 'Higher Order Saddle Point Structure'
end if
if (nLambda > nConst) Point_Desc = 'Constrained '//trim(Point_Desc)
if (iPrint >= 5) then
  if (SlStop) then
    if (Conv1) then
      write(Lu,'(A,I3,A)') ' Geometry is converged in ',kIter-iOff_iter,' iterations to a '//trim(Point_Desc)
    else
      write(Lu,'(A)') ' No convergence after max iterations'
      if (Lu /= u6) write(u6,'(/A)') ' No convergence after max iterations'
    end if
  else
    write(Lu,'(A)') ' Convergence not reached yet!'
  end if
end if
if (FindTS .and. SlStop .and. Conv1) then
  if (.not. TSReg) then
    if (iPrint >= 5) then
      write(Lu,*)
      write(Lu,'(A)') ' FindTS was requested, but the TS regime was not reached.'
      write(Lu,'(A)') ' The converged structure is probably not the desired TS.'
    end if
    iStop = 16
  end if
end if
!if (iPrint == 7) then
!  write(Lu,*)
!  write(Lu,'(A)') '********************************* Geometry Statistics for Geometry Optimization '// &
!                  '*********************************'
!  nPrint(118) = 7
!  call List(' Internal coordinates ',Lbl,qInt,nInter,iter+1)
!  call List(' Internal forces    ',Lbl,dqInt,nInter,iter)
!end if

! The energy change should not be too large
Maxed = 1.0e2_wp
if (abs(E_Delta) > Maxed) then
  write(u6,*) 'The predicted energy change is too large: ',E_Delta
  write(u6,'(A)') ' This can''t be right!'
  write(u6,'(A)') ' This job will be terminated.'
  iStop = 8
  SlStop = .true.
end if
if (iPrint >= 5) then
  write(Lu,*)
  write(Lu,'(A)') '********************************************************'// &
                  '*********************************************************'
  write(Lu,'(A)') '********************************************************'// &
                  '*********************************************************'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
Terminate = .false.
IRCRestart = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Write summary of conical intersection characterization data

if (Conv1 .and. NADC .and. EDiffZero) then
  if (iPrint >= 5) then
    if (.not. ApproxNADC) call CI_Summary(Lu)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Book keeping for Saddle optimization for a TS. Note that this will
! have to be done on both of the runfiles!

if (Conv1 .and. Saddle) then

  ! Here if a macro iteration in the Saddle TS optimization is completed.

  !ENew = Energy(iter)+E_Delta
  call mma_allocate(Tmp,nSaddle,Label='Tmp')

  ! Store the info for later generation of MOLDEN formated files

  call Get_dArray('Saddle',Tmp,nSaddle)
  E_Reac = Tmp(6*nAtom+1)
  E_Prod = Tmp(6*nAtom+2)

  call mma_allocate(E_S,nSaddle_Max,Label='E_S')
  call mma_allocate(C_S,3,nAtom,nSaddle_Max,Label='C_S')
  call mma_allocate(G_S,3,nAtom,nSaddle_Max,Label='G_S')
  if (iSaddle == 0) then

    ! Initiate with data from the starting points

    E_S(:) = Zero
    C_S(:,:,:) = Zero
    G_S(:,:,:) = Zero
    iSaddle = 1
    if (E_Reac <= E_Prod) then
      E_S(iSaddle) = E_Reac
      C_S(:,:,iSaddle) = reshape(Tmp(1:3*nAtom),[3,nAtom])
    else
      E_S(iSaddle) = E_Prod
      C_S(:,:,iSaddle) = reshape(Tmp(3*nAtom+1:6*nAtom),[3,nAtom])
    end if

  else

    call Get_dArray('MEP-Energies',E_S,nSaddle_Max)
    call Get_dArray('MEP-Coor',C_S,3*nAtom*nSaddle_Max)
    call Get_dArray('MEP-Grad',G_S,3*nAtom*nSaddle_Max)

  end if

  ! Add the new data

  iSaddle = iSaddle+1
  E_S(iSaddle) = Energy(iter)
  C_S(:,:,iSaddle) = Cx(:,:,iter)
  G_S(:,:,iSaddle) = Gx(:,:,iter)

  ! Put data on RUNFILE

  call Put_dArray('MEP-Energies',E_S,nSaddle_Max)
  call Put_dArray('MEP-Coor',C_S,3*nAtom*nSaddle_Max)
  call Put_dArray('MEP-Grad',G_S,3*nAtom*nSaddle_Max)
  call Put_iScalar('nMEP',iSaddle)

  call mma_deallocate(G_S)
  call mma_deallocate(C_S)
  call mma_deallocate(E_S)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now update the "Saddle" field on both runfiles.

  do iFile=1,2

    if (iFile == 1) then
      call NameRun('RUNREAC')
    else
      call NameRun('RUNPROD')
    end if

    ! Update info on the runfile.

    call Get_dArray('Saddle',Tmp,nSaddle)
    E1 = Tmp(6*nAtom+1)
    E2 = Tmp(6*nAtom+2)
    !write(u6,*) 'ENew=',ENew
    !write(u6,*) 'E1,E2=',E1,E2
    if (E1 <= E2) then
      !write(u6,*) 'Update reactant'
      Tmp(6*nAtom+1) = Energy(iter)
      E1 = Energy(iter)
      Tmp(1:3*nAtom) = reshape(Cx(:,:,iter),[3*nAtom])
    else
      !write(u6,*) 'Update product'
      Tmp(6*nAtom+2) = Energy(iter)
      E2 = Energy(iter)
      Tmp(3*nAtom+1:6*nAtom) = reshape(Cx(:,:,iter),[3*nAtom])
    end if
    ! Set flag that seward should process the info! This should not
    ! be done for the final macro iteration.
    if (.not. FindTS) Tmp(6*nAtom+5) = One
    call Put_dArray('Saddle',Tmp,nSaddle)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Converged or not, create the saddle.molden file
  ! after each macro iteration

  call NameRun('RUNREAC')
  call mma_allocate(E_r,nSaddle_Max,Label='E_r')
  call mma_allocate(C_r,3*nAtom,nSaddle_Max,Label='C_r')
  call mma_allocate(G_r,3*nAtom,nSaddle_Max,Label='G_r')
  call Qpg_iScalar('nMEP',Found)
  if (Found) then
    call Get_dArray('MEP-Energies',E_r,nSaddle_Max)
    call Get_dArray('MEP-Coor',C_r,3*nAtom*nSaddle_Max)
    call Get_dArray('MEP-Grad',G_r,3*nAtom*nSaddle_Max)
    call Get_iScalar('nMEP',iSaddle_r)
  else
    E_r(1) = E_Reac
    C_r(:,1) = Tmp(1:3*nAtom)
    G_r(:,:) = Zero
    iSaddle_r = 1
  end if

  call NameRun('RUNPROD')
  call mma_allocate(E_p,nSaddle_Max,Label='E_p')
  call mma_allocate(C_p,3*nAtom,nSaddle_Max,Label='C_p')
  call mma_allocate(G_p,3*nAtom,nSaddle_Max,Label='G_p')
  call Qpg_iScalar('nMEP',Found)
  if (Found) then
    call Get_dArray('MEP-Energies',E_p,nSaddle_Max)
    call Get_dArray('MEP-Coor',C_p,3*nAtom*nSaddle_Max)
    call Get_dArray('MEP-Grad',G_p,3*nAtom*nSaddle_Max)
    call Get_iScalar('nMEP',iSaddle_p)
  else
    E_p(1) = E_Prod
    C_p(:,1) = Tmp(3*nAtom+1:6*nAtom)
    G_p(:,:) = Zero
    iSaddle_p = 1
  end if

  ! Merge the two lists

  jSaddle = iSaddle_r
  do iSaddle=iSaddle_p,1,-1
    jSaddle = jSaddle+1
    E_r(jSaddle) = E_p(iSaddle)
    C_r(:,jSaddle) = C_p(:,iSaddle)
    G_r(:,jSaddle) = G_p(:,iSaddle)
  end do
  call mma_deallocate(G_p)
  call mma_deallocate(C_p)
  call mma_deallocate(E_p)

  ! Align the structures sequentially, only for visualization
  ! (gradients are not changed, though)
  ! TODO: Rotate the gradients too

  do iSaddle=1,iSaddle_r+iSaddle_p-1
    call Align(C_r(:,iSaddle+1),C_r(:,iSaddle),nAtom)
  end do

  call Intergeo('MD_SADDLE',E_r,C_r,G_r,nAtom,iSaddle_r+iSaddle_p)

  call mma_deallocate(G_r)
  call mma_deallocate(C_r)
  call mma_deallocate(E_r)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If the Saddle TS optimization is not yet completed set up
  ! data for the next macro iteration.

  if (.not. FindTS) then
    call mma_deallocate(Tmp)
    !write(u6,*) 'Reset $SubProject'

    ! Reset $SubProject for the next macro iteration.

    LuInput = 11
    LuInput = IsFreeUnit(LuInput)
    call StdIn_Name(StdIn)
    call Molcas_Open(LuInput,StdIn)
    if (E1 <= E2) then
      write(LuInput,'(A)') '> EXPORT SubProject=.Reac'
      !write(u6,*) 'SubProject=.Reac'
      call NameRun('RUNREAC')
    else
      write(LuInput,'(A)') '> EXPORT SubProject=.Prod'
      !write(u6,*) 'SubProject=.Prod'
      call NameRun('RUNPROD')
    end if

    ! Signal whether next iteration will be the first in the branch

    call Qpg_iScalar('nMEP',Found)
    if (.not. Found) write(LuInput,'(A)') '> EXPORT SADDLE_FIRST=1'
    write(LuInput,'(A,I3)') '> EXIT ',_RC_CONTINUE_LOOP_
    close(LuInput)

    ! Set flags to request yet another macro iteration.

    Terminate = .false.
    iStop = 6
    SlStop = .false.
  else
    call NameRun('RUNFILE')
    call mma_deallocate(Tmp)
    nSaddle = 0
    call Put_dArray('Saddle',[Zero],nSaddle)
    call Put_iScalar('nMEP',nSaddle)
  end if

  ! Update the active runfile wrt the total number of micro
  ! iterations done in all macro iterations of this branch.

  call NameRun('RUNFILE')
  call Put_iScalar('iOff_Iter',iter)
end if

! Disable first iteration signal right after the first iteration
! (in each branch)

if ((.not. Conv1) .and. Saddle) then
  if ((iter == 1) .and. (iStop == 1)) then
    LuInput = 11
    LuInput = IsFreeUnit(LuInput)
    call StdIn_Name(StdIn)
    call Molcas_Open(LuInput,StdIn)
    write(LuInput,'(A)') '> EXPORT SADDLE_FIRST=0'
    write(LuInput,'(A,I3)') '> EXIT ',_RC_CONTINUE_LOOP_
    close(LuInput)
    iStop = 6
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Book keeping for minimum energy path search

IRC = 0
call Qpg_iScalar('IRC',Found)
if (Found) call Get_iScalar('IRC',IRC)

TurnBack = .false.
if (MEP .or. rMEP) then
  if (Conv1) then

    ! Is this the first iteration or not?

    iMEP = iMEP+1

    ! Save information for the current step

    call mma_allocate(E_MEP,nMEP+1,Label='E_MEP')
    call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
    call mma_allocate(G_MEP,3,nAtom,nMEP+1,Label='G_MEP')
    if (iMEP > 1) then
      call Get_dArray('MEP-Energies',E_MEP,nMEP+1)
      call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
      call Get_dArray('MEP-Grad',G_MEP,3*nAtom*(nMEP+1))
    else
      E_MEP(:) = Zero
      C_MEP(:,:,:) = Zero
      G_MEP(:,:,:) = Zero
      E_MEP(iMEP) = Energy(iOff_iter+1)
      C_MEP(:,:,iMEP) = Cx(:,:,iOff_iter+1)
      G_MEP(:,:,iMEP) = Gx(:,:,iOff_iter+1)
    end if

    E_MEP(iMEP+1) = Energy(iter)
    C_MEP(:,:,iMEP+1) = Cx(:,:,iter)
    G_MEP(:,:,iMEP+1) = Gx(:,:,iter)
    call Put_dArray('MEP-Energies',E_MEP,nMEP+1)
    call Put_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
    call Put_dArray('MEP-Grad',G_MEP,3*nAtom*(nMEP+1))
    call Put_iScalar('nMEP',iMEP)

    ! Save the path so far (energies, coordinates and forces)

    call Intergeo('MD_MEP',E_MEP,C_MEP,G_MEP,nAtom,iMEP+1)

    ! Compute energy difference and RMS between last two structures

    eDiffMEP = E_MEP(iMEP+1)-E_MEP(iMEP)
    call mma_allocate(Coor1,3,mTtAtm,Label='Coor1')
    call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')
    call AtmLst(C_MEP(:,:,iMEP),nAtom,Coor1,mTtAtm)
    call AtmLst(C_MEP(:,:,iMEP+1),nAtom,Coor2,mTtAtm)
    call OptRMS_Slapaf(Coor1,Coor2,mTtAtm,RMS,RMSMax)
    call mma_deallocate(Coor1)
    call mma_deallocate(Coor2)

    call mma_deallocate(G_MEP)
    call mma_deallocate(C_MEP)
    call mma_deallocate(E_MEP)

  else

    ! Test for "turn back", i.e. when the trial structure
    ! is getting too close to the previous converged structure,
    ! this may be an indication of an ill-behaved constraint
    if (iMEP >= 1) then
      call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
      call mma_allocate(OfRef,3,nAtom,Label='OfRef')
      call mma_Allocate(Tmp,3*nAtom,Label='Tmp')
      call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
      OfRef(:,:) = C_MEP(:,:,iMEP+1)

      ! Using hypersphere measure, even with "transverse" MEPs,
      ! this should not be a problem
      call SphInt(Cx(:,:,iter),nAtom,OfRef,.false.,refDist,Tmp,.false.,'dummy   ',rDum,.false.)
      call SphInt(Cx(:,:,iter),nAtom,OfRef,.true.,prevDist,Tmp,.false.,'dummy   ',rDUm,.false.)
      if (prevDist < Half*refDist) then
        TurnBack = .true.
        Conv1 = .true.
        SlStop = .true.
        iStop = 0
        Terminate = .true.
      end if
      call mma_deallocate(Tmp)
      call mma_deallocate(OfRef)
      call mma_deallocate(C_MEP)
    end if

  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! List internal coordinates and gradients

kkIter = iter+1
if (iPrint >= 8) then
  call List(' Internal coordinates ',Lbl,qInt,nInter,kkIter)
  call List(' Internal forces    ',Lbl,dqInt,nInter,iter)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put out the new reference structure and the new starting
! structure to be used for the next MEP point.
! For rMEP keep the reference structure!
! Note that this is done in weighted Cartesian coordinates!

if ((Conv1 .or. (iter == 1)) .and. (MEP .or. rMEP)) then
  if ((iMEP >= 1) .and. (iPrint >= 5)) then
    write(u6,*)
    call CollapseOutput(1,'IRC/Minimum Energy Path Information')
  end if

  ResGrad = huge(ResGrad)
  if (.not. Terminate) then
    BadConstraint = .false.
    call MEP_Dir(Cx,Gx,nAtom,iMEP,iOff_iter,iPrint,IRCRestart,ResGrad,BadConstraint)
    Coor(:,:) = Cx(:,:,iter+1)
    call Put_iScalar('iOff_Iter',iter)
  end if

  if (MEP) then
    if (IRC == 0) then
      MEP_Text = 'MEP'
    else if (IRC == 1) then
      MEP_Text = 'IRC(forward)'
    else
      MEP_Text = 'IRC(backward)'
    end if
  else if (rMEP) then
    MEP_Text = 'rMEP'
  else
    MEP_Text = ''
  end if

  ! Should we terminate or not? Not done on the first iteration.

  if (iMEP > 0) then

    ! Test for energy increase (optionally disabled).
    eTest = eMEPTest .and. (eDiffMEP > Zero)
    if ((MEP .and. eTest) .and. (.not. Terminate)) then
      Terminate = .true.
      if (iPrint >= 5) then
        write(u6,*)
        write(u6,'(A)') ' '//trim(MEP_Text)//'-search terminated due to energy increase!'
        write(u6,*)
      end if
    end if

    ! Test for energy decrease (optionally disabled).
    eTest = eMEPTest .and. (eDiffMEP < Zero)
    if ((rMEP .and. eTest) .and. (.not. Terminate)) then
      Terminate = .true.
      if (iPrint >= 5) then
        write(u6,*)
        write(u6,'(A)') ' '//trim(MEP_Text)//'-search terminated due to energy decrease!'
        write(u6,*)
      end if
    end if

    ! Test for small gradient.
    if ((iMEP > 1) .or. (IRC == 0)) then
      if ((ResGrad < ThrMEP) .and. (.not. Terminate)) then
        Terminate = .true.
        if (iPrint >= 5) then
          write(u6,*)
          write(u6,'(A)') ' '//trim(MEP_Text)//'-search terminated due to small gradient!'
          write(u6,*)
        end if
      end if
    end if

    ! Test for small step.
    if ((RMS < ThrGrd*Four) .and. (.not. Terminate)) then
      Terminate = .true.
      if (iPrint >= 5) then
        write(u6,*)
        write(u6,'(A)') ' '//trim(MEP_Text)//'-search terminated due to small geometry change!'
        write(u6,*)
      end if
    end if

    ! Test for max number of points.
    if ((iMEP >= nMEP) .and. (.not. Terminate)) then
      Terminate = .true.
      if (iPrint >= 5) then
        write(u6,*)
        write(u6,'(A)') ' '//trim(MEP_Text)//'-search terminated due to max number of path points!'
        write(u6,*)
      end if
    end if

    ! Test for constraint misbehavior.
    if ((BadConstraint .and. (.not. Terminate)) .or. TurnBack) then
      Terminate = .true.
      if (iPrint >= 5) then
        write(u6,*)
        write(u6,'(A)') ' '//trim(MEP_Text)//'-search terminated due to problematic constraint!'
        write(u6,*)
      end if
    end if

    ! If IRC reset for backward IRC search.

    if (Terminate .and. (IRC == 1)) IRCRestart = .true.

  end if

  if (Conv1) call Chkpnt_update_MEP(.not. TurnBack,IRCRestart)

  if (Conv1 .and. Terminate) then
    if (IRC /= 0) then
      call mma_allocate(E_MEP,nMEP+1,Label='E_MEP')
      call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
      call mma_allocate(G_MEP,3,nAtom,nMEP+1,Label='G_MEP')
      call Get_dArray('MEP-Energies',E_MEP,nMEP+1)
      call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
      call Get_dArray('MEP-Grad',G_MEP,3*nAtom*(nMEP+1))
      if (IRC == 1) then
        IRCRestart = .true.
        IRC = -1
        call Put_iScalar('IRC',IRC)

        ! Store away data for IRC molden file. Forward part.

        call Put_dArray('IRC-Energies',E_MEP,iMEP+1)
        call Put_dArray('IRC-Coor',C_MEP,3*nAtom*(iMEP+1))
        call Put_dArray('IRC-Grad',G_MEP,3*nAtom*(iMEP+1))
        call Put_dArray('Ref_Geom',Cx,3*nAtom)

        ! Write a temporary file
        ! (will be overwritten when the backward part is done)

        call Intergeo('MD_IRC',E_MEP,C_MEP,G_MEP,nAtom,iMEP+1)

        Terminate = .false.

      else if (IRC == -1) then

        ! Assemble molden file for IRC

        nBackward = iMEP+1
        call qpg_dArray('IRC-Energies',Found,nForward)
        nIRC = nForward+nBackward-1
        call mma_allocate(E_IRC,nIRC,Label='E_IRC')
        call mma_allocate(C_IRC,3,nAtom,nIRC,Label='C_IRC')
        call mma_allocate(G_IRC,3,nAtom,nIRC,Label='G_IRC')

        j = 0
        do i=nBackward,1,-1
          j = j+1
          E_IRC(j) = E_MEP(i)
          C_IRC(:,:,j) = C_MEP(:,:,i)
          G_IRC(:,:,j) = G_MEP(:,:,i)
        end do

        call Get_dArray('IRC-Energies',E_IRC(nBackward),nForward)
        call Get_dArray('IRC-Coor',C_IRC(:,:,nBackward:nIRC),nForward*3*nAtom)
        call Get_dArray('IRC-Grad',G_IRC(:,:,nBackward:nIRC),nForward*3*nAtom)

        call Intergeo('MD_IRC',E_IRC,C_IRC,G_IRC,nAtom,nIRC)

        call mma_deallocate(G_IRC)
        call mma_deallocate(C_IRC)
        call mma_deallocate(E_IRC)
      end if

      call mma_deallocate(G_MEP)
      call mma_deallocate(C_MEP)
      call mma_deallocate(E_MEP)
    end if
  end if

  if (.not. Terminate) then
    iStop = 1
    SlStop = .false.
  end if

  ! Print out the path so far

  if ((iMEP >= 1) .and. (iPrint >= 5)) then
    call mma_allocate(E_MEP,nMEP+1,Label='E_MEP')
    call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
    call mma_allocate(L_MEP,nMEP+1,Label='L_MEP')
    call mma_allocate(Cu_MEP,nMEP+1,Label='Cu_MEP')
    call Get_dArray('MEP-Energies',E_MEP,nMEP+1)
    call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
    call Get_dArray('MEP-Lengths',L_MEP,nMEP+1)
    call Get_dArray('MEP-Curvatures',Cu_MEP,nMEP+1)
    write(u6,*)
    CumLen = Zero
    if (Cu_MEP(1+iMEP) >= Zero) then
      write(u6,*) '         Cumul.'
      write(u6,*) 'Point  Length (bohr)       Energy  Curvature'
      write(u6,*) '--------------------------------------------'
      do i=0,iMEP
        CumLen = CumLen+L_MEP(1+i)
        write(u6,200) i,CumLen,E_MEP(1+i),Cu_MEP(1+i)
      end do
    else
      write(u6,*) '         Cumul.'
      write(u6,*) 'Point  Length (bohr)       Energy'
      write(u6,*) '---------------------------------'
      do i=0,iMEP
        CumLen = CumLen+L_MEP(1+i)
        write(u6,200) i,CumLen,E_MEP(1+i)
      end do
    end if
    if (iPrint > 6) then
      write(u6,*)
      do i=0,iMEP
        call RecPrt(' Coordinates',' ',C_MEP(:,:,i+1),3,nAtom)
      end do
    end if
    call CollapseOutput(0,'IRC/Minimum Energy Path Information')
    write(u6,*)
    call mma_deallocate(E_MEP)
    call mma_deallocate(C_MEP)
    call mma_deallocate(L_MEP)
    call mma_deallocate(Cu_MEP)
  end if

end if

if (IRCRestart) then

  ! Prepare the runfile to start from Scratch

  iMEP = 0
  call Put_iScalar('nMEP',iMEP)
  call mma_allocate(Information,7,Label='Information')
  Information(:) = 0
  Information(1) = -99     ! Deactivate the record
  call Put_iArray('Slapaf Info 1',Information,7)
  call mma_deallocate(Information)
  iOff_Iter = 0
  call Put_iScalar('iOff_Iter',iOff_Iter)

  ! Restore data

  call Put_dScalar('Last Energy',Energy(1))
  call Put_dArray('GRAD',Gx,3*nAtom)
  call Put_dArray('Unique Coordinates',Cx,3*nAtom)
  call Put_Coord_New(Cx,nAtom)
  Coor(:,:) = Cx(:,:,1)
  Cx(:,:,iter+1) = Cx(:,:,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Figure out if the last energy should be computed!

Last_Energy = SlStop .and. (iStop /= 16) .and. (iStop /= 8)
Last_Energy = Last_Energy .and. (.not. MEP) .and. (.not. rMEP)
Last_Energy = Last_Energy .and. (.not. (Numerical .and. (kIter == 1)))
if (Last_Energy) iStop = 2
!                                                                      *
!***********************************************************************
!                                                                      *
! Write geometry file for MOLDEN. Note that this file should not be
! generated if Slapaf is running new geometries for any numerical
! procedure!

if ((.not. Just_Frequencies) .and. (.not. Numerical)) then
  call Write_QMMM(Cx,nAtom,iter)
  call Intergeo('MD_GEO',Energy,Cx,Gx,nAtom,iter+1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

4 format(A,2(ES11.4,A),A,A)
5 format(A,2(2(ES11.4,A),A,A))
200 format(1X,I5,1X,F10.6,1X,F16.8,1X,F10.6)

end subroutine Convrg
