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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  Thomas Bondo Pedersen, Nov. 28, 2012:
!  This code was obtained from S. Ten-no; conditions (license or other
!  terms) are unknown at present.
!  It implements the minimax approximation for the numerical
!  quadrature for Laplace transformation of 1/x. Weights and
!  points are returned in Coeff:
!     Coeff(1) = Omega(1)   Coeff(2) = Alpha(1)
!     Coeff(3) = Omega(2)   Coeff(4) = Alpha(2)
!     etc.
!  Reference: Takatsuka, Ten-no, Hackbusch; JCP 129, 044112 (2008)

subroutine Remez(Vrbse,K_Lap,EMin,EMax,Coeff,Demand,Inf)
!======================================================================*
!----------------------------------------------------------------------*
!                                                                      *
!     Function : Main program of Remez algorithm                       *
!                                                                      *
!                1      K_Lap                                          *
!             ------- := SUM Omega(N)*EXP(-Alpha(N)*DELTA)             *
!              DELTA     N=1                                           *
!                                                                      *
!              DELTA  := EA+EB-EI-EJ                                   *
!                                                                      *
!     Author(s): Akio Takatsuka (2007)                                 *
!                                                                      *
!     K_Lap : The number of quadrature points (K_Lap=1,...,20)         *
!             K_Lap=0 : K_Lap is automatically determined.             *
!                                                                      *
!----------------------------------------------------------------------*
!======================================================================*

use ReMez_mod

implicit real*8(A-H,O-Z)
logical Vrbse
logical Inf
parameter(ItrEnd=50,ZERO=0.0D+00,ONE=1.0D+00,MxList=30)
character*8 Demand
logical Verbose
!tbp  dimension Alpha(20), Omega(20), Xi(40), Coeff(40), V(40), T(40), DD(82), CofOld(40), RList(MxList), IRMax(20), RMin(20), &
!               RMax(20)
real*8 Coeff(40), T(40), DD(82), CofOld(40), RList(MxList), RMin(20), RMax(20)
integer IRMax(20)
data RList/2.0D+00,5.0D+00,1.0D+01,2.0D+01,3.0D+01,4.0D+01,5.0D+01,6.0D+01,7.0D+01,8.0D+01,9.0D+01,1.0D+02,2.0D+02,3.0D+02, &
           4.0D+02,5.0D+02,6.0D+02,7.0D+02,8.0D+02,9.0D+02,1.0D+03,2.0D+03,3.0D+03,4.0D+03,5.0D+03,6.0D+03,7.0D+03,8.0D+03, &
           9.0D+03,1.0D+04/
data RMin/2.0D+00,2.0D+00,2.0D+00,2.0D+00,2.0D+00,2.0D+00,2.0D+00,1.0D+01,1.0D+01,1.0D+01,1.0D+01,1.0D+01,1.0D+01,1.0D+01,1.0D+01, &
          1.0D+01,2.0D+01,4.0D+01,5.0D+01,6.0D+01/
data RMax/8.667D+00,4.154D+01,1.468D+02,4.361D+02,1.154D+03,2.802D+03,6.373D+03,1.375D+04,2.839D+04,5.650D+04,1.090D+05,2.045D+05, &
          3.738D+05,6.692D+05,1.175D+06,2.080D+06,3.450D+06,5.754D+06,9.494D+06,1.550D+07/
data IRMax/3,7,13,16,22,23,27,31,31,31,31,31,31,31,31,31,31,31,31,31/
logical Change, Dbg, SkpRem, StopBA

call Untested('Laplace quadrature generation (subroutine remez)')

Dbg = .false.
SkpRem = .false.
Change = .false.
IMes = -1

Verbose = Vrbse .or. Dbg

call Remez_SetupPrint(Verbose)

write(IW,'(/A)') ' Remez: minimax approximation for Laplace grid '
write(IW,'(A/)') ' ============================================= '

! ===== Check input variables =====

if (EMin > EMax) then
  write(IW,'(/A/)') ' Input values of energy are unsuitable. '
  write(IW,*) 'EMin, EMax',EMin,EMax
  K_Lap = -1
  call Remez_ShutDownPrint(Verbose)
  return
else if ((K_Lap < 0) .or. (K_Lap > 20)) then
  write(IW,'(/A)') ' Input value of K is unsuitable. '
  write(IW,*) 'K_Lap',K_Lap
  K_Lap = -2
  call Remez_ShutDownPrint(Verbose)
  return
end if

! ===== Set initial values =====
!       InitR : initial R

EMinIv = ONE/EMin
R = EMax*EMinIv

J = MxList
do I=1,MxList
  EDiff = R-RList(J)
  if (EDiff > ZERO) then
    InitR = J
    goto 100
  end if
  J = J-1
end do
InitR = 1
100 continue

! ===== Determine K_Lap value (If K_Lap isn't determined.) =====

if (K_Lap == 0) call DfineK(K_Lap,R,InitR,Demand)
IDim = 2*K_Lap

RMax0 = RMax(K_Lap)
write(IW,'(1X,A,I3/)') '# of quadrature points =',K_Lap
write(IW,'(1X,A,F20.10)') 'Maximum R =',RMax0
write(IW,'(A)') ' MO energy'
write(IW,'(4X,A,2X,F14.8)') 'Min =',EMin
write(IW,'(4X,A,2X,F14.8)') 'Max =',EMax
write(IW,'(4X,A,F16.8)') 'R   =',R

! ===== Check the R value =====
!   R <= RMin(K_Lap) or R > RMax(K_Lap)

if (R <= RMin(K_Lap)) then
  SkpRem = .true.
  InitR = 1
  IMes = 0
  RLim = RMin(K_Lap)
else if ((R > RList(IRMax(K_Lap)-1)) .and. (R < RMax(K_Lap))) then
  InitR = IRMax(K_Lap)-1
else if (R > RMax(K_Lap)) then
  SkpRem = .true.
  InitR = IRMax(K_Lap)
  IMes = 1
end if
if (Inf) then
  SkpRem = .true.
  InitR = IRMax(K_Lap)
  IMes = 2
end if
if (SkpRem) then
  write(IW,'(/A)') ' Remez step skipped.'
  if (IMes == 0) write(IW,'(A,F8.3,A)') ' Because R is smaller than',RLim,'.'
  if (IMes == 1) write(IW,'(A,F8.3,A)') ' Because R is larger than',RMax(K_Lap),'.'
  goto 999
end if

! ===== Set initial values (Omega,Alpha,T) =====

call SetExp(K_Lap,InitR,R,RIni,Coeff,T)

if (Dbg) then
  write(IW,'(/A/)') 'Initial data for Remez algorithm'
  do I=1,K_Lap
    Idx = 2*I-1
    write(IW,'(A,I3,A,F21.18,2X,A,I3,A,F20.16)') ' Omega(',I,') = ',Coeff(Idx),'Alpha(',I,') = ',Coeff(Idx+1)
  end do
end if

! ===== Start Remez step =====

Theta = ZERO
Theta2 = ONE

write(IW,'(A,F5.0/)') ' Remez step starts from ',RIni
call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
if (Dbg) then
  write(IW,'(A/)') ' Error output (initial)'
  do I=1,2*(2*K_Lap+1)
    write(IW,*) 'dd',I,DD(I)
  end do
end if
if (StopBA) then
  Change = .true.
  goto 999
end if

! ===== Iteration =====

do Iter=1,ItrEnd
  call DCOPY_(IDim,Coeff,1,CofOld,1)
  call SlvNt2(K_Lap,R,Coeff,T,Theta2,VVMax,StopBA)
  if (StopBA) then
    Change = .true.
    goto 999
  end if
  !call TStat(" Newton2",1)
  call SlvNt1(K_Lap,100,Coeff,T)
  !call TStat(" Newton1",1)

  if (Dbg) then
    write(IW,'(/A,I3)') 'Iter =',Iter
    call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
    do J=1,IDim+1
      write(IW,'(A,I2,A,F20.17,2X,A,F20.14)') ' Error(',J,') = ',DD(J),' X = ',DD(2*K_Lap+1+J)
    end do
    write(IW,*)
  end if

  ! ===== Check convergence =====
  ! CofOld is updated by the mean of Coeff

  do J=1,IDim
    CofOld(J) = abs(CofOld(J)-Coeff(J))
  end do
  call FindAM(IDim+1,DD,DDMax,IMax)
  DCofMx = FindMx(IDim,CofOld)

  if (Dbg) then
    call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
    call FindAM(IDim+1,DD,DDMax,IMax)
    if (Iter == 1) write(IW,'(A,5X,A,15X,A,17X,A)') ' Iter','Max Change','Max DD','Max error'
    write(IW,'(I3,3(2X,E23.15E3),A,I2,A)') Iter,DCofMx,VVMax,DDMax,' (',IMax,')'
  else
    if (Iter == 1) write(IW,'(A,5X,A,15X,A)') ' Iter','Max Change','Max DD'
    write(IW,'(I3,2(2X,E23.15E3))') Iter,DCofMx,VVMax
  end if

  if ((VVMax < 5.0D-16) .and. (DCofMx < 1.0D-05)) goto 777

  !call TStat(" Check T",1)
end do

write(IW,'(A,I3,A)') ' The number of iterations exceeded ',ItrEnd,' .'

! ===== Check the result =====

777 continue
call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
if (StopBA) then
  Change = .true.
  goto 999
end if
call FindAM(IDim+1,DD,DDMax,IMax)
call ChkAcc(K_Lap,InitR,DDMax,R,Change)

! ===== Change parameters =====
! If the convergence isn't good, the larger R's parameters are used.

write(IW,'(/A,I6)') ' Iteration =',Iter
999 continue
if (Change .or. SkpRem) then
  write(IW,'(1X,A)') 'Change!!'
  call SetExp(K_Lap,InitR,R,RIni,Coeff,T)
end if

! ===== Determine Xi =====

!call TStat(" Check T",1)
call SlvNt1(K_Lap,1,Coeff,T)
if (Change)  R = RIni
!call TStat(" Newton1",1)

write(IW,'(/A,F10.4/)') ' Optimized solution in R = ',R
do I=1,K_Lap
  Idx = 2*I-1
  write(IW,'(A,I3,A,F21.18,2X,A,I3,A,F20.17)') ' Omega(',I,') = ',Coeff(Idx),'Alpha(',I,') = ',Coeff(Idx+1)
end do
write(IW,*)
call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
do I=1,IDim+1
  write(IW,'(A,I2,A,F20.17,2X,A,F20.14)') ' Error(',I,') = ',DD(I),' X = ',DD(2*K_Lap+1+I)
end do
write(IW,*)
do I=1,IDim
  write(IW,'(A,I2,A,F20.15)') ' Xi(',I,') = ',T(I)
end do
write(IW,'()')

! ===== Maltiple the scaling parameter (EMinIv = 1/EMin) =====

call DSCAL_(IDim,EMinIv,Coeff,1)
write(IW,'(/A/)') ' Final solution '
do I=1,K_Lap
  Idx = 2*I-1
  write(IW,'(A,I3,A,F21.18,2X,A,I3,A,F20.17)') ' Omega(',I,') = ',Coeff(Idx),'Alpha(',I,') = ',Coeff(Idx+1)
end do

call Remez_ShutDownPrint(Verbose)

return

end subroutine Remez
