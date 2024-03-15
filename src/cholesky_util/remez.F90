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

use ReMez_mod, only: IW
use Constants, only: Zero, One, Two, Ten
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: Vrbse, Inf
integer(kind=iwp), intent(inout) :: K_Lap
real(kind=wp), intent(in) :: EMin, EMax
real(kind=wp), intent(out) :: Coeff(40)
character(len=8), intent(in) :: Demand
integer(kind=iwp) :: I, I_Dim, Idx, IMax, IMes, InitR, Iter, J
real(kind=wp) :: CofOld(40), DCofMx, DD(82), DDMax, EDiff, EMinIv, R, RIni, RLim, RMax0, T(40), Theta, Theta2, VVMax
logical(kind=iwp) :: Change, Conv, Dbg, SkpRem, StopBA, Verbose
integer(kind=iwp), parameter :: IRMax(20) = [3,7,13,16,22,23,27,31,31,31,31,31,31,31,31,31,31,31,31,31], ItrEnd = 50, MxList = 30
real(kind=wp), parameter :: RList(MxList) = [2.0e0_wp,5.0e0_wp,1.0e1_wp,2.0e1_wp,3.0e1_wp,4.0e1_wp,5.0e1_wp,6.0e1_wp,7.0e1_wp, &
                                             8.0e1_wp,9.0e1_wp,1.0e2_wp,2.0e2_wp,3.0e2_wp,4.0e2_wp,5.0e2_wp,6.0e2_wp,7.0e2_wp, &
                                             8.0e2_wp,9.0e2_wp,1.0e3_wp,2.0e3_wp,3.0e3_wp,4.0e3_wp,5.0e3_wp,6.0e3_wp,7.0e3_wp, &
                                             8.0e3_wp,9.0e3_wp,1.0e4_wp], &
                            RMin(20) = [Two,Two,Two,Two,Two,Two,Two,Ten,Ten,Ten,Ten,Ten,Ten,Ten,Ten,Ten,2.0e1_wp,4.0e1_wp, &
                                        5.0e1_wp,6.0e1_wp], &
                            RMax(20) = [8.667e0_wp,4.154e1_wp,1.468e2_wp,4.361e2_wp,1.154e3_wp,2.802e3_wp,6.373e3_wp,1.375e4_wp, &
                                        2.839e4_wp,5.650e4_wp,1.090e5_wp,2.045e5_wp,3.738e5_wp,6.692e5_wp,1.175e6_wp,2.080e6_wp, &
                                        3.450e6_wp,5.754e6_wp,9.494e6_wp,1.550e7_wp]
real(kind=wp), external :: FindMx

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
! InitR : initial R

EMinIv = One/EMin
R = EMax*EMinIv

InitR = 1
J = MxList
do I=1,MxList
  EDiff = R-RList(J)
  if (EDiff > Zero) then
    InitR = J
    exit
  end if
  J = J-1
end do

! ===== Determine K_Lap value (If K_Lap isn't determined.) =====

if (K_Lap == 0) call DfineK(K_Lap,R,InitR,Demand)
I_Dim = 2*K_Lap

RMax0 = RMax(K_Lap)
write(IW,'(1X,A,I3/)') '# of quadrature points =',K_Lap
write(IW,'(1X,A,F20.10)') 'Maximum R =',RMax0
write(IW,'(A)') ' MO energy'
write(IW,'(4X,A,2X,F14.8)') 'Min =',EMin
write(IW,'(4X,A,2X,F14.8)') 'Max =',EMax
write(IW,'(4X,A,F16.8)') 'R   =',R

! ===== Check the R value =====
! R <= RMin(K_Lap) or R > RMax(K_Lap)

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
else

  ! ===== Set initial values (Omega,Alpha,T) =====

  call SetExp(K_Lap,InitR,RIni,Coeff,T)

  if (Dbg) then
    write(IW,'(/A/)') 'Initial data for Remez algorithm'
    do I=1,K_Lap
      Idx = 2*I-1
      write(IW,'(A,I3,A,F21.18,2X,A,I3,A,F20.16)') ' Omega(',I,') = ',Coeff(Idx),'Alpha(',I,') = ',Coeff(Idx+1)
    end do
  end if

  ! ===== Start Remez step =====

  Theta = Zero
  Theta2 = One

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
  else

    ! ===== Iteration =====

    Conv = .false.
    do Iter=1,ItrEnd
      CofOld(1:I_Dim) = Coeff(1:I_Dim)
      call SlvNt2(K_Lap,R,Coeff,T,Theta2,VVMax,StopBA)
      if (StopBA) then
        Change = .true.
        exit
      end if
      !call TStat(" Newton2",1)
      call SlvNt1(K_Lap,100,Coeff,T)
      !call TStat(" Newton1",1)

      if (Dbg) then
        write(IW,'(/A,I3)') 'Iter =',Iter
        call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
        do J=1,I_Dim+1
          write(IW,'(A,I2,A,F20.17,2X,A,F20.14)') ' Error(',J,') = ',DD(J),' X = ',DD(2*K_Lap+1+J)
        end do
        write(IW,*)
      end if

      ! ===== Check convergence =====
      ! CofOld is updated by the mean of Coeff

      CofOld(1:I_Dim) = abs(CofOld(1:I_Dim)-Coeff(1:I_Dim))
      call FindAM(I_Dim+1,DD,DDMax,IMax)
      DCofMx = FindMx(I_Dim,CofOld)

      if (Dbg) then
        call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
        call FindAM(I_Dim+1,DD,DDMax,IMax)
        if (Iter == 1) write(IW,'(A,5X,A,15X,A,17X,A)') ' Iter','Max Change','Max DD','Max error'
        write(IW,'(I3,3(2X,ES23.15E3),A,I2,A)') Iter,DCofMx,VVMax,DDMax,' (',IMax,')'
      else
        if (Iter == 1) write(IW,'(A,5X,A,15X,A)') ' Iter','Max Change','Max DD'
        write(IW,'(I3,2(2X,ES23.15E3))') Iter,DCofMx,VVMax
      end if

      if ((VVMax < 5.0e-16_wp) .and. (DCofMx < 1.0e-5_wp)) then
        Conv = .true.
        exit
      end if

      !call TStat(" Check T",1)
    end do

    if (.not. Change) then
      if (.not. Conv) write(IW,'(A,I3,A)') ' The number of iterations exceeded ',ItrEnd,' .'

      ! ===== Check the result =====

      call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
      if (StopBA) then
        Change = .true.
      else
        call FindAM(I_Dim+1,DD,DDMax,IMax)
        call ChkAcc(K_Lap,InitR,DDMax,R,Change)

        ! ===== Change parameters =====
        ! If the convergence isn't good, the larger R's parameters are used.

        write(IW,'(/A,I6)') ' Iteration =',Iter
      end if
    end if
  end if
end if
if (Change .or. SkpRem) then
  write(IW,'(1X,A)') 'Change!!'
  call SetExp(K_Lap,InitR,RIni,Coeff,T)
end if

! ===== Determine Xi =====

!call TStat(" Check T",1)
call SlvNt1(K_Lap,1,Coeff,T)
if (Change) R = RIni
!call TStat(" Newton1",1)

write(IW,'(/A,F10.4/)') ' Optimized solution in R = ',R
do I=1,K_Lap
  Idx = 2*I-1
  write(IW,'(A,I3,A,F21.18,2X,A,I3,A,F20.17)') ' Omega(',I,') = ',Coeff(Idx),'Alpha(',I,') = ',Coeff(Idx+1)
end do
write(IW,*)
call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
do I=1,I_Dim+1
  write(IW,'(A,I2,A,F20.17,2X,A,F20.14)') ' Error(',I,') = ',DD(I),' X = ',DD(2*K_Lap+1+I)
end do
write(IW,*)
do I=1,I_Dim
  write(IW,'(A,I2,A,F20.15)') ' Xi(',I,') = ',T(I)
end do
write(IW,'()')

! ===== Multiply the scaling parameter (EMinIv = 1/EMin) =====

Coeff(1:I_Dim) = EMinIV*Coeff(1:I_Dim)
write(IW,'(/A/)') ' Final solution '
do I=1,K_Lap
  Idx = 2*I-1
  write(IW,'(A,I3,A,F21.18,2X,A,I3,A,F20.17)') ' Omega(',I,') = ',Coeff(Idx),'Alpha(',I,') = ',Coeff(Idx+1)
end do

call Remez_ShutDownPrint(Verbose)

return

end subroutine Remez
