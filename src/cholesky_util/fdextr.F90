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

subroutine FdExtr(K_Lap,T,Coeff,R,Theta,DD,StpBA)
!-----------------------------------------------------------------------
! Function : Find extrema in each interval
!-----------------------------------------------------------------------

use ReMez_mod, only: IW
use Constants, only: Zero, One, Two, Half
use Definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: T(40), Coeff(40), R
real(kind=wp), intent(out) :: Theta, DD(82)
logical(kind=iwp), intent(out) :: StpBA
integer(kind=iwp) :: I, IDimEnd, IDr, IDrEnd, Iter, IX, J
real(kind=wp) :: DFX, DifX, DrDif, DrInv, Dum, FF, FFF, FMax, FNew, FX, X, X1, X2, XErr1, XErr2, XErr3, XM1, XM2, XM3, XMx, XNew, &
                 XXMax
logical(kind=iwp) :: Conv, Dbg
integer(kind=iwp), parameter :: MxIter = 10000
real(kind=wp), parameter :: Thrs = 1.0e-9_wp
real(kind=wp), external :: GetDr1, GetDr2, QuadErr

IDimEnd = 2*K_Lap+1
Dbg = .false.
StpBA = .false.
do I=1,IDimEnd

  ! ===== End points =====

  if (I == 1) then
    XXMax = One
  else if (I == IDimEnd) then
    XXMax = R
  else

    ! ===== Initial values (midpoint)

    X1 = T(I-1)
    X2 = T(I)
    X = (X1+X2)*Half

    ! ===== Solve Equations =====

    Conv = .false.
    Theta = One
    outer: do Iter=1,MxIter
      FX = GetDr1(K_Lap,X,Coeff)
      DFX = GetDr2(K_Lap,X,Coeff)
      do
        XNew = X-Theta*FX/DFX
        DifX = abs(XNew-X)
        if (Dbg) write(IW,*) Iter,XNew,Difx
        if (DifX < Thrs) then
          Conv = .true.
          exit outer
        end if
        FNew = GetDr1(K_Lap,XNew,Coeff)
        FFF = (One-Theta*Half)*FX
        if (abs(FNew) < abs(FFF)) exit
        if (Dbg) write(IW,*) FNew,FFF
        Theta = Theta*Half
      end do
      X = XNew
    end do outer

    if (Conv) then
      XXMax = XNew
    else
      write(IW,'(A)') '*************** Max Iteration in FdExtr'
      write(IW,'(A,I3,A,ES23.15E3)') 'I =',I,' Max DifX. =',DifX
      !StpBA = .true.
      !return

      FMax = Zero
      XXMax = Zero

      IX = -1
      XErr2 = 9.99e2_wp
      XErr3 = XErr2
      IDr = 1000
      IDrEnd = IDr+1
      DrInv = One/real(IDr,kind=wp)
      DrDif = (X2-X1)*DrInv

      do J=1,IDrEnd
        X = X1+(J-1)*DrDif
        FF = QuadErr(K_Lap,X,Coeff)
        XErr1 = XErr2
        XErr2 = XErr3
        XErr3 = FF
        if (abs(FF) > abs(FMax)) then
          FMax = FF
          XXMax = X
          IX = J
        end if
      end do
      XXMax = -XXMax

      if ((IX /= 1) .and. (IX == IDrEnd)) then

        Dum = X1+(IX-2)*DrDif
        XM1 = QuadErr(K_Lap,Dum,Coeff)

        Dum = X1+(IX-1)*DrDif
        XM2 = QuadErr(K_Lap,Dum,Coeff)

        Dum = X1+IX*DrDif
        XM3 = QuadErr(K_Lap,Dum,Coeff)

        XErr1 = (XM3-XM1)/(Two*(XM3-Two*XM2+XM1))
        XErr2 = XErr1*DrInv-XXMax
        XMx = max(abs(XM1),abs(XM3))
        if (abs(XM2) > XMx) then
          XXMax = XErr2
          XXMax = XNew
        end if
      end if
    end if
  end if

  ! ===== Save data =====

  FMax = QuadErr(K_Lap,XXMax,Coeff)
  DD(I) = FMax
  DD(I+IDimEnd) = XXMax
end do

return

end subroutine FdExtr
