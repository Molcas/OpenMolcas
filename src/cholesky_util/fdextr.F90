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

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(ZERO=0.0D+00,ONE=1.0D+00,PT5=5.0D-01,MxIter=10000,Thrs=1.0D-09,TWO=2.0D+00)
real*8 Coeff(40), T(40), DD(82)
logical Dbg, StpBA
IDim = 2*K_Lap
IDimEnd = IDim+1
Dbg = .false.
StpBA = .false.
do I=1,IDimEnd

  ! ===== End points =====

  if (I == 1) then
    XXMax = ONE
    goto 888
  else if (I == IDimEnd) then
    XXMax = R
    goto 888
  end if

  ! ===== Initial values (midpoint)

  X1 = T(I-1)
  X2 = T(I)
  X = (X1+X2)*PT5

  ! ===== Solve Equations =====

  Theta = ONE
  do Iter=1,MxIter
    FX = GetDr1(K_Lap,X,Coeff)
    DFX = GetDr2(K_Lap,X,Coeff)
100 continue
    XNew = X-Theta*FX/DFX
    DifX = abs(XNew-X)
    if (Dbg) write(IW,*) Iter,XNew,Difx
    if (DifX < Thrs) goto 777
    FNew = GetDr1(K_Lap,XNew,Coeff)
    FFF = (ONE-Theta*PT5)*FX
    if (abs(FNew) >= abs(FFF)) then
      if (Dbg) write(IW,*) FNew,FFF
      Theta = Theta*PT5
      goto 100
    end if
    X = XNew
  end do

  write(IW,'(A)') '*************** Max Iteration in FdExtr'
  write(IW,'(A,I3,A,E23.15E3)') 'I =',I,' Max DifX. =',DifX
  !StpBA = .true.
  !return

  FMax = ZERO
  XXMax = ZERO

  IX = -1
  XErr2 = 9.99D+02
  XErr3 = XErr2
  IDr = 1000
  IDrEnd = IDr+1
  DrInv = ONE/dble(IDr)
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

  if (IX == 1) goto 888
  if (IX == IDrEnd) goto 888

  Dum = X1+(IX-2)*DrDif
  XM1 = QuadErr(K_Lap,Dum,Coeff)

  Dum = X1+(IX-1)*DrDif
  XM2 = QuadErr(K_Lap,Dum,Coeff)

  Dum = X1+IX*DrDif
  XM3 = QuadErr(K_Lap,Dum,Coeff)

  XErr1 = (XM3-XM1)/(TWO*(XM3-TWO*XM2+XM1))
  XErr2 = XErr1*DrInv-XXMax
  XMx = max(abs(XM1),abs(XM3))
  if (abs(XM2) <= XMx) goto 888
  XXMax = XErr2

777 continue
  XXMax = XNew

  ! ===== Save data =====

888 continue
  FMax = QuadErr(K_Lap,XXMax,Coeff)
  DD(I) = FMax
  DD(I+IDimEnd) = XXMax
end do

return

end subroutine FdExtr
