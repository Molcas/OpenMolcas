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

subroutine SlvNt2(K_Lap,R,Coeff,T,Theta2,VVMax,StopBA)
!-----------------------------------------------------------------------
! Function : Newton-Raphson method (type 2)
!
!     R      : The size of approximation range
!     T      : clossing points
!     TOld   : saved data of the first T
!     Coeff  : Omega & Alpha
!     CofOld : saved data of the first Coeff
!     DD     : maximum error & its coordinate
!     VVMax  : Maximum mean of VV array
!-----------------------------------------------------------------------

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(TOL=1.0D-22,ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,PT5=0.5D+00,TLim=2.0D-05)
logical Error, NG, Dbg, StopBA
real*8 Coeff(40), CofOld(40), T(40), TOld(40), VV(40), W(40), DD(82), A(40,40)

IDim = 2*K_Lap
New2 = 10
Theta2Mx = ONE
NG = .false.
Dbg = .false.
if (Dbg) write(IW,*) 'theta2',Theta2

call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
if (StopBA) return
call DCOPY_(IDim,Coeff,1,CofOld,1)

do I=1,IDim
  VV(I) = DD(I)+DD(I+1)
end do
Eps0 = FindMx(IDim,VV)
!-tbp: avoid the risk of using eps1 later:
Eps1 = Eps0

if (Eps0 > TOL) then
  Eps = 1.0D-03
  do J=1,IDim
    Eprel = T(J)*Eps
    EprInv = ONE/Eprel
    TCpy = T(J)
    T(J) = T(J)*(ONE+Eps)
    call SlvNt1(K_Lap,New2,Coeff,T)
    call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
    if (StopBA) return
    do I=1,IDim
      Temp = DD(I)+DD(I+1)-VV(I)
      A(I,J) = Temp*EprInv
    end do
    call DCOPY_(IDim,CofOld,1,Coeff,1)
    T(J) = TCpy
  end do
  call SlvEqs(IDim,A,W,VV,Error)
  if (.not. Error) goto 555
  call DCOPY_(IDim,T,1,TOld,1)

222 continue
  do I=1,IDim
    T(I) = TOld(I)-Theta2*W(I)
  end do
  call CkAltT(K_Lap,R,T,NG)
  if (NG) then
    write(IW,'(A)') '!! wrong T-values !!'
    call AbortG()
    New2 = 100
    call DCOPY_(IDim,TOld,1,T,1)
    goto 333
  end if
  call DCOPY_(IDim,CofOld,1,Coeff,1)
  call SlvNt1(K_Lap,New2,Coeff,T)
  call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
  if (StopBA) return
  do I=1,IDim
    VV(I) = DD(I)+DD(I+1)
  end do
  Eps1 = FindMx(IDim,VV)
  if (Dbg) write(IW,*) 'eps1',Eps1
  if (Eps1 < Eps0) then
    Theta2 = TWO*Theta2
    if (Theta2 > Theta2Mx) Theta2 = Theta2Mx
  end if
  goto 999

333 continue
  if (Theta2 < TLim) then
    write(IW,'(A)') ' Theta2 becomes too small.'
  else
    Theta2 = Theta2*PT5
    goto 222
  end if

555 continue
end if

999 continue
VVMax = Eps1
call SlvNt1(K_Lap,New2,Coeff,T)

if (Dbg) then
  write(IW,*) 'eps1',VVMax
  do I=1,IDim
    write(IW,*) 'coeff',i,coeff(i)
  end do
  do I=1,IDim
    write(IW,*) 'T',i,t(i)
  end do
end if

return

end subroutine SlvNt2
