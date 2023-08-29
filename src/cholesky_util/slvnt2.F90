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

use ReMez_mod, only: IW
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: R
real(kind=wp), intent(inout) :: Coeff(40), T(40), Theta2
real(kind=wp), intent(out) :: VVMax
logical(kind=iwp), intent(out) :: StopBA
integer(kind=iwp) :: I, I_Dim, J, New2
real(kind=wp) :: A(40,40), CofOld(40), DD(82), Eprel, EprINv, Eps, Eps0, Eps1, TCpy, Temp, Theta, Theta2Mx, TOld(40), VV(40), W(40)
logical(kind=iwp) :: Error, NG, Dbg
real(kind=wp), parameter :: TLim = 2.0e-5_wp, TOL = 1.0e-22_wp
real(kind=wp), external :: FindMx

I_Dim = 2*K_Lap
New2 = 10
Theta2Mx = One
NG = .false.
Dbg = .false.
if (Dbg) write(IW,*) 'theta2',Theta2

call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
if (StopBA) return
CofOld(1:I_Dim) = Coeff(1:I_Dim)

VV(1:I_Dim) = DD(1:I_Dim)+DD(2:I_Dim+1)
Eps0 = FindMx(I_Dim,VV)
!-tbp: avoid the risk of using eps1 later:
Eps1 = Eps0

if (Eps0 > TOL) then
  Eps = 1.0e-3_wp
  do J=1,I_Dim
    Eprel = T(J)*Eps
    EprInv = One/Eprel
    TCpy = T(J)
    T(J) = T(J)*(One+Eps)
    call SlvNt1(K_Lap,New2,Coeff,T)
    call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
    if (StopBA) return
    do I=1,I_Dim
      Temp = DD(I)+DD(I+1)-VV(I)
      A(I,J) = Temp*EprInv
    end do
    Coeff(1:I_Dim) = CofOld(1:I_Dim)
    T(J) = TCpy
  end do
  call SlvEqs(I_Dim,A,W,VV,Error)
  if (Error) then
    TOld(1:I_Dim) = T(1:I_Dim)

    do
      T(1:I_Dim) = TOld(1:I_Dim)-Theta2*W(1:I_Dim)
      call CkAltT(K_Lap,R,T,NG)
      if (NG) then
        write(IW,'(A)') '!! wrong T-values !!'
        call WarningMessage(2,'Remez aborting!')
        call Abend()
        New2 = 100
        T(1:I_Dim) = TOld(1:I_Dim)
      else
        Coeff(1:I_Dim) = CofOld(1:I_Dim)
        call SlvNt1(K_Lap,New2,Coeff,T)
        call FdExtr(K_Lap,T,Coeff,R,Theta,DD,StopBA)
        if (StopBA) return
        VV(1:I_Dim) = DD(1:I_Dim)+DD(2:I_Dim+1)
        Eps1 = FindMx(I_Dim,VV)
        if (Dbg) write(IW,*) 'eps1',Eps1
        if (Eps1 < Eps0) then
          Theta2 = Two*Theta2
          if (Theta2 > Theta2Mx) Theta2 = Theta2Mx
        end if
        exit
      end if

      if (Theta2 < TLim) then
        write(IW,'(A)') ' Theta2 becomes too small.'
        exit
      else
        Theta2 = Theta2*Half
      end if
    end do

  end if
end if

VVMax = Eps1
call SlvNt1(K_Lap,New2,Coeff,T)

if (Dbg) then
  write(IW,*) 'eps1',VVMax
  do I=1,I_Dim
    write(IW,*) 'coeff',i,coeff(i)
  end do
  do I=1,I_Dim
    write(IW,*) 'T',i,t(i)
  end do
end if

return

end subroutine SlvNt2
