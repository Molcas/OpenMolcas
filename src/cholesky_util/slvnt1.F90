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

subroutine SlvNt1(K_Lap,NTimes,Coeff,T)
!-----------------------------------------------------------------------
! Function : Newton-Raphson method (type 1)
!            obtain omega & alpha
!-----------------------------------------------------------------------

use ReMez_mod, only: IW
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K_Lap, NTimes
real(kind=wp), intent(inout) :: Coeff(40)
real(kind=wp), intent(in) :: T(40)
integer(kind=iwp) :: I, I_Dim, IRes, Iter
real(kind=wp) :: A(40,40), CoeffOld(40), Eps0, Eps1, Theta, VV(40), W(40)
logical(kind=iwp) :: Dbg, Error
real(kind=wp), parameter :: THTMIN = 1.0e-4_wp, TOL = 1.0e-22_wp

Theta = One
Dbg = .false.
I_Dim = 2*K_Lap

outer: do Iter=1,NTimes
  IRes = 1
  call AltErr(K_Lap,Coeff,T,VV,Eps0)
  if (Eps0 > Tol) then
    call PtDiff(I_Dim,Coeff,T,A)

    if (Dbg) then
      write(IW,'(A)') 'Check values of derivatives'
      call Laplace_PRSQ(A,I_Dim,I_Dim,I_Dim)
    end if

    call SlvEqs(I_Dim,A,W,VV,Error)
    if (Dbg) then
      write(IW,*) 'Check solutions of equations',Error
      do I=1,I_Dim
        write(IW,*) I,W(I)
      end do
    end if
    if (Error) then
      CoeffOld(1:I_Dim) = Coeff(1:I_Dim)

      do
        Coeff(1:I_Dim) = CoeffOld(1:I_Dim)-Theta*W(1:I_Dim)
        call AltErr(K_Lap,Coeff,T,VV,Eps1)
        if (Eps1 < Eps0) then
          Theta = Two*Theta
          if (Theta > One) Theta = One
          cycle outer
        end if
        if (Theta < THTMIN) then
          cycle outer
        else
          Theta = Theta*Half
        end if
      end do
    end if
  end if
  IRes = 0
  if (IRes == 0) exit
end do outer

return

end subroutine SlvNt1
