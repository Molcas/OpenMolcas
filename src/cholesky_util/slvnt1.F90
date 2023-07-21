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

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(TOL=1.0D-22,ONE=1.0D+00,TWO=2.0D+00,THTMIN=1.0D-04,PT5=0.5D+00)
real*8 Coeff(40), CoeffOld(40), W(40), VV(40), T(40), A(40,40)
logical Error, Dbg

Theta = 1.0D+00
Dbg = .false.
IDim = 2*K_Lap

do Iter=1,NTimes
  IRes = 1
  call AltErr(K_Lap,Coeff,T,VV,Eps0)
  if (Eps0 > Tol) then
    call PtDiff(IDim,Coeff,T,A)

    if (Dbg) then
      write(IW,'(A)') 'Check values of derivatives'
      call Laplace_PRSQ(A,IDim,IDim,IDim)
    end if

    call SlvEqs(IDim,A,W,VV,Error)
    if (Dbg) then
      write(IW,*) 'Check solutions of equations',Error
      do I=1,IDim
        write(IW,*) I,W(I)
      end do
    end if
    if (.not. Error) goto 555
    call DCOPY_(IDim,Coeff,1,CoeffOld,1)

222 continue
    do I=1,IDim
      Coeff(I) = CoeffOld(I)-Theta*W(I)
    end do
    call AltErr(K_Lap,Coeff,T,VV,Eps1)
    if (Eps1 < Eps0) then
      Theta = TWO*Theta
      if (Theta > ONE) Theta = ONE
      goto 888
    end if
    if (Theta < THTMIN) then
      goto 888
    else
      Theta = Theta*PT5
      goto 222
    end if
555 continue
  end if
  IRes = 0
888 continue
  if (IRes == 0) goto 999
end do
999 continue

return

end subroutine SlvNt1
