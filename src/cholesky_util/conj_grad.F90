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

subroutine Conj_Grad(Done,lVector,Prec,X,XTemp,R,RTemp,P,PTemp,Z,ZTemp,AP,Tolerance,res)
!***********************************************************************
!
! Purpose: Takes one step in a PCG algorithm. The preconditioner has
!          to be diagonal at the moment.
!
! Algorithm: (Ref: Kolla upp referensen)
!
! Initial values:         suggestions to make it work:
! x1 = x1                 x1 = 0
! r1 = b-A*x_1            r1 = b
! p1 = p1                 p1 = r1
! z1 = (A_pre)^-1 * r1    z1 = (diag(A))^-1 * r1
!
! CG-iterations
! alfa  = (r_k * z_k) / (p_k * [A * p_k])
! x_k+1 = x_k + alfa * p_k
! r_k+1 =  r_k - alfa * [A * p_k]
! beta  = (r_k+1 * z_k+1) / (r_k * z_k)
! p_k+1 =  z_k+1 + Beta * p_k
!
! Apart from a check for convergence, that is all of it.
!
!***********************************************************************
!
! Parameters:
!
! Done:      A logical parameter that decides if there is convergence
!
! X:         A vector with a value for x_k which will be overwritten
!            with x_k+1.
!
! Xtemp:     A temporary scratch vector which will be used to replace X.
!            If its defined outside the routine a lot of allocating and
!            deallocating is avoided.
!
! R:         A vector with a value for r_k which will be overwritten
!            with r_k+1.
!
! Xtemp:     A temporary scratch vector which will be used to replace R.
!
! P:         A vector with a value for p_k which will be overwritten
!            with p_k+1.
!
! Ptemp:     A temporary scratch vector which will be used to replace P.
!
! lVector:   The vectorial dimension of the problem.
!            A is assumed to be a square matrix.
!
! Prec:      A diagonal preconditioner.
!
! AP:        A vector containing the product of A and P.
!            The reason for not specifying A is that the dimension may
!            be too big to keep in the workspace.
!
! Res:       Returns the norm of the residual. Can be useful for the
!            user to make a judgement call if the procedure does not
!            converge below the tolerance.
!
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: Done
integer(kind=iwp) :: lVector
real(kind=wp) :: AP(lVector), P(lVector), Prec(lVector), PTemp(lVector), R(lVector), res, RTemp(lVector), Tolerance, X(lVector), &
                 XTemp(lVector), Z(lVector), ZTemp(lVector)
real(kind=wp) :: Alfa, Beta
real(kind=wp), external :: ddot_

! Copy the vectors into the temporary ones. ipXTemp, ipRTemp etc. will
! then contain x_k, r_k etc.
XTemp(:) = X(:)
RTemp(:) = R(:)
PTemp(:) = P(:)
ZTemp(:) = Z(:)

! Now we calculate Alfa = (r_k * z_k) / (p_k * [A * p_k])
Alfa = ddot_(lVector,RTemp,1,ZTemp,1)/ddot_(lVector,PTemp,1,AP,1)

! Now we calculate x_k+1 = x_k + alfa * p_k

X(:) = X(:)+Alfa*PTemp(:)

! Now we calculate r_k+1 = r_k - alfa * [A * p_k]
R(:) = R(:)-Alfa*AP(:)

! At this point we can both see if we are done and
Res = sqrt(ddot_(lVector,R,1,R,1))
if (Res < Tolerance) then
  Done = .true.
else

  ! Now we calculate z_k+1 =
  Z(:) = R(:)*Prec(:)

  ! Now we calculate beta = (r_k+1 * z_k+1) / (r_k * z_k)
  Beta = ddot_(lVector,R,1,Z,1)/ddot_(lVector,RTemp,1,ZTemp,1)

  ! Now we calculate P_k+1
  P(:) = Z(:)+Beta*PTemp(:)

end if

end subroutine Conj_Grad
