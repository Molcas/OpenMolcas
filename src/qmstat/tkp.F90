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
!
!-- Compute some auxiliary numbers.
!
      Subroutine TKP(Tau,dKappa,Rho,RhoA,RhoB,EA,EB,R                   &
     &              ,dNeigh,lTooSmall)
      Implicit Real*8 (a-h,o-z)


      Logical lTooSmall

      Tau=(EA-EB)/(EA+EB)
      Rho=0.5d0*(EA+EB)*R
      RhoA=(1+Tau)*Rho
      RhoB=(1-Tau)*Rho
      If(abs(Tau).gt.dNeigh) then
        dKappa=0.5d0*(Tau+1.0d0/Tau)
        lTooSmall=.false.
      Else
        lTooSmall=.true.
      Endif

      Return
      End
