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
      SubRoutine DMInvKap_td(DigPrec,rIn,rout)
!
!     Prec  Diagonal Preconditioner from Prec_td
!     rIn   nDensC long orb RHS
!     rout  Trail vector rout = T^-1B^x
!
      use MCLR_Data, only: nDensC
      Implicit None
      Real*8 rOut(*),rin(*),DigPrec(*)

      Integer k
!
!-------------------------------------------------------------------------
! Multiply the 1/precond in vector form, rTemp, with RHS in vector form, rIn
! ---> New trail vector.
!-------------------------------------------------------------------------
!
!
      Do k=1, nDensC
         Rout(k) = rIn(k)/DigPrec(k)
      End Do
!
!
      End SubRoutine DMInvKap_td
