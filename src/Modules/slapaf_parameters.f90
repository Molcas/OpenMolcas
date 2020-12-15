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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
Module Slapaf_Parameters
implicit none
Private
Public:: iRow, iRow_c, iInt, nFix, ddV_Schlegel, HWRS, iOptH, HUpMet, HrmFrq_Show, IRC, &
         nBVec, nDimBC, Curvilinear, Redundant, FindTS, User_Def, Analytic_Hessian
Integer:: iRow=0
Integer:: iRow_c=0
Integer:: iInt=0
Integer:: nFix=0
Integer:: IRC=0
Integer:: nBVec=0
Integer:: nDimBC=0

Logical:: Curvilinear=.True.
Logical:: Redundant=.False.
Logical:: FindTS=.False.
Logical:: HrmFrq_Show=.False.
Logical:: ddV_Schlegel=.False.
Logical:: HWRS=.True.
Logical:: User_Def=.False.
Logical:: Analytic_Hessian=.False.
!                                                                      *
!***********************************************************************
!                                                                      *
!     Hessian update
! 1   iOptH=00000001 (  1) Meyer (disabled)
! 2   iOptH=00000010 (  2) BP (disabled)
! 3   iOptH=00000100 (  4) BFGS
! 4   iOptH=00001000 (  8) None
! 5   iOptH=00010000 ( 16) MPS, for TS search
! 6   iOptH=-.1..... ( 32) Not used
! 7   iOptH=01000000 ( 64) EU, for TS search
! 8   iOptH=10000000 (128) TS-BFGS, for TS search
!
Integer:: iOptH=4
Character(LEN=6):: HUpMet=' None '

End Module Slapaf_Parameters
