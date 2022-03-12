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
!--------------------------------------------------------------------------*
! A small internal subroutine to compute the Nemo exchange repulsion.      *
!--------------------------------------------------------------------------*
      Real*8 Function ExNemo(i,j,a)
      Implicit Real*8 (a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
      Real*8 a
      Integer i,j
!The function
      ExNemo=Exp(-Sexrep(i,j)/a)*Sexre1(i,j)+Sexre2(i,j)*(a**20)
      Return
      End
