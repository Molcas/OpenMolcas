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

!----------------------------------------------------------------------*
! A small function to compute the Nemo exchange repulsion.             *
!----------------------------------------------------------------------*
real*8 function ExNemo(i,j,a)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
real*8 a
integer i, j

!The function
ExNemo = exp(-Sexrep(i,j)/a)*Sexre1(i,j)+Sexre2(i,j)*(a**20)

return

end function ExNemo
