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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine CiOvlp(jRoot,S1,S2,CI_vec)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the overlap of the CI_vector, CI_vec, with the           *
!     a set of test vectors.                                           *
!                                                                      *
!     calling arguments:                                               *
!     jRoot   : integer                                                *
!               root identifier                                        *
!     S1      : array of real*8, input/output                          *
!               overlap matrix with test vectors                       *
!     S2      : array of real*8, input/output                          *
!               norm of the test configurations in the CI vector       *
!     CI_vec  : array of real*8, input                                 *
!               CI_vector                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
integer(kind=iwp), intent(in) :: jRoot
real(kind=wp), intent(inout) :: S1(lRoots,lRoots), S2(lRoots,lRoots)
real(kind=wp), intent(in) :: CI_vec(nConf)
integer(kind=iwp) :: iConf, kRef, kRoot
real(kind=wp) :: Sum1, Sum2

if (ITER == 1) return

do kRoot=1,nRoots
  Sum1 = Zero
  Sum2 = Zero
  do kRef=1,mxRef
    iConf = jCj(kRoot,kRef)
    if (iConf /= 0) then
      Sum1 = Sum1+cCI(kRoot,kRef)*CI_vec(iConf)
      Sum2 = Sum2+CI_vec(iConf)*CI_vec(iConf)
    end if
  end do
  S1(jRoot,kRoot) = abs(Sum1)
  S2(jRoot,kRoot) = Sum2
end do

return

end subroutine CiOvlp
