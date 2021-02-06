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

module Vibrot_globals

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: nRot_Max = 200, nobs = 10, npin = 500, npoint = 5000
integer(kind=iwp) :: J1A, J2A, lambda, ispc, iobs, nop, Vibwvs, iadvib, Vibwvs1, Vibwvs2, n0, nvib1, n02, nvib21, J1B, J2B, &
                     IfPrWf, iscale, iallrot
integer(kind=iwp) :: iadrsp(nRot_Max), iad12(nRot_Max), iad13(nRot_Max), iplot(nobs), npobs(nobs)
real(kind=wp) :: R0o(nobs), R1o(nobs), dRo(nobs), RinO(npin,nobs), Obsin(npin,nobs), EoutO(npoint+4,nobs)
character(len=2) :: Atom1, Atom2
character(len=80) :: Titobs(nobs)

public :: Atom1, Atom2, dRo, EoutO, iad12, iad13, iadrsp, iadvib, iallrot, ifPrWf, iobs, iplot, iscale, ispc, J1A, J1B, J2A, J2B, &
          lambda, n0, n02, nop, npin, npobs, npoint, nRot_Max, nvib1, nvib21, Obsin, R0o, R1o, RinO, Titobs, Vibwvs, Vibwvs1, &
          Vibwvs2

end module Vibrot_globals
