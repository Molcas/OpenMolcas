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

subroutine coeff(ralpha,rbetaa,rbetas)

use MCLR_Data, only: mS2
use input_mclr, only: iSpin
use Constants, only: Zero, One, Two, Six, Eight, Half
use Definitions, only: wp, u6

implicit none
real(kind=wp), intent(out) :: ralpha, rbetaa, rbetas
real(kind=wp) :: rcg10, rcg11, rcg20, rcg21, rgamma, RMS, Spin
real(kind=wp), external :: clebsch_gordan_mclr

Spin = real(ispin-1,kind=wp)*Half
rms = real(ms2,kind=wp)*Half

if ((rms == Zero) .or. (rms /= spin)) then
  write(u6,*)
  write(u6,*) '====='
  write(u6,*)
  write(u6,*) 'Sorry, I am just able to calculate the'
  write(u6,*) 'Spin polariztion for high spin states'
  write(u6,*) 'Welcome back after you have recalculated'
  write(u6,*) 'your wave function'
  write(u6,*)
  write(u6,*)
  call Quit_OnUserError()
end if

rcg21 = clebsch_gordan_mclr(Two,One,spin,rms-one,spin,rms)
rcg11 = clebsch_gordan_mclr(One,One,spin,rms-one,spin,rms)
rcg20 = clebsch_gordan_mclr(Two,Zero,spin,rms,spin,rms)
rcg10 = clebsch_gordan_mclr(One,Zero,spin,rms,spin,rms)
rgamma = sqrt(spin*(spin+One)-rms*(rms-One))

ralpha = rMS**2
rBetaa = rMS/sqrt(Eight)*rgamma*rcg11/rcg10
rbetaS = Zero

if (abs(Two-spin) <= spin) then
  rAlpha = ralpha-rMS*rgamma*rcg21/(rcg20*sqrt(Six))
  rBetas = -rMS*rgamma*rcg21/(Two*sqrt(Six)*rcg20)
end if

end subroutine coeff
