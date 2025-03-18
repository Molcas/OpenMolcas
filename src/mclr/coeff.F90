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

use Constants, only: Zero, One, Two, Four, Six, Eight
use MCLR_Data, only: mS2
use input_mclr, only: iSpin

implicit none
real*8 ralpha, rbetaa, rbetas
real*8, parameter :: twentyfour = Four*Six
real*8 Spin, RMS, rcg21, rcg11, rcg20, rcg10, rgamma
real*8, external :: clebsch_gordan_mclr

Spin = dble(ispin-1)/Two
rms = dble(ms2)/Two

if ((rms == 0.0d0) .or. (rms /= spin)) then
  write(6,*)
  write(6,*) '====='
  write(6,*)
  write(6,*) 'Sorry, I am just able to calculate the'
  write(6,*) 'Spin polariztion for high spin states'
  write(6,*) 'Welcome back after you have recalculated'
  write(6,*) 'your wave function'
  write(6,*)
  write(6,*)
  call Quit_OnUserError()
end if

rcg21 = clebsch_gordan_mclr(Two,one,spin,rms-one,spin,rms)
rcg11 = clebsch_gordan_mclr(One,One,spin,rms-one,spin,rms)
rcg20 = clebsch_gordan_mclr(Two,Zero,spin,rms,spin,rms)
rcg10 = clebsch_gordan_mclr(One,Zero,spin,rms,spin,rms)
rgamma = sqrt(spin*(spin+One)-rms*(rms-One))

ralpha = rMS**2
rBetaa = rMS/sqrt(Eight)*rgamma*rcg11/rcg10
rbetaS = Zero

if (abs(Two-spin) <= spin) then
  rAlpha = ralpha-rMS*rgamma*rcg21/(rcg20*sqrt(Six))
  rBetas = -rMS*rgamma*rcg21/(sqrt(Twentyfour)*rcg20)
end if

end subroutine coeff
