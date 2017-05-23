************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*
      subroutine coeff(ralpha,rbetaa,rbetas)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "cstate_mclr.fh"
      Two=2.0d0
      rOne=1.0d0
      Zero=0.0d0
      Spin=DBLE(ispin-1)/Two
      rms=DBLE(ms2)/Two
      If (rms.eq.0.0d0.or.rms.ne.spin) then
      Write(6,*)
      Write(6,*) '====='
      Write(6,*)
      Write(6,*) 'Sorry, I am just able to calculate the'
      Write (6,*)'Spin polariztion for high spin states'
      Write(6,*) 'Welcome back after you have recalculated'
      write(6,*) 'your wave function'
      Write(6,*)
      Write(6,*)
      Call Quit_OnUserError()
      end if
      rcg21=clebsch_gordon(Two,rone,spin,rms-rone,spin,rms)
      rcg11=clebsch_gordon(rOne,rOne,spin,rms-rone,spin,rms)
      rcg20=clebsch_gordon(Two,Zero,spin,rms,spin,rms)
      rcg10=clebsch_gordon(rOne,Zero,spin,rms,spin,rms)
      rgamma=sqrt(spin*(spin+1.0D0)-rms*(rms-1.0D0))
*
      ralpha=rMS**2
      rBetaa=rMS/sqrt(8.0d0)*rgamma*rcg11/rcg10
      rbetaS=0.0d0
      If (abs(2.0d0-spin).le.spin) Then
      rAlpha=ralpha-rMS*rgamma*rcg21/(rcg20*sqrt(6.0d0))
      rBetas=-rMS*rgamma*rcg21/(sqrt(24.0d0)*rcg20)
      End If
*
      return
      end
