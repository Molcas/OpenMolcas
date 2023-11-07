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
      real*8 function  gammaf(x)
!***********************************************************************
!                                                                      *
! Object: to compute the angular contribution to the multipole integral*
!         between continuum basis functions within an R-matrix sphere  *
!         (phi integration)                                            *
!                                                                      *
!***********************************************************************
      use rmat, only: lcosf, lsinf, n_gam, m_gam
      use Constants, only: Zero, One, Two
      Implicit None
      Real*8 x

      Integer k1, k2
      Real*8, external:: dgamma_molcas
      Real*8 arg1, arg2, arg3
!
      lcosf=n_gam
      lsinf=m_gam
      k1=(-1)**lsinf
      k2=(-1)**lcosf
      if(k1.eq.(-1).or.k2.eq.(-1)) then
       gammaf=Zero
      else
       arg1=(DBLE(lcosf)+One)/Two
       arg2=(DBLE(lsinf)+One)/Two
       arg3=(DBLE(lsinf)+DBLE(lcosf)+Two)/Two
       gammaf=Two*dgamma_molcas(arg1)*dgamma_molcas(arg2)/
     >                                  dgamma_molcas(arg3)
      Endif
!
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_real(x)
      End
