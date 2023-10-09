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
      REAL*8 FUNCTION gammln(xx)
      Implicit None
      REAL*8 xx

      INTEGER j
      REAL*8 ser,tmp,x,y
      Real*8, Parameter:: stp=2.5066282746310005d0
      Real*8, Parameter:: cof(6)=[76.18009172947146d0,
     &                           -86.50532032941677d0,
     &                            24.01409824083091d0,
     &                            -1.231739572450155d0,
     &                              .1208650973866179d-2,
     &                             -.5395239384953d-5]
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END FUNCTION gammln
