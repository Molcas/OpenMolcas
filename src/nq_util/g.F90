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
      Function G(Arg)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 G
      g=-1000.0D0
!
      Arg_=DBLE(Int(Arg))
      If (Abs(Arg-Arg_).lt.Half/Two) Then
!        Integer argument
         G=One
         rG=One
      Else
!        fractional argument
         G=Sqrt(Pi)
         rG=Half
      End If
!
 99   Continue
         If (Abs(rG-Arg).lt.Half/Two) goto 666
         G=rG*G
         rG=rG+One
      Go To 99
666   continue
      return
!
      End
