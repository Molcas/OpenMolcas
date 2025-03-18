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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************
      Function Fact(R)
      Real*8 R,Fact
      integer n,i,j
      n=nint(R)
      i=1
      If (n.eq.0) Then
        Fact=1.0d0
        return
      end if
      Do j=1,n
      i=i*j
      end do
      Fact=DBLE(i)
      return
      end
