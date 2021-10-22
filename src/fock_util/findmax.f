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

       Real*8 Function FindMax(LaJ,numA)

       Implicit Real*8 (a-h,o-z)
       Integer numA
       Real*8 LaJ(NumA)

       Real*8 XMax
       Integer ia

          XMax = Abs(LaJ(1))

          Do ia=2,numA

             XMax = Max(XMax,Abs(LaJ(ia)))

          End Do

          FindMax = XMax

       Return
       End
