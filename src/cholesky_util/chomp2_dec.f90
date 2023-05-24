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
!
! Stuff for decomposing (ai|bj) integrals or amplitudes in MP2:
!
Module chomp2_dec
Public:: InCore, iOption_MP2CD, NowSym, ip_EOc, ip_EVir
Logical InCore(8)
Integer iOption_MP2CD, NowSym
Integer ip_EOc, ip_EVir
End Module chomp2_dec
