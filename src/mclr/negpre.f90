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
Module negpre
Integer, Parameter:: MXSTATE=10
Integer, Parameter:: luciv=31

Real*8 ERAS(MXSTATE)
Real*8 P1(MXSTATE*(MXSTATE+1)/2)
Real*8 P1INV(MXSTATE*(MXSTATE+1)/2)
Logical ngp
Real*8, Allocatable:: SS(:)
End Module negpre
