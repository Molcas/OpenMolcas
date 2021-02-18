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
Module KSDFT_Info
Character(LEN=16) KSDFA
Real*8 :: funcaa=0.0D0, funcbb=0.0D0 ,funccc=0.0D0
Real*8, Allocatable:: F_xca(:), F_xcb(:), tmpB(:)
Integer :: LuMC,LuMT
End Module KSDFT_Info
