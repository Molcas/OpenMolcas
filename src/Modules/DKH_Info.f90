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
Module DKH_Info
Private
Public :: nCtrLD, iCtrLD, radiLD, DKroll, LDKroll, BSS
Integer :: i
Integer ::  nCtrLD=0, iCtrLD(10)=[(0,i=1,10)]
Real*8  :: radiLD=0.0D0
Logical :: DKroll=.False.
Logical :: LDKroll=.False.
Logical :: BSS   =.False.
End Module DKH_Info
