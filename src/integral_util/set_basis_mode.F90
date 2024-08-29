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

subroutine Set_Basis_Mode(Label)

use BasisMode, only: All_Mode, Atomic, Auxiliary_Mode, Basis_Mode, Fragment_Mode, kCnttp, Valence_Mode, With_Auxiliary_Mode, &
                     With_Fragment_Mode

implicit none
character(len=*), intent(in) :: Label
character(len=7) :: Lbl

Atomic = .false.
kCnttp = 0
Lbl = Label(1:7)
call UpCase(Lbl)
select case (Lbl)
  case ('VALENCE')
    Basis_Mode = Valence_Mode
  case ('AUXILIA')
    Basis_Mode = Auxiliary_Mode
  case ('FRAGMEN')
    Basis_Mode = Fragment_Mode
  case ('WITHAUX')
    Basis_Mode = With_Auxiliary_Mode
  case ('WITHFRA')
    Basis_Mode = With_Fragment_Mode
  case ('ALL')
    Basis_Mode = All_Mode
  case default
    call WarningMessage(2,'Set_Basis_Mode: illegal mode, Label='//Lbl)
    call Abend()
end select

end subroutine Set_Basis_Mode
