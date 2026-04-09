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
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************

subroutine Nr_Shells(nSkal)
!***********************************************************************
!                                                                      *
!     Object: to compute the number of unique shells in the input.     *
!                                                                      *
!     Author: Roland Lindh, Chemical Physics, University of Lund,      *
!             Sweden. January '98.                                     *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use BasisMode, only: All_Mode, Atomic, Auxiliary_Mode, Basis_Mode, Fragment_Mode, kCnttp, lCnttp, Valence_Mode, &
                     With_Auxiliary_Mode, With_Fragment_Mode
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: nSkal
integer(kind=iwp) :: iCnttp, nTest, iCnt, iAng, iShll, nExpi, nBasisi

!                                                                      *
!***********************************************************************
!                                                                      *
! Determine the number of shells

nSkal = 0
if ((Basis_Mode /= Valence_Mode) .and. (Basis_Mode /= Auxiliary_Mode) .and. (Basis_Mode /= Fragment_Mode) .and. &
    (Basis_Mode /= With_Auxiliary_Mode) .and. (Basis_Mode /= With_Fragment_Mode) .and. (Basis_Mode /= All_Mode)) then
  call WarningMessage(2,'Nr_Shells: illegal Basis_Mode')
  call Abend()
end if

if (.not. Atomic) then
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Molecular set up                                                   *
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *

  do iCnttp=1,nCnttp
    nTest = dbsc(iCnttp)%nVal-1
    do iCnt=1,dbsc(iCnttp)%nCntr

      do iAng=0,nTest
        iShll = dbsc(iCnttp)%iVal+iAng
        nExpi = Shells(iShll)%nExp
        if (nExpi == 0) cycle
        nBasisi = Shells(iShll)%nBasis
        if (nBasisi == 0) cycle

        if ((Basis_Mode == Valence_Mode) .and. (Shells(iShll)%Aux .or. Shells(iShll)%Frag)) cycle
        if ((Basis_Mode == Auxiliary_Mode) .and. (.not. Shells(iShll)%Aux)) cycle
        if ((Basis_Mode == Fragment_Mode) .and. (.not. Shells(iShll)%Frag)) cycle
        if ((Basis_Mode == With_Auxiliary_Mode) .and. Shells(iShll)%Frag) cycle
        if ((Basis_Mode == With_Fragment_Mode) .and. Shells(iShll)%Aux) cycle
        nSkal = nSkal+1
      end do  ! iAng
    end do    ! iCnt
  end do      ! iCnttp

else
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Atomic set up                                                      *
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *

  do iCnttp=kCnttp,lCnttp
    nTest = dbsc(iCnttp)%nVal-1
    do iAng=0,nTest
      iShll = dbsc(iCnttp)%iVal+iAng
      nExpi = Shells(iShll)%nExp
      if (nExpi == 0) cycle
      nBasisi = Shells(iShll)%nBasis
      if (nBasisi == 0) cycle

      if (Shells(iShll)%Frag) cycle
      nSkal = nSkal+1

    end do  ! iAng
  end do
  if (dbsc(kCnttp)%Aux) nSkal = nSkal+1 ! Add dummy shell
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Nr_Shells
