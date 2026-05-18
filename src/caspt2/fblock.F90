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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine FBLOCK(FIFA,NO,NI,NA,NS,FIT,FTI,FIA,FAI,FTA,FAT)
! Extract three rectangular submatrices FIT, FIA and FTA from the
! triangular matrix FIFA.
! SVC: add transposed submatrices to avoid complicated strides in the
! low-level sgm subroutines

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NO, NI, NA, NS
real(kind=wp), intent(in) :: FIFA((NO*(NO+1))/2)
real(kind=wp), intent(out) :: FIT(NI,NA), FTI(NA,NI), FIA(NI,NS), FAI(NS,NI), FTA(NA,NS), FAT(NS,NA)
integer(kind=iwp) :: IA, IAI, IAT, IATOT, II, IT, ITI, ITTOT

do IT=1,NA
  ITTOT = NI+IT
  do II=1,NI
    ITI = (ITTOT*(ITTOT-1))/2+II
    FIT(II,IT) = FIFA(ITI)
    FTI(IT,II) = FIFA(ITI)
  end do
end do
do IA=1,NS
  IATOT = NI+NA+IA
  do II=1,NI
    IAI = (IATOT*(IATOT-1))/2+II
    FIA(II,IA) = FIFA(IAI)
    FAI(IA,II) = FIFA(IAI)
  end do
end do
do IA=1,NS
  IATOT = NI+NA+IA
  do IT=1,NA
    ITTOT = NI+IT
    IAT = (IATOT*(IATOT-1))/2+ITTOT
    FTA(IT,IA) = FIFA(IAT)
    FAT(IA,IT) = FIFA(IAT)
  end do
end do

end subroutine FBLOCK
