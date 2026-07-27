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

use Index_Functions, only: iTri, nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NO, NI, NA, NS
real(kind=wp), intent(in) :: FIFA(nTri_Elem(NO))
real(kind=wp), intent(out) :: FIT(NI,NA), FTI(NA,NI), FIA(NI,NS), FAI(NS,NI), FTA(NA,NS), FAT(NS,NA)
integer(kind=iwp) :: IA, IAI, IAT, IATOT, II, IT, ITI, ITTOT

do IT=1,NA
  ITTOT = NI+IT
  do II=1,NI
    ITI = iTri(ITTOT,II)
    FIT(II,IT) = FIFA(ITI)
    FTI(IT,II) = FIFA(ITI)
  end do
end do
do IA=1,NS
  IATOT = NI+NA+IA
  do II=1,NI
    IAI = iTri(IATOT,II)
    FIA(II,IA) = FIFA(IAI)
    FAI(IA,II) = FIFA(IAI)
  end do
end do
do IA=1,NS
  IATOT = NI+NA+IA
  do IT=1,NA
    ITTOT = NI+IT
    IAT = iTri(IATOT,ITTOT)
    FTA(IT,IA) = FIFA(IAT)
    FAT(IA,IT) = FIFA(IAT)
  end do
end do

end subroutine FBLOCK
