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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine GetPDFTFocks(FMO1t,FMO2t,nTri)

use MCLR_Data, only: nAcPr2
use input_mclr, only: nRoots
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTri
real(kind=wp), intent(out) :: FMO1t(nTri,nRoots), FMO2t(nAcPr2,nRoots)

call Get_DArray('F1_PDFT',FMO1t,nTri*nRoots)
call Get_DArray('F2_PDFT',FMO2t,nAcPr2*nRoots)

end subroutine GetPDFTFocks
