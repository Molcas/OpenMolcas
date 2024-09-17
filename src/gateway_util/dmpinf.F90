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
! Copyright (C) 1992, Roland Lindh                                     *
!               Markus P. Fuelscher                                    *
!***********************************************************************

subroutine DmpInf()
!***********************************************************************
!                                                                      *
! Object: to dump all input information on the file INFO.              *
!                                                                      *
! Called from: Seward                                                  *
!                                                                      *
! Calling    :                                                         *
!              Put_dArray                                              *
!              Get_dArray                                              *
!              Put_iArray                                              *
!              Get_iArray                                              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1992                                             *
!                                                                      *
! modified by M.P. Fuelscher                                           *
! - changed to used communication file                                 *
!***********************************************************************

use External_Centers, only: External_Centers_Dmp
use Basis_Info, only: Basis_Info_Dmp
use Center_Info, only: Center_Info_Dmp
use Symmetry_Info, only: Symmetry_Info_Dmp
use SOAO_Info, only: SOAO_Info_Dmp
use Sizes_of_Seward, only: Size_Dmp
use DKH_Info, only: DKH_Info_Dmp
use Gateway_Info, only: Gateway_Info_Dmp
use RICD_Info, only: RICD_Info_Dmp
use nq_Info, only: NQ_Info_Dmp
use rctfld_module, only: PCM_Info_Dmp
use Gateway_Info, only: Gateway_Info_Dmp

implicit none

!                                                                      *
!***********************************************************************
!                                                                      *
call SOAO_Info_Dmp()
call Basis_Info_Dmp()
call Center_Info_Dmp()
call Symmetry_Info_Dmp()
call Size_Dmp()
call DKH_Info_Dmp()
call Gateway_Info_Dmp()
call RICD_Info_Dmp()
!                                                                      *
!***********************************************************************
!                                                                      *
! Reaction field parameters

call PCM_Info_Dmp()
!                                                                      *
!***********************************************************************
!                                                                      *
call NQ_Info_Dmp()
!                                                                      *
!***********************************************************************
!                                                                      *
call External_Centers_Dmp()
!                                                                      *
!***********************************************************************
!                                                                      *
call DMP_EFP()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine DmpInf
