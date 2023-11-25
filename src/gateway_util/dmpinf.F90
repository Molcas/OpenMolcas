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

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
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
use rctfld_module, only: CRFEnd, CRFStrt, iRFEnd, iRFStrt, lRFEnd, lRFStrt, rRFEnd, rRFStrt
use Gateway_Info, only: Gateway_Info_Dmp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: Length
integer(kind=iwp), external :: ip_of_iWork, ip_of_Work

call DmpInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt)

! This is to allow type punning without an explicit interface
contains

subroutine DmpInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt)

  integer(kind=iwp), target, intent(inout) :: cRFStrt, iRFStrt, lRFStrt
  real(kind=wp), target, intent(inout) :: rRFStrt
  integer(kind=iwp), pointer :: p_cRF(:), p_iRF(:), p_lRF(:)
  real(kind=wp), pointer :: p_rRF(:)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call SOAO_Info_Dmp()
  call Basis_Info_Dmp()
  call Center_Info_Dmp()
  call Symmetry_Info_Dmp()
  call Size_Dmp()
  call DKH_Info_Dmp()
  call Gateway_Info_Dmp()
  call RICD_Info_Dmp()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Reaction field parameters

  Length = ip_of_iWork(lRFEnd)-ip_of_iWork(lRFStrt)+1
  call c_f_pointer(c_loc(lRFStrt),p_lRF,[Length])
  call Put_iArray('RFlInfo',p_lRF,Length)

  Length = ip_of_Work(rRFEnd)-ip_of_Work(rRFStrt)+1
  call c_f_pointer(c_loc(rRFStrt),p_rRF,[Length])
  call Put_dArray('RFrInfo',p_rRF,Length)

  Length = ip_of_iWork(iRFEnd)-ip_of_iWork(iRFStrt)+1
  call c_f_pointer(c_loc(iRFStrt),p_iRF,[Length])
  call Put_iArray('RFiInfo',p_iRF,Length)

  Length = ip_of_iWork(cRFEnd)-ip_of_iWork(cRFStrt)+1
  call c_f_pointer(c_loc(cRFStrt),p_cRF,[Length])
  call Put_iArray('RFcInfo',p_cRF,Length)

  nullify(p_lRF,p_rRF,p_iRF,p_cRF)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call NQ_Info_Dmp()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call External_Centers_Dmp()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DMP_EFP()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  return

end subroutine DmpInf_Internal

end subroutine DmpInf
