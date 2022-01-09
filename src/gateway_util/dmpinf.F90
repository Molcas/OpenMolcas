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

use iso_c_binding
use External_Centers
use Basis_Info, only: Basis_Info_Dmp
use Center_Info, only: Center_Info_Dmp
use Symmetry_Info, only: Symmetry_Info_Dmp
use SOAO_Info, only: SOAO_Info_Dmp
use Sizes_of_Seward, only: Size_Dmp
use DKH_Info, only: DKH_Info_Dmp
use Real_Info, only: Real_Info_Dmp
use RICD_Info, only: RICD_Info_Dmp
use Logical_Info, only: Logical_Info_Dmp
use nq_Info

implicit real*8(A-H,O-Z)
#include "stdalloc.fh"
#include "real.fh"
#include "rctfld.fh"

call DmpInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt,cQStrt,iQStrt,rQStrt)

! This is to allow type punning without an explicit interface
contains

subroutine DmpInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt,cQStrt,iQStrt,rQStrt)

  integer, target :: cRFStrt, iRFStrt, lRFStrt, cQStrt, iQStrt
  real*8, target :: rRFStrt, rQStrt
  integer, pointer :: p_cRF(:), p_iRF(:), p_lRF(:), p_cQ(:), p_iQ(:)
  real*8, pointer :: p_rRF(:), p_rQ(:)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call SOAO_Info_Dmp()
  call Basis_Info_Dmp()
  call Center_Info_Dmp()
  call Symmetry_Info_Dmp()
  call Size_Dmp()
  call DKH_Info_Dmp()
  call Real_Info_Dmp()
  call RICD_Info_Dmp()
  call Logical_Info_Dmp()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Reaction field parameters

  Len = ip_of_iWork(lRFEnd)-ip_of_iWork(lRFStrt)+1
  call c_f_pointer(c_loc(lRFStrt),p_lRF,[Len])
  call Put_iArray('RFlInfo',p_lRF,Len)

  Len = ip_of_Work(rRFEnd)-ip_of_Work(rRFStrt)+1
  call c_f_pointer(c_loc(rRFStrt),p_rRF,[Len])
  call Put_dArray('RFrInfo',p_rRF,Len)

  Len = ip_of_iWork(iRFEnd)-ip_of_iWork(iRFStrt)+1
  call c_f_pointer(c_loc(iRFStrt),p_iRF,[Len])
  call Put_iArray('RFiInfo',p_iRF,Len)

  Len = ip_of_iWork(cRFEnd)-ip_of_iWork(cRFStrt)+1
  call c_f_pointer(c_loc(cRFStrt),p_cRF,[Len])
  call Put_iArray('RFcInfo',p_cRF,Len)

  nullify(p_lRF,p_rRF,p_iRF,p_cRF)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Numerical integration information and parameters

  Len = ip_of_Work(rQEnd)-ip_of_Work(rQStrt)+1
  call c_f_pointer(c_loc(rQStrt),p_rQ,[Len])
  call Put_dArray('Quad_r',p_rQ,Len)

  Len = ip_of_iWork(iQEnd)-ip_of_iWork(iQStrt)+1
  call c_f_pointer(c_loc(iQStrt),p_iQ,[Len])
  call Put_iArray('Quad_i',p_iQ,Len)

  Len = ip_of_iWork(cQEnd)-ip_of_iWork(cQStrt)+1
  call c_f_pointer(c_loc(cQStrt),p_cQ,[Len])
  call Put_iArray('Quad_c',p_cQ,Len)

  nullify(p_rQ,p_iQ,p_cQ)
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
