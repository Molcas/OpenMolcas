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
! Copyright (C) 1992,2020, Roland Lindh                                *
!***********************************************************************

subroutine GetInf(DoRys,nDiff)
!***********************************************************************
!                                                                      *
! Object: to read all input information on the file INFO.              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1992                                             *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Real_Spherical, only: lMax_Internal, Sphere
use Her_RW, only: nPrp
use External_Centers, only: nOrdEF
use Gateway_global, only: Test
use DKH_Info, only: DKroll
use Sizes_of_Seward, only: S
use rctfld_module, only: lMax, cRFStrt, iRFStrt, lRFStrt, rRFStrt, cRFEnd, iRFEnd, lRFEnd, rRFEnd

implicit none
integer nDiff
logical DoRys
#include "SysDef.fh"
integer Len
integer, external :: ip_of_work, ip_of_iWork

call GetInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt)

! This is to allow type punning without an explicit interface
contains

subroutine GetInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt)

  integer, target :: cRFStrt, iRFStrt, lRFStrt
  real*8, target :: rRFStrt
  integer, pointer :: p_cRF(:), p_iRF(:), p_lRF(:)
  real*8, pointer :: p_rRF(:)

  ! Load the dynamic input area.

  call Get_Info_Dynamic()

  ! Load the static input area.

  call Get_Info_Static()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Reaction field parameters

  Len = ip_of_iWork(lRFEnd)-ip_of_iWork(lRFStrt)+1
  call c_f_pointer(c_loc(lRFStrt),p_lRF,[Len])
  call Get_iArray('RFlInfo',p_lRF,Len)

  Len = ip_of_Work(rRFEnd)-ip_of_Work(rRFStrt)+1
  call c_f_pointer(c_loc(rRFStrt),p_rRF,[Len])
  call Get_dArray('RFrInfo',p_rRF,Len)

  Len = ip_of_iWork(iRFEnd)-ip_of_iWork(iRFStrt)+1
  call c_f_pointer(c_loc(iRFStrt),p_iRF,[Len])
  call Get_iArray('RFiInfo',p_iRF,Len)

  Len = ip_of_iWork(cRFEnd)-ip_of_iWork(cRFStrt)+1
  call c_f_pointer(c_loc(cRFStrt),p_cRF,[Len])
  call Get_iArray('RFcInfo',p_cRF,Len)

  nullify(p_lRF,p_rRF,p_iRF,p_cRF)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate the transformation matrices

  if (S%iAngMx-1 >= lMax) then
    call Sphere(S%iAngMx)
    lmax_internal = S%iAngMx
  else
    call Sphere(lMax)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Set the highest number of differentiations which will be
  ! applied to the basis functions. In this case 2 + 1 ( the
  ! kinetic energy operator and a differentiaion with respect
  ! to a nuclear coordinate.

  nPrp = max(lMax,3)

  ! Setup of tables for coefficients of the Rys roots and weights.

  if (S%iAngMx == 0) nDiff = 2
  if (DKroll .and. (nOrdEF > 0)) nDiff = nDiff+nOrdEF
  if (.not. Test) call Setup_RW(DoRys,nDiff)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Set up for contracted calculation

  call Flip_Flop(.false.)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Get_EFP()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  return

end subroutine GetInf_Internal

end subroutine GetInf
