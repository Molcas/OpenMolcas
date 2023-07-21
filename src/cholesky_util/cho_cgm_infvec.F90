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
! Copyright (C) Thomas Bondo Pedersen                                  *
!               2020,2021, Roland Lindh                                *
!***********************************************************************

subroutine Cho_CGM_InfVec(InfVcT,NVT,n)

implicit none
integer, pointer :: InfVcT(:,:,:)
integer n
integer NVT(n)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Cho_X_GetIP_InfVec(InfVcT)
    integer, pointer :: InfVct(:,:,:)
  end subroutine Cho_X_GetIP_InfVec
end interface

call Cho_X_GetIP_InfVec(InfVcT)
call Cho_X_GetTotV(NVT,n)

end subroutine Cho_CGM_InfVec
