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
!***********************************************************************

subroutine CloseP()
!***********************************************************************
!                                                                      *
! Object: to close the handling of the 2nd order density matrix.       *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use pso_stuff, only: B_PT2, Bin, Case_MP2, CMO, CMOPT2, D0, DS, DSVar, DVar, G1, G2, G_Toc, Gamma_On, iOffAO, lPSO, LuGam, &
                     LuGamma, LuGamma2, LuGamma_PT2, SO2cI, SSDM, WRK1, WRK2
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: DoCholesky
character(len=8) :: Method_chk

call Get_cArray('Relax Method',Method_chk,8)
if (Method_chk == 'CASPT2  ') then
  if (allocated(iOffAO)) then
    ! Conventional ERIs
    call mma_deallocate(iOffAO)
    close(LuGamma_PT2)
    call mma_deallocate(CMOPT2)
    if (allocated(SSDM)) call mma_deallocate(SSDM)
    call mma_deallocate(WRK1)
    call mma_deallocate(WRK2)
  else
    ! RI/Cholesky ERIs
    call mma_deallocate(B_PT2)
    close(LuGamma2)
  end if
end if

if (case_mp2) then
  call DecideOnCholesky(DoCholesky)
  if (.not. DoCholesky) call DaClos(LuGam)
end if
!*********** columbus interface ****************************************
if (Gamma_On) then
  call DaClos(LuGamma)
  call mma_deallocate(Bin)
  call mma_deallocate(G_Toc)
  call mma_deallocate(SO2cI)
end if

if (lPSO) then
  call mma_deallocate(G2)
  call mma_deallocate(G1)
end if

call mma_deallocate(CMO)
call mma_deallocate(DSVar)
call mma_deallocate(DS)
call mma_deallocate(DVar)
call mma_deallocate(D0)

end subroutine CloseP
