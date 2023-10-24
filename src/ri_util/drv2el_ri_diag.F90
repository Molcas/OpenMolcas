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
! Copyright (C) 1990,1991,1993,1998, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El_RI_Diag(ThrAO,TInt,nTInt)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!***********************************************************************

use SOAO_Info, only: iOffSO
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use RI_glob, only: nSkal_Valence
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ThrAO
integer(kind=iwp), intent(in) :: nTInt
real(kind=wp), intent(out) :: TInt(nTInt)
integer(kind=iwp) :: iIrrep, nAcc, nSkal
logical(kind=iwp) :: DoFock, DoGrad, Indexation

!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do not generate
! tables for indexation.

DoFock = .false.
DoGrad = .false.
Indexation = .false.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
nSkal_Valence = nSkal
!                                                                      *
!***********************************************************************
!                                                                      *
! Update iOffSO and call the Cholesky code which does this.

nAcc = 0
do iIrrep=0,nIrrep-1
  iOffSO(iIrrep) = nAcc
  nAcc = nAcc+nBas(iIrrep)
end do
call RI_XDiag(TInt,nTInt)
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drv2El_RI_Diag
