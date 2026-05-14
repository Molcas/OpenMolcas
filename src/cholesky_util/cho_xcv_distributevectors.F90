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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!               2012,2014, Victor P. Vysotskiy                         *
!***********************************************************************

subroutine Cho_XCV_DistributeVectors(irc,SP_BatchDim,nSP_Batch,idSP,n_idSP,NVT,l_NVT)
!
! Thomas Bondo Pedersen, April 2010.
!
! Parallel execution: distribute vectors across nodes.
! Serial execution: reorder vectors on tmp files and write them to
! permanent vector files.
!
! Victor P. Vysotskiy, 2012:
! Number of 'ga_put' has been remarkably reduced.
! Victor P. Vysotskiy, 2014:
! Number of 'ga_get' has been remarkably reduced by using the stripped mode

use Cholesky, only: Cho_Real_Par, INF_PASS, IPRINT, LuPri
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSP_Batch, SP_BatchDim(nSP_Batch), n_idSP, idSP(n_idSP), l_NVT, NVT(l_NVT)
real(kind=wp) :: C0, C1, W0, W1

irc = 0
if (Cho_Real_Par) then
  if (iPrint >= Inf_Pass) call CWTime(C0,W0)
  call Cho_XCV_DV_P(irc,SP_BatchDim,nSP_Batch,idSP,n_idSP,NVT,l_NVT)
  if (iPrint >= Inf_Pass) then
    call CWTime(C1,W1)
    write(LuPri,'(/,1X,A)') 'Timing of vector distribution:'
    call Cho_PrtTim(' ',C1,C0,W1,W0,-1)
  end if
else
  if (iPrint >= Inf_Pass) call CWTime(C0,W0)
  call Cho_XCV_DV_S(irc,SP_BatchDim,nSP_Batch,idSP,n_idSP)
  if (iPrint >= Inf_Pass) then
    call CWTime(C1,W1)
    write(LuPri,'(/,1X,A)') 'Timing of vector write:'
    call Cho_PrtTim(' ',C1,C0,W1,W0,-1)
  end if
end if

end subroutine Cho_XCV_DistributeVectors
