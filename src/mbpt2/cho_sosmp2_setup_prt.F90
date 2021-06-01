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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Setup_Prt(irc)
! Francesco Aquilante  May 2007
!
! Purpose: print setup for SOS-MP2.

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
#include "chomp2_cfg.fh"
#include "chomp2.fh"

irc = 0

call Cho_Head('Cholesky SOS-MP2 Setup','=',80,u6)
write(u6,*)

if (nBatch > 1) then
  write(u6,'(A,I6,A,I6,A)') 'The list of',nOccT,' occupied orbitals has been split in',nBatch,' batches:'
  write(u6,*) 'Batching is not allowed in SOS-MP2 : I stop here! '
  call Abend()
else if (nBatch == 1) then
  write(u6,'(A,I6,A)') 'The list of',nOccT,' occupied orbitals is not split:'
else
  write(u6,*) 'Oops, #batches over occupied orbitals is non-positive: ',nBatch
  irc = -101
  return
end if

write(u6,'(//,A)') 'The following tasks will be performed:'
write(u6,'(A)') ' * AO-to-MO transformation of original Cholesky vectors.'
if (DecoMP2) then
  write(u6,'(A)') ' * Cholesky decomposition of M=(ai|bj)^2 matrix.'
end if
write(u6,*) ' * Calculation of SOS-MP2 correlation energy.'

call xFlush(u6)

end subroutine Cho_SOSmp2_Setup_Prt
