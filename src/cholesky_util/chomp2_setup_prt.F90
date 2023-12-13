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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_Setup_Prt(irc)
!
! Thomas Bondo Pedersen, Nov. 2004 / Feb. 2005.
!
! Purpose: print setup for Cholesky MP2.

use Cholesky, only: nSym
use ChoMP2, only: ChoAlg, DecoMP2, ForceBatch, iFirst, Laplace, Laplace_nGridPoints, LnBatOrb, LnOcc, nBatch, nBatOrbT, nDel, &
                  nFro, nOcc, NumBatOrb, NumOcc, nVir, SOS_MP2
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: iBatch, iCount(8), iSym

irc = 0

iCount(:) = 0

call Cho_Head('Cholesky MP2 Setup','=',80,u6)
! The values but not the names 'occupied' are updated to work
! also for batching over all orbitals
if (nBatch > 1) then
  write(u6,'(A,I6,A,I6,A)') 'The list of',nBatOrbT,' occupied orbitals has been split in',nBatch,' batches:'
else if (nBatch == 1) then
  write(u6,'(A,I6,A)') 'The list of',nBatOrbT,' occupied orbitals is not split:'
else
  write(u6,*) 'Oops, #batches over occupied orbitals is non-positive: ',nBatch
  irc = -101
  return
end if

write(u6,'(/,3X,A)') ' Batch  First   Last #Occ/irrep'
if (nSym == 1) then
  write(u6,'(3X,A)') '-------------------------------'
else if (nSym == 2) then
  write(u6,'(3X,A)') '-----------------------------------'
else if (nSym == 4) then
  write(u6,'(3X,A)') '-------------------------------------------------'
else if (nSym == 8) then
  write(u6,'(3X,A)') '-----------------------------------------------------------------------------'
else
  write(u6,*) 'Oops, #irreps is out of bounds: ',nSym
  irc = -102
  return
end if
do iBatch=1,nBatch
  if (.false.) then
    write(u6,'(3X,I6,1X,I6,1X,I6,8(1X,I6))') iBatch,iFirst(iBatch),iFirst(iBatch)+NumBatOrb(iBatch)-1, &
                                             (LnBatOrb(iSym,iBatch),iSym=1,nSym)
  else
    write(u6,'(3X,I6,1X,I6,1X,I6,8(1X,I6))') iBatch,iFirst(iBatch),iFirst(iBatch)+NumOcc(iBatch)-1,(LnOcc(iSym,iBatch),iSym=1,nSym)
  end if
  do iSym=1,nSym
    if (.false.) then
      iCount(iSym) = iCount(iSym)+LnBatOrb(iSym,iBatch)
    else
      iCount(iSym) = iCount(iSym)+LnOcc(iSym,iBatch)
    end if
  end do
end do
if (nSym == 1) then
  write(u6,'(3X,A)') '-------------------------------'
else if (nSym == 2) then
  write(u6,'(3X,A)') '-----------------------------------'
else if (nSym == 4) then
  write(u6,'(3X,A)') '-------------------------------------------------'
else if (nSym == 8) then
  write(u6,'(3X,A)') '-----------------------------------------------------------------------------'
end if
write(u6,'(3X,A,14X,8(1X,I6))') 'Total:',(iCount(iSym),iSym=1,nSym)
if (nSym == 1) then
  write(u6,'(3X,A)') '-------------------------------'
else if (nSym == 2) then
  write(u6,'(3X,A)') '-----------------------------------'
else if (nSym == 4) then
  write(u6,'(3X,A)') '-------------------------------------------------'
else if (nSym == 8) then
  write(u6,'(3X,A)') '-----------------------------------------------------------------------------'
end if
do iSym=1,nSym
  if (.false.) then
    if (iCount(iSym) /= nOcc(iSym)+nVir(iSym)+nFro(iSym)+nDel(iSym)) then
      write(u6,*) 'Oops, #Occ/irrep is incorrect....'
      irc = -103
      return
    end if
  else
    if (iCount(iSym) /= nOcc(iSym)) then
      write(u6,*) 'Oops, #Occ/irrep is incorrect....'
      irc = -103
      return
    end if
  end if
end do
if ((nBatch > 1) .and. ForceBatch) write(u6,'(/,A)') 'Notice: batching has been requested by user.'

write(u6,'(//,A)') 'The following tasks will be performed:'
write(u6,'(A)') ' * AO-to-MO transformation of original Cholesky vectors.'
if (DecoMP2) write(u6,'(A)') ' * Cholesky decomposition of (ai|bj) integrals.'
if (nBatch > 1) then
  if (DecoMP2) then
    write(u6,'(A)') ' * Presort of Cholesky vectors from (ai|bj) decomposition.'
  else
    write(u6,'(A)') ' * Presort of MO Cholesky vectors.'
  end if
end if
if (Laplace .and. SOS_MP2) then
  write(u6,'(A)') ' * Calculation of Laplace-SOS-MP2 correlation energy.'
  if (Laplace_nGridPoints == 0) then
    write(u6,'(A)') '   Numerical Laplace integration quadrature: default'
  else
    write(u6,'(A,I8)') '   Numerical Laplace integration quadrature:',Laplace_nGridPoints
  end if
else
  write(u6,'(A)') ' * On-the-fly assembly of (ai|bj) integrals and calculation of MP2 energy correction.'
  write(u6,'(A,I3,A)') '   [Cholesky algorithm:',ChoAlg,']'
end if

call xFlush(u6)

end subroutine ChoMP2_Setup_Prt
