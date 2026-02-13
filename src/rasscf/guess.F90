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
! Copyright (C) 1998, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Guess(CMO)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Diagonalize core Hamiltonian to get starting orbitals.           *
!                                                                      *
!     calling arguments:                                               *
!     CMO     : real*8, output                                         *
!               starting vectors                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1998                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoNuc, sNoOri
use general_data, only: NSYM, NBAS, NTOT1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
real*8 CMO(*)
character(len=8) Label
real*8, allocatable :: Tmp1(:)
integer iRC, i1, i2, iBas, iComp, iOpt, iSyLbl, iSym
#include "warnings.h"

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! allocate work space
call mma_allocate(Tmp1,nTot1,Label='Tmp1')

! load bare nuclei Hamiltonian

iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'OneHam  '
call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
if (iRc /= 0) then
  write(u6,*) ' RASSCF tried to construct start orbitals from'
  write(u6,*) ' diagonalization of core Hamiltonian, but ran into'
  write(u6,*) ' a severe error: Failed to read the Hamiltonian'
  write(u6,*) ' from the ONEINT file. Something may be wrong with'
  write(u6,*) ' the file.'
  call Quit(_RC_IO_ERROR_READ_)
end if

! diagonalize bare nuclei Hamiltonian

i1 = 1
i2 = 1
do iSym=1,nSym
  iBas = nBas(iSym)
  call dCopy_(iBas*iBas,[zero],0,CMO(i2),1)
  call dCopy_(iBas,[one],0,CMO(i2),iBas+1)
  call Jacob(Tmp1(i1),CMO(i2),iBas,iBas)
  call JacOrd(Tmp1(i1),CMO(i2),iBas,iBas)
  i1 = i1+(iBas*iBas+iBas)/2
  i2 = i2+iBas*iBas
end do

! deallocate work space
call mma_deallocate(Tmp1)

end subroutine Guess
