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
! Copyright (C) 1999, Anders Bernhardsson                              *
!               1999, Coen de Graaf                                    *
!***********************************************************************

subroutine Dens_IF_SCF(COEFF,CMO,Mode)
! A small stupid interface for creating the alpha and beta
! occupation numbers and corresponding molecular orbitals,
! used in MOLDEN.
!
! EAW 990118
!
! For SCF it gets even more 'stupid', just reading the MO coefficients and
! dumping it into a large matrix of dimension nTot x nTot
!   (Coen de Graaf, Oct. 99)

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: COEFF(*), CMO(*)
character, intent(in) :: Mode
integer(kind=iwp) :: i, ip1, ip2, iS, nBas(0:7), nIrrep, nTot, nTot2

! COEN WANTED IT AS A BLOCKED MATRIX, SO HERE THEY COME...

call Get_iScalar('nSym',nIrrep)
call Get_iArray('nBas',nBas,nIrrep)
nTot = 0
nTot2 = 0
do iS=0,nIrrep-1
  nTot = nTot+nBas(iS)
  nTot2 = nTot2+nBas(iS)**2
end do

ip1 = 1
ip2 = 1
if (Mode == 'F') COEFF(:nTot**2) = Zero
!call RecPrt('COEFF',' ',Coeff,nTot,nTot)
do iS=0,nIrrep-1
  if (Mode == 'B') CMO(ip1:ip1+nbas(is)**2-1) = Zero
  !call RecPrt('CMO',' ',CMO(ip1),nbas(is),nbas(is))
  do i=1,nbas(is)
    if (Mode == 'F') call dcopy_(nbas(is),CMO(ip1),1,COEFF(ip2),1)
    if (Mode == 'B') call dcopy_(nbas(is),COEFF(ip2),1,CMO(ip1),1)
    ip1 = ip1+nbas(is)
    ip2 = ip2+nTot
  end do
  !call RecPrt('CMO',' ',CMO(ip1-nbas(is)**2),nbas(is),nbas(is))
  ip2 = ip2+nbas(is)
end do
!call RecPrt('COEFF',' ',Coeff,nTot,nTot)

end subroutine Dens_IF_SCF
