!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Get_Can_Lorb(Ene,Fock,nO,nX,jOrb,Umat,iSym)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Ene(*), Fock(*), Umat(*)
integer(kind=iwp) :: nO, nX, jOrb(nO), iSym
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ii, ip_eta, ip_Z, ip_ZZ, j, nOx, nXx

if (nO < 1) return

call GetMem('eta_ik','Allo','Real',ip_eta,2*nX**2+1)
ip_Z = ip_eta+nX**2
ip_ZZ = ip_Z+nX
call FZero(Work(ip_eta),nX**2)
do i=1,nX
  ii = ip_eta+nX*(i-1)+i-1
  Work(ii) = Ene(i)
end do

nXx = max(1,nX)
nOx = max(1,nO)
call DGEMM_('N','N',nX,nO,nX,One,Work(ip_eta),nXx,Umat(1),nXx,Zero,Work(ip_Z),nXx)
call DGEMM_('T','N',nO,nO,nX,One,Umat(1),nXx,Work(ip_Z),nXx,Zero,Work(ip_eta),nOx)

call Eigen_Molcas(nO,Work(ip_eta),Work(ip_Z),Work(ip_ZZ))

call dcopy_(nO**2,Work(ip_eta),1,Umat,1)
do i=1,nO
  ii = ip_Z+i-1
  j = jOrb(i)
  Fock(j) = Work(ii)
end do
call GetMem('eta_ik','Free','Real',ip_eta,2*nX**2+1)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iSym)

end subroutine Get_Can_Lorb
