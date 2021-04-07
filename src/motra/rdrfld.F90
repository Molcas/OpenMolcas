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

subroutine RdRfld(ipHOne)

implicit real*8(A-H,O-Z)
#include "files_motra.fh"
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
logical Found

!----------------------------------------------------------------------*
! If this is a perturbative reaction field calculation then            *
! modifiy the one-electron Hamiltonian by the reaction field and       *
! the nuclear attraction by the cavity self-energy                     *
!----------------------------------------------------------------------*
nTemp = 0
do iSym=1,nSym
  nTemp = nTemp+nBas(iSym)*(nBas(iSym)+1)/2
end do
call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
call f_Inquire('RUNOLD',Found)
if (Found) call NameRun('RUNOLD')
call get_dscalar('RF Self Energy',ERFself)
PotNuc = PotNuc+ERFself
call get_darray('Reaction field',Work(lTemp),nTemp)
if (Found) call NameRun('RUNFILE')
call Daxpy_(nTemp,1.0d0,Work(lTemp),1,Work(ipHone),1)
call GetMem('RFFLD','Free','Real',lTemp,nTemp)
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*
return

end subroutine RdRfld
