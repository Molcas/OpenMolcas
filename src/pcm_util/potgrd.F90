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

subroutine PotGrd(Temp,nGrad)

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
external PCMGrd1, PCMMmg
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "wldata.fh"
#include "rctfld.fh"
character Method*8, Label*80
real*8 Temp(nGrad)
logical DiffOp
real*8, allocatable :: D_Var(:)
! Statement function
nElem(i) = (i+1)*(i+2)/2

! Prologue
iRout = 131
iPrint = nPrint(iRout)
call CWTime(TCpu1,TWall1)

! Allocate memory for density

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

! Get the method label
!write(6,*) ' Read Method label'
call Get_cArray('Relax Method',Method,8)

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!write(6,*) ' Read density matrix'
call mma_allocate(D_Var,nDens,Label='D_Var')
call Get_D1ao_Var(D_var,nDens)

if (iPrint >= 99) then
  write(6,*) 'variational 1st order density matrix'
  ii = 1
  do iIrrep=0,nIrrep-1
    write(6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',D_Var(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! nOrdOp: order/rank of the operator
! Work(ip1): lOper of each component of the operator

nOrdOp = 0
nComp = nElem(nOrdOp)
call GetMem('Coor','Allo','Real',ipC,3*nComp)
call GetMem('lOper','Allo','Inte',ip1,nComp)
call dcopy_(nComp*3,[Zero],0,Work(ipC),1)
iWork(ip1) = 1
DiffOp = .true.
call dZero(Temp,nGrad)
call OneEl_g_mck(PCMGrd1,PCMMmG,Temp,nGrad,DiffOp,Work(ipC),D_Var,nDens,iWork(ip1),nComp,nOrdOp,Label)
call PrGrad_mck(' TEST (PCM) contribution',Temp,nGrad,ChDisp,5)
!                                                                      *
!***********************************************************************
!                                                                      *
call GetMem('lOper','Free','Inte',ip1,nComp)
call GetMem('Coor','Free','Real',ipC,3*nComp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end

call mma_deallocate(D_Var)

call CWTime(TCpu2,TWall2)
call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)

return

end subroutine PotGrd
