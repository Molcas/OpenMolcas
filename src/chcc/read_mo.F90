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

subroutine read_mo(CMO,nfro,no,nv,ndel,nbas,nOrb)

use Data_Structures, only: DSBA_Type

implicit real*8(A-H,O-Z)
! declaration of calling arguments
type(DSBA_Type) CMO
integer lthCMO
integer nfro_scf(8)
integer nfro
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: CMO_t(:,:)
#include "SysDef.fh"

!... Read nSym, Energy, nBas, nOrb, nOcc, nFro, CMO and orbital energies from COMFILE

call Get_iArray('nFro',nFro_scf,1)
if (nFro_scf(1) /= 0) then
  write(6,*) 'Some orbitals were frozen in SCF!'
  call Abend()
end if

lthCMO = nBas*nBas
call mma_allocate(CMO_t,nBas,nBas,Label='CMO_t')
call Get_dArray_chk('Last orbitals',CMO_t,lthCMO)

! - transpose MO matrix, skip the frozen occupied orbitals

call mo_transp(CMO%A0,CMO_t(:,1+nfro:nOrb),no,nv,ndel,nbas)

call mma_deallocate(CMO_t)

return

end subroutine read_mo
