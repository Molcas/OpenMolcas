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
! Copyright (C) 1997, Luis Serrano-Andres                              *
!***********************************************************************

subroutine SUPSCH(SMAT,CMOO,CMON)
! Program RASSCF
!
! Objective: To check the order of the input of natural orbitals
!            to obtain the right labels for the Supersymmetry matrix.
!
! Called from ortho, neworb, fckpt2, and natorb.
!
! Luis Serrano-Andres
! University of Lund, Sweden, 1997
! **** Molcas-4 *** Release 97 04 01 **********

use general_data, only: NSYM, NBAS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp) SMAT(*), CMOO(*), CMON(*)
integer(kind=iwp) :: iSym, nOrb_Tot, nOrbMx
integer(kind=iwp), allocatable :: IxSym2(:)
real(kind=wp), allocatable :: Temp1(:), Temp2(:)

nOrbMX = 0
nOrb_tot = 0
do iSym=1,nSym
  nOrbMX = max(nOrbMX,nBas(iSym))
  nOrb_tot = nOrb_tot+nBas(iSym)
end do

call mma_allocate(Temp1,nOrbMX*nOrbMX,Label='Temp1')
call mma_allocate(Temp2,nOrbMX*nOrbMX,Label='Temp2')
call mma_allocate(IxSym2,nOrb_tot,Label='IxSym2')

call SUPSCH_INNER(SMAT,CMOO,CMON,Temp1,Temp2,nOrbMX,IxSym2,nOrb_tot)

call mma_deallocate(IxSym2)
call mma_deallocate(Temp2)
call mma_deallocate(Temp1)

end subroutine SUPSCH
