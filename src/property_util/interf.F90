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

subroutine Interf(i_root,Ene,isuseene,iscasvb)
!***********************************************************************
!                                                                      *
!     Object: Driver toward MOLDEN interface                           *
!                                                                      *
!***********************************************************************

use casvb_global, only: ifvb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: i_root, isuseene, iscasvb
real(kind=wp), intent(in) :: Ene(*)
integer(kind=iwp) :: iDum(7,8), iS, iUHF, LuTmp, nB, nB2
character(len=80) :: Note
character(len=10) :: Filename
real(kind=wp), allocatable :: AdCMOA(:), AdCMOB(:), CA(:,:), CB(:,:), EAB(:,:), OccA(:), OccB(:)
integer(kind=iwp), external :: isFreeUnit
#include "rasdim.fh"
#include "general.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute memory requirements and allocate memory

nB = 0
nB2 = 0
do iS=1,nSym
  nB = nB+nBas(iS)
  nB2 = nB2+nBas(iS)**2
end do

call mma_allocate(OccA,nB,label='OCCA')
call mma_allocate(OccB,nB,label='OCCB')
call mma_allocate(EAB,nB,2,label='ENERGY')
call mma_allocate(CA,nB,nB,label='CMOA')
call mma_allocate(CB,nB,nB,label='CMOB')
call mma_allocate(AdCMOA,nB2,label='AdCMOA')
call mma_allocate(AdCMOB,nB2,label='AdCMOB')
!                                                                      *
!***********************************************************************
!                                                                      *
! For the moment: Orbital energies just zero
if (isuseene /= 0) then
  EAB(:,1) = Ene(1:nB)
  EAB(:,2) = Ene(1:nB)
else
  EAB(:,:) = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the coeff. of sym adapted basis functions (CA, CB) and
! the spin orbital occupations (OccA, OccB)

call Dens_IF(i_root,CA,CB,OccA,OccB)
call Dens_IF_SCF(CA,AdCMOA,'B')
call Dens_IF_SCF(CB,AdCMOB,'B')
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out info on a temporary vector file.

Note = 'Temporary orbital file for the MOLDEN interface.'
LuTmp = 50
LuTmp = IsFreeUnit(LuTmp)
iUHF = IFVB
if (i_root /= 0) iUHF = 1
call WrVec_('TMPORB',LuTmp,'COE',iUHF,nSym,nBas,nBas,AdCMOA,AdCMOB,OccA,OccB,EAB(:,1),EAB(:,2),iDum,Note,0)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(OccA)
call mma_deallocate(OccB)
call mma_deallocate(EAB)
call mma_deallocate(CA)
call mma_deallocate(CB)
call mma_deallocate(AdCMOA)
call mma_deallocate(AdCMOB)
!                                                                      *
!***********************************************************************
!                                                                      *
if (i_root /= 0) then
  if (i_root <= 9) then
    write(filename,'(A7,I1)') 'MD_CAS.',i_root
  else if (i_root <= 99) then
    write(filename,'(A7,I2)') 'MD_CAS.',i_root
  else if (i_root <= 999) then
    write(filename,'(A7,I3)') 'MD_CAS.',i_root
  else
    filename = 'MD_CAS.x'
  end if
else
  filename = 'MD_CAS'
end if
if (iscasvb == 1) filename = 'MD_VB'
!                                                                      *
!***********************************************************************
!                                                                      *
! Call the generic MOLDEN interface

call Molden_Interface(iUHF,'TMPORB',filename)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Interf
