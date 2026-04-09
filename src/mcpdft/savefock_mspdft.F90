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
! Copyright (C) 2018, Andrew M. Sand                                   *
!               2019, Thais R. Scott                                   *
!               2021, Jie J. Bao                                       *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine SaveFock_msPDFT(cmo,h1e,D1Act,NQ,p2d,state)

use PrintLevel, only: DEBUG
use mcpdft_output, only: iPrLoc
use wadr, only: fockocc
use rasscf_global, only: ISTORP, nacpar, nacpr2, nFint, ntot4
use mspdftgrad, only: F1MS, F2MS, FocMS, FxyMS
use general_data, only: nbas, nfro, norb, nsym, ntot1, ntot2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cmo(ntot2), h1e(ntot1), D1Act(nacpar), P2D(nacpr2)
integer(kind=iwp), intent(in) :: NQ, state
integer(kind=iwp) :: iprlev
real(kind=wp), allocatable :: dm2(:), FA_V(:), FI_V(:), fock(:), h1e_mo(:), ONTOPO(:), ONTOPT(:), Q(:)

call mma_allocate(fock,ntot4,label='fock')

fock(:) = Zero

write(u6,'(2X,A)') 'Calculating potentials for analytic gradients for MS-PDFT'

IPRLEV = IPRLOC(3)

call mma_allocate(h1e_mo,ntot1,label='h1e_mo')
call ao2mo_1e(cmo,h1e,h1e_mo,nsym,nbas,norb,nfro)

! loading one-electron potential and two-electron potential
! Used as F1 and F2 in equations 58 and 59 in Ref1.
call mma_allocate(ONTOPT,nfint,Label='OnTopT')
call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
call Get_dArray('ONTOPT',OnTopT,NFINT)
call Get_dArray('ONTOPO',OnTopO,NTOT1)

! Store for later...
call Get_TUVX(OnTopT,f2ms(:,state))

call mma_allocate(FI_V,Ntot1,Label='FI_V')
call mma_allocate(FA_V,Ntot1,Label='FA_V')
! Note that these are stored in MO basis
call Get_dArray('FI_V',FI_V,NTOT1)
call Get_dArray('FA_V',FA_V,NTOT1)

fi_v(:) = h1e_mo(:)+OnTopO(:)+FI_V(:)
F1MS(:,state) = fi_v(:)

! _______________________________________________________________________
! This next part is to generate the MC-PDFT generalized fock operator.
! The corrections (from the potentials) to FI and FA are built in the NQ
! part of the code, for efficiency's sake.  It still needs to be
! debugged.

! Reordering of the two-body density matrix.
if (ISTORP(NSYM+1) > 0) then
  call mma_allocate(dm2,ISTORP(NSYM+1),label='dm2')
  call PMAT_RASSCF(P2d,dm2)
else
  call mma_allocate(dm2,1,label='dm2')
  dm2(:) = Zero
end if

! Must add to existing fock operator (occ/act).
call mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
call fock_update(fock,fi_v,fa_v,D1Act,dm2,Q,NQ,OnTopT,CMO)
call mma_deallocate(Q)
call mma_deallocate(dm2)
call mma_deallocate(OnTopO)
call mma_deallocate(OnTopT)
call mma_deallocate(FI_V)
call mma_deallocate(FA_V)

focms(:,state) = fockocc(:)
FxyMS(:,state) = fock(:)
if (IPRLEV >= DEBUG) then
  write(u6,*) 'FOCC_OCC'
  call wrtmat(fockocc,1,ntot1,1,ntot1)
end if

call Put_iScalar('SA ready',1)

call mma_deallocate(fock)
call mma_deallocate(h1e_mo)

end subroutine SaveFock_msPDFT
