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
! Copyright (C) 2026, Andrew M. Sand                                   *
!***********************************************************************

subroutine molel_write(CMO)

use onedat, only: snoori
use rasscf_global, only: lRoots
use general_data, only: nactel, nfro, nish, nash, nbas, ntot, ntot1, ntot2, nsym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u0
use Constants, only: One, Zero

implicit none
#include "warnings.h"
integer(kind=iwp) :: comp, opt, rc, sylbl, counter, i, ibas, tot_el
integer(kind=iwp) :: isym, j, nactorb, nAOs, offset, molet_lu
integer(kind=iwp), external :: isFreeUnit
character(len=8) :: label
character(len=80) :: FMT1
real(kind=wp), intent(in) :: cmo(ntot2)
real(kind=wp), allocatable :: int1e_ovlp(:), OvMat(:), F1_PDFT(:), tmp1(:), &
tmp3(:), tmp4(:), tmp5(:)


!overlap matrix - AO
rc = -1
opt = ibset(0,sNoOri)
comp = 1
sylbl = 1
label = 'Mltpl  0'

molet_lu = IsFreeUnit(40)

call mma_allocate(int1e_ovlp,ntot1+4,label='int1e_ovlp')
call mma_allocate(OvMat,ntot**2,label='OvMat')
call rdone(rc,opt,label,comp,int1e_ovlp,sylbl)
call molcas_open(molet_lu,'MolEl.dat')
nAOs = 0
do isym=1,nsym
  nAOs = nAOs + nBas(iSym)
end do

Tot_El = 2*sum(nfro+nish)+nactel


FMT1 = "(I4,1X,I4,1X,I4,1X,I4,1X,I4)"
write(molet_lu,*) 'Number of states,orbitals,electrons,ActOrb,ActEl'
write(molet_lu,FMT1) lroots,nAOs,Tot_El,nActOrb,nActEl
write(molet_lu,*) 'Overlap Matrix (AO)'

FMT1 = "(I4,1X,I4,1X,F20.14)"
counter = 0
offset = 0
do isym=1,nsym
   iBas = nBas(iSym)
   do i=1,iBas
      do j=1,i
         write(molet_lu,FMT1) i+offset,j+offset,int1e_ovlp(1+counter)
         OvMat(1+(i+offset-1)*ntot+(j+offset-1)) = int1e_ovlp(1+counter)
         OvMat(1+(j+offset-1)*ntot+(i+offset-1)) = int1e_ovlp(1+counter)
         counter = counter + 1
      end do
   end do
   offset = offset + iBas
end do
call mma_deallocate(int1e_ovlp)



write(molet_lu,*) 'Molecular orbital coefficients'
counter = 0
offset = 0
do isym=1,nsym
   iBas = nBas(iSym)
   do i=1,iBas
      do j=1,iBas
         write(molet_lu,FMT1) i+offset,j+offset,CMO(1+counter)
         counter = counter +1
      end do
   end do
   offset = offset + iBas
end do

!Read this in, but it is in MO - convert to AO
!call mma_allocate(F1_PDFT,NTOT1,Label='F1p')
!call get_darray('F1_PDFT',FMO1t,ntot1)
      write(molet_lu,*) 'Effective Hamiltonian(s) (AO)'

!S * C * F_mo * C^T * S
call mma_allocate(F1_PDFT,NTOT1,Label='F1p')
call mma_allocate(tmp1,ntot**2,Label='tmp1')
call mma_allocate(tmp3,ntot**2,Label='tmp3')
call mma_allocate(tmp4,ntot**2,Label='tmp4')
call mma_allocate(tmp5,ntot**2,Label='tmp5')


call DGEMM_('N','N',ntot,ntot,ntot, &
              One,OvMat,ntot, &
              CMO,ntot, &
              Zero,tmp1,ntot)

call get_darray('F1_PDFT',F1_PDFT,ntot1)

call Square(F1_PDFT,tmp3,1,ntot,ntot)

call DGEMM_('N','N',ntot,ntot,ntot, &
              One,Tmp1,ntot, &
              Tmp3,ntot, &
              Zero,Tmp4,ntot)
call DGEMM_('N','T',ntot,ntot,ntot, &
              One,Tmp4,ntot, &
              Tmp1,ntot, &
              Zero,Tmp5,ntot)

counter = 1
do i=1,ntot
  do j=1,ntot
    write(molet_lu,FMT1) i,j,tmp5(counter)
    counter = counter + 1
  end do
end do

Close(molet_lu)

call mma_deallocate(OvMat)
call mma_deallocate(tmp1)
call mma_deallocate(tmp3)
call mma_deallocate(tmp4)
call mma_deallocate(tmp5)
call mma_deallocate(F1_PDFT)



end subroutine molel_write
