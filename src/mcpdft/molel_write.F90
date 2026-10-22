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
! Copyright (C) 2025, Andrew M. Sand                                   *
!***********************************************************************

subroutine molel_write()

use onedat, only: snoori
use rasscf_global, only: lRoots
use general_data, only: nactel, nfro, nish, nash, nbas, ntot, ntot1, ntot2, nsym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u0
use Constants, only: One, Zero

implicit none
#include "warnings.h"
integer(kind=iwp) :: comp, opt, rc, sylbl, counter, i, ibas, tot_el
integer(kind=iwp) :: isym, j, nactorb, nAOs, offset
character(len=8) :: label
character(len=80) :: FMT1
real(kind=wp), allocatable :: int1e_ovlp(:), OvMat(:), CMO(:), F1_PDFT(:), tmp1(:), tmp2(:), &
tmp3(:), tmp4(:)


!overlap matrix - AO
rc = -1
opt = ibset(0,sNoOri)
comp = 1
sylbl = 1
label = 'Mltpl  0'

call mma_allocate(int1e_ovlp,ntot1+4,label='int1e_ovlp')
call mma_allocate(OvMat,ntot**2,label='OvMat')
call rdone(rc,opt,label,comp,int1e_ovlp,sylbl)
!Open(unit=1099,file='MolEl.dat', action='write',iostat=ios)
call molcas_open(1099,'MolEl.dat')
nAOs = 0
do isym=1,nsym
  nAOs = nAOs + nBas(iSym)
end do
Tot_El = 0
nActOrb = 0
Do iSym=1,nSym
   Tot_El=Tot_El+ 2*(nFro(iSym)+nIsh(iSym))
   nActOrb = nActOrb + nAsh(iSym)
End Do
Tot_El = Tot_El + nActEl

FMT1 = "(I4,1X,I4,1X,I4,1X,I4,1X,I4)"
write(1099,*) 'Number of states,orbitals,electrons,ActOrb,ActEl'
write(1099,FMT1) lroots,nAOs,Tot_El,nActOrb,nActEl
write(1099,*) 'Overlap Matrix (AO)'

FMT1 = "(I4,1X,I4,1X,F20.14)"
counter = 0
offset = 0
do isym=1,nsym
   iBas = nBas(iSym)
   do i=1,iBas
      do j=1,i
         write(1099,FMT1) i+offset,j+offset,int1e_ovlp(1+counter)
         OvMat(1+(i+offset-1)*ntot+(j+offset-1)) = int1e_ovlp(1+counter)
         OvMat(1+(j+offset-1)*ntot+(i+offset-1)) = int1e_ovlp(1+counter)
         counter = counter + 1
      end do
   end do
   offset = offset + iBas
end do
call mma_deallocate(int1e_ovlp)



call mma_allocate(CMO,NTOT2,Label='CMO')
call ReadVc_m(CMO)
write(1099,*) 'Molecular orbital coefficients'
counter = 0
offset = 0
do isym=1,nsym
   iBas = nBas(iSym)
   do i=1,iBas
      do j=1,iBas
         write(1099,FMT1) i+offset,j+offset,CMO(1+counter)
         counter = counter +1
      end do
   end do
   offset = offset + iBas
end do



!Read this in, but it is in MO - convert to AO
!call mma_allocate(F1_PDFT,NTOT1,Label='F1p')
!call get_darray('F1_PDFT',FMO1t,ntot1)
      write(1099,*) 'Effective Hamiltonian(s) (AO)'

!S * C * F_mo * C^T * S
call mma_allocate(F1_PDFT,NTOT1,Label='F1p')
call mma_allocate(tmp1,ntot**2,Label='tmp1')
call mma_allocate(tmp2,ntot**2,Label='tmp2')
call mma_allocate(tmp3,ntot**2,Label='tmp3')
call mma_allocate(tmp4,ntot**2,Label='tmp4')


  call DGEMM_('N','N',ntot,ntot,ntot, &
              One,OvMat,ntot, &
              CMO,ntot, &
              Zero,tmp1,ntot)
  call DGEMM_('T','N',ntot,ntot,ntot, &
              One,CMO,ntot, &
              OvMat,ntot, &
              Zero,Tmp2,ntot)

call get_darray('F1_PDFT',F1_PDFT,ntot1)

call Square(F1_PDFT,tmp3,1,ntot,ntot)

  call DGEMM_('N','N',ntot,ntot,ntot, &
              One,Tmp1,ntot, &
              Tmp3,ntot, &
              Zero,Tmp4,ntot)
  call DGEMM_('N','N',ntot,ntot,ntot, &
              One,Tmp4,ntot, &
              Tmp2,ntot, &
              Zero,Tmp1,ntot)
counter = 1
do i=1,ntot
  do j=1,ntot
    write(1099,FMT1) i,j,tmp1(counter)
    counter = counter + 1
  end do
end do

Close(1099)

call mma_deallocate(OvMat)
call mma_deallocate(CMO)
call mma_deallocate(tmp1)
call mma_deallocate(tmp2)
call mma_deallocate(tmp3)
call mma_deallocate(tmp4)
call mma_deallocate(F1_PDFT)

end subroutine molel_write
