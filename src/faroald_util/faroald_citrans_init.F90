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

!> nsym, osym (the symmetries of the norb active orbitals) and lsym (the
!> symmetry of the CSF space); C1 is nsym = lsym = 1, osym = 1. RAS: nrsprt
!> partitions, oras the partition of each orbital, nrasel(k) the minimum
!> number of electrons in the partitions 1..k; the orbitals must be the
!> levels of the split graph (by partition, within each by symmetry). CAS
!> is nrsprt = 1, oras = 1, nrasel = nel.
subroutine Faro_CITrans_Init(nel,norb,Faro_mult,nsym,osym,lsym,nrsprt,oras,nrasel,nrasmx)

use Faroald_SecQuant, only: Faro_Binom
use Faroald_CITrans, only: Faro_graph_init, Faro_sym_init, Faro_ndet_group, Faro_ncsf_group, Faro_ndo_max, Faro_ndo_min, &
                           Faro_ndoc_group, Faro_nsoc_group, &
                   Faro_spintable_create, &
                   Faro_spintabs, Faro_spintabs_allocate
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nel, norb, Faro_mult
integer(kind=iwp), intent(in) :: nsym, osym(norb), lsym, nrsprt, oras(norb), nrasel(nrsprt), nrasmx(nrsprt)
integer(kind=iwp) :: ido, iso, Faro_nela, Faro_nelb

! compute alpha/beta subsets
Faro_nela = (nel+(Faro_mult-1))/2
Faro_nelb = (nel-(Faro_mult-1))/2

! determine the range of the configuration groups
if (nel > norb) then
  Faro_ndo_min = nel-norb
else
  Faro_ndo_min = 0
end if
Faro_ndo_max = Faro_nelb

! compute the various sizes per group
call mma_allocate(Faro_ndoc_group,[Faro_ndo_min,Faro_ndo_max],label='ndoc_group')
call mma_allocate(Faro_nsoc_group,[Faro_ndo_min,Faro_ndo_max],label='nsoc_group')
call mma_allocate(Faro_ndet_group,[Faro_ndo_min,Faro_ndo_max],label='ndet_group')
call mma_allocate(Faro_ncsf_group,[Faro_ndo_min,Faro_ndo_max],label='ncsf_group')

call Faro_spintabs_allocate()
! loop over configurations
do ido=Faro_ndo_min,Faro_ndo_max
  iso = nel-2*ido
  ! compute different block sizes
  Faro_ndoc_group(ido) = Faro_Binom(ido,norb)
  Faro_nsoc_group(ido) = Faro_Binom(iso,norb-ido)
  Faro_ndet_group(ido) = Faro_Binom(Faro_nela-ido,iso)
  Faro_ncsf_group(ido) = Faro_ndet_group(ido)-Faro_Binom(Faro_nela-ido+1,iso)
  ! compute+store spin table for this configuration
  Faro_spintabs(ido)%ndet = Faro_ndet_group(ido)
  Faro_spintabs(ido)%ncsf = Faro_ncsf_group(ido)
  call Faro_spintable_create(iso,Faro_nelb-ido,Faro_spintabs(ido))
end do

! the file order of the CSFs (Faro_sort)
call Faro_sym_init(nel,norb,nsym,osym,lsym,nrsprt,oras,nrasel,nrasmx)
call Faro_graph_init(nel,norb,Faro_mult)

end subroutine Faro_CITrans_Init
