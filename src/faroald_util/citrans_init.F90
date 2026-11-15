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
subroutine citrans_init(nel,norb,mult,nsym,osym,lsym,nrsprt,oras,nrasel)

use second_quantization, only: binom_coef
use citrans, only: citrans_graph_init, citrans_sym_init, ndet_group, ncsf_group, ndo_max, ndo_min, ndoc_group, nsoc_group, &
                   spintable_create, &
                   spintabs, spintabs_allocate
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nel, norb, mult
integer(kind=iwp), intent(in) :: nsym, osym(norb), lsym, nrsprt, oras(norb), nrasel(nrsprt)
integer(kind=iwp) :: ido, iso, nela, nelb

! compute alpha/beta subsets
nela = (nel+(mult-1))/2
nelb = (nel-(mult-1))/2

! determine the range of the configuration groups
if (nel > norb) then
  ndo_min = nel-norb
else
  ndo_min = 0
end if
ndo_max = nelb

! compute the various sizes per group
call mma_allocate(ndoc_group,[ndo_min,ndo_max],label='ndoc_group')
call mma_allocate(nsoc_group,[ndo_min,ndo_max],label='nsoc_group')
call mma_allocate(ndet_group,[ndo_min,ndo_max],label='ndet_group')
call mma_allocate(ncsf_group,[ndo_min,ndo_max],label='ncsf_group')

call spintabs_allocate()
! loop over configurations
do ido=ndo_min,ndo_max
  iso = nel-2*ido
  ! compute different block sizes
  ndoc_group(ido) = binom_coef(ido,norb)
  nsoc_group(ido) = binom_coef(iso,norb-ido)
  ndet_group(ido) = binom_coef(nela-ido,iso)
  ncsf_group(ido) = ndet_group(ido)-binom_coef(nela-ido+1,iso)
  ! compute+store spin table for this configuration
  spintabs(ido)%ndet = ndet_group(ido)
  spintabs(ido)%ncsf = ncsf_group(ido)
  call spintable_create(iso,nelb-ido,spintabs(ido))
end do

! the file order of the CSFs (citrans_sort)
call citrans_sym_init(nel,norb,nsym,osym,lsym,nrsprt,oras,nrasel)
call citrans_graph_init(nel,norb,mult)

end subroutine citrans_init
