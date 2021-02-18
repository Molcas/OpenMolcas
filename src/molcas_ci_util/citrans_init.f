************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine citrans_init(nel,norb,mult)
      use second_quantization
      use citrans
      implicit none

      integer, intent(in) :: nel, norb, mult
      integer :: nela, nelb

      integer :: ido, iso

      ! compute alpha/beta subsets
      nela = (nel+(mult-1))/2
      nelb = (nel-(mult-1))/2

      ! determine the range of the configuration groups
      if (nel.gt.norb) then
        ndo_min = nel - norb
      else
        ndo_min = 0
      end if
      ndo_max = nelb

      ! compute the various sizes per group
      allocate (ndoc_group(ndo_min:ndo_max))
      allocate (nsoc_group(ndo_min:ndo_max))
      allocate (ndet_group(ndo_min:ndo_max))
      allocate (ncsf_group(ndo_min:ndo_max))

      allocate (spintabs(ndo_min:ndo_max))
      ! loop over configurations
      do ido = ndo_min, ndo_max
        iso = nel - 2 * ido
        ! compute different block sizes
        ndoc_group(ido) = binom_coef(ido,norb)
        nsoc_group(ido) = binom_coef(iso,norb-ido)
        ndet_group(ido) = binom_coef(nela-ido,iso)
        ncsf_group(ido) = ndet_group(ido) - binom_coef(nela-ido+1,iso)
        ! compute+store spin table for this configuration
        spintabs(ido)%ndet = ndet_group(ido)
        spintabs(ido)%ncsf = ncsf_group(ido)
        call spintable_create(iso,nelb-ido,spintabs(ido))
      end do

      end subroutine
