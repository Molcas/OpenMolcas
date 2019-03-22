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
      ! The initialization subroutine lives outside of the
      ! faroald module so that it can be called separately.
      ! It fills in the necessary global data in the module,
      ! in preparation for the sigma_update routine.
      subroutine faroald_init(nactel,nasht,ispin)
      use faroald
      use second_quantization
      implicit none

      integer, intent(in) :: nactel, nasht, ispin

      my_nel = nactel
      my_norb = nasht
      mult = ispin

      ! compute alpha/beta subsets
      nela = (my_nel+(mult-1))/2
      nhoa = my_norb - nela
      nelb = (my_nel-(mult-1))/2
      nhob = my_norb - nelb
      ndeta = binom_coef(nela,my_norb)
      ndetb = binom_coef(nelb,my_norb)
      my_ndet = ndeta * ndetb

      ! initialize binomial tables needed for ranking
      call rank_init

      ! number of possible non-desctructive single excitations
      max_ex1a = nela * (nhoa + 1)
      max_ex1b = nelb * (nhob + 1)

!     The number of non-vanishing elements in F is the number of
!     non-vanashing single and double excitations from a string, and
!     computed as the number of strings which differ in two, one, or zero
!     spin orbitals from a given string.
      max_ex2a = (nela * (nela - 1) * (nhoa) * (nhoa - 1))/4 +
     & nela * nhoa + 1
      max_ex2b = (nelb * (nelb - 1) * (nhob) * (nhob - 1))/4 +
     & nelb * nhob + 1

!     For sigma1, we use an excitation table for the beta strings.
      allocate (ex1_b(max_ex1b,ndetb))
      call ex1_init(nelb, my_norb, ex1_b)

!     For sigma2, we use an excitation table for the alpha strings, unless
!     Ms = 0 (mult = 1), as for singlet spins sigma2 is not computed and
!     just taken as the transpose of sigma1.
      if (mult.ne.1) then
        allocate (ex1_a(max_ex1a,ndeta))
        call ex1_init(nela, my_norb, ex1_a)
      end if

!     For sigma3, we need to construct a list of determinant couples L(i)
!     = E_pq R(i) for a given p and q for the alpha strings. This is done
!     by a subroutine that takes L, R, and sgn arrays. The maximum size of
!     these arrays is when p=q, at which point the number of couples is:
      max_LRs = binom_coef(nela-1,my_norb-1)

      end subroutine
