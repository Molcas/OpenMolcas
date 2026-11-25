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

subroutine Faro_Init(nactel,nasht,ispin,full)
! The initialization subroutine lives outside of the
! Faroald module so that it can be called separately.
! It fills in the necessary global data in the module,
! in preparation for the Faro_sigma_update routine.
! full: also the excitation tables and occupation patterns of the full
! (CAS) string and determinant spaces (Faro_sigma_update, densities, traci);
! not needed for the string-driven RAS code (Faroald_RAS), which only uses
! the sizes, the rank tables and the integrals.

use Faroald, only: Faro_ex1_a, Faro_ex1_b, Faro_ex1_init, Faro_Full, Faro_gtuvx, Faro_htu, Faro_max_ex1a, Faro_max_ex1b, &
                   Faro_max_ex2a, Faro_max_ex2b, Faro_max_LRs, Faro_mult, &
                   Faro_nDet, Faro_nEl, Faro_nOrb, Faro_ndeta, Faro_ndetb, Faro_nela, Faro_nelb, Faro_nhoa, Faro_nhob
use Faroald, only: mma_allocate ! with extensions for ex1_struct
!use Faroald, only: Faro_verify_occ_patterns, Faro_build_patterns, Faro_analyse_patterns, build_conf_arcw
use Faroald, only: Faro_verify_occ_patterns, Faro_build_patterns, Faro_analyse_patterns
use Faroald_SecQuant, only: Faro_Binom, Faro_rank_init
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nactel, nasht, ispin
logical(kind=iwp), intent(in) :: full

Faro_nEl = nactel
Faro_nOrb = nasht
Faro_mult = ispin

! compute alpha/beta subsets
Faro_nela = (Faro_nEl+(Faro_mult-1))/2
Faro_nhoa = Faro_nOrb-Faro_nela
Faro_nelb = (Faro_nEl-(Faro_mult-1))/2
Faro_nhob = Faro_nOrb-Faro_nelb
Faro_ndeta = Faro_Binom(Faro_nela,Faro_nOrb)
Faro_ndetb = Faro_Binom(Faro_nelb,Faro_nOrb)
Faro_nDet = Faro_ndeta*Faro_ndetb

! initialize binomial tables needed for ranking
call Faro_rank_init()

! number of possible non-destructive single excitations
Faro_max_ex1a = Faro_nela*(Faro_nhoa+1)
Faro_max_ex1b = Faro_nelb*(Faro_nhob+1)

! The number of non-vanishing elements in F is the number of
! non-vanashing single and double excitations from a string, and
! computed as the number of strings which differ in two, one, or zero
! spin orbitals from a given string.
Faro_max_ex2a = (Faro_nela*(Faro_nela-1)*(Faro_nhoa)*(Faro_nhoa-1))/4+Faro_nela*Faro_nhoa+1
Faro_max_ex2b = (Faro_nelb*(Faro_nelb-1)*(Faro_nhob)*(Faro_nhob-1))/4+Faro_nelb*Faro_nhob+1

Faro_Full = full
if (full) then

! For SG_Epq_Psi, we use an excitation table for the beta strings.
call mma_allocate(Faro_ex1_b,Faro_max_ex1b,Faro_ndetb,label='ex1_b')
call Faro_ex1_init(Faro_nelb,Faro_nOrb,Faro_ex1_b)

! For sigma2, we use an excitation table for the alpha strings, unless
! Ms = 0 (Faro_mult = 1), as for singlet spins sigma2 is not computed and
! just taken as the transpose of SG_Epq_Psi.
!if (Faro_mult /= 1) then
call mma_allocate(Faro_ex1_a,Faro_max_ex1a,Faro_ndeta,label='ex1_a')
call Faro_ex1_init(Faro_nela,Faro_nOrb,Faro_ex1_a)
!end if

call Faro_verify_occ_patterns()
call Faro_build_patterns()
call Faro_analyse_patterns()

end if

! For sigma3, we need to construct a list of determinant couples L(i)
! = E_pq R(i) for a given p and q for the alpha strings. This is done
! by a subroutine that takes L, R, and sgn arrays. The maximum size of
! these arrays is when p=q, at which point the number of couples is:
Faro_max_LRs = Faro_Binom(Faro_nela-1,Faro_nOrb-1)

call mma_allocate(Faro_htu,Faro_nOrb,Faro_nOrb,label='htu')
call mma_allocate(Faro_gtuvx,Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb,label='gtuvx')

end subroutine Faro_Init
