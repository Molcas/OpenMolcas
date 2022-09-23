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

module ChT3_global

! no              - # of occupied orbitals
! nv              - # of virtual orbitals
! nc              - # of Cholesky vectors
! nfr             - # if frozen orbitals
! DimGrpaR        - Dimensions of group (a,b) for routines evaluating o3v3 contributions
! [LT][12]Name    - Names of all used files
! TWall, TCpu     - pomocne
! TWall0, TCpu0   - zaciatocne hodnoty
! TWall_l, TCpu_l - lokalne hodnoty

! NBLOCK is used for direct-access unformatted I/O via multi_*.
! Count a record size in real(wp) words. Should be 2**k.

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: nblock = 2**11
character, parameter :: ICH(3) = ['A','B','C']

integer(kind=iwp) :: IOPT(2), IT, LunAux, maxdim, nc, nfr, NNOAB(3), NNUAB(3), no, NOAB(2), NUAB(2), nv, NvGrp, printkey, &
                     t3_starta, t3_startb, t3_stopa, t3_stopb
real(kind=wp) :: TCpu, TCpu_l, TCpu0, TWall, TWall_l, TWall0
logical(kind=iwp) :: gen_files, run_triples
integer(kind=iwp), allocatable :: DimGrpaR(:)
character(len=6), allocatable :: L1Name(:), L2Name(:,:), T2Name(:,:)

public :: DimGrpaR, gen_files, ICH, IOPT, IT, L1Name, L2Name, LunAux, maxdim, nblock, nc, nfr, NNOAB, NNUAB, no, NOAB, NUAB, nv, &
          NvGrp, printkey, run_triples, T2Name, t3_starta, t3_startb, t3_stopa, t3_stopb, TCpu, TCpu_l, TCpu0, TWall, TWall_l, &
          TWall0

end module ChT3_global
