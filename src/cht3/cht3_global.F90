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

! maxGrp          - maximal # of Groups
! no              - # of occupied orbitals
! nv              - # of virtual orbitals
! nc              - # of Cholesky vectors
! nfr             - # if frozen orbitals
! DimGrpaR        - Dimensions of group (a,b) for routines evaluating o3v3 contributions
! [LT][12]Name    - Names of all used files
! NChLoc          - local number of Cholesky vectors on each node
! TWall, TCpu     - pomocne
! TWall0, TCpu0   - zaciatocne hodnoty
! TWall_l, TCpu_l - lokalne hodnoty

! NBLOCK is used for direct-access unformatted I/O via multi_*.
! Count a record size in real(wp) words. Should be 2**k.

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: maxGrp = 32, MaxNod = 512, nblock = 2**11, ndupmx = 20
character, parameter :: ICH(3) = ['A','B','C']

integer(kind=iwp) :: DimGrpaR(maxGrp), dupblk(ndupmx), IOPT(2), IT, LunAux, maxdim, NChLoc(0:MaxNod-1), nc, ndup, nfr, NNOAB(3), &
                     NNUAB(3), no, NOAB(2), NUAB(2), nv, NvGrp, printkey, t3_starta, t3_startb, t3_stopa, t3_stopb
real(kind=wp) :: TCpu, TCpu_l, TCpu0, TWall, TWall_l, TWall0
logical(kind=iwp) :: gen_files, run_triples
character(len=132) :: dupfil(ndupmx)
character(len=6) :: L1Name(maxGrp), L2Name(maxGrp,maxGrp), T2Name(maxGrp,maxGrp)

public :: DimGrpaR, dupblk, dupfil, gen_files, ICH, IOPT, IT, L1Name, L2Name, LunAux, maxdim, maxGrp, nblock, nc, NChLoc, ndup, &
          ndupmx, nfr, NNOAB, NNUAB, no, NOAB, NUAB, nv, NvGrp, printkey, run_triples, T2Name, t3_starta, t3_startb, t3_stopa, &
          t3_stopb, TCpu, TCpu_l, TCpu0, TWall, TWall_l, TWall0

end module ChT3_global
