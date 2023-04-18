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

!0      maximal # of Groups, SubGroups
!       maxGrp, maxSGrp
!
!1.1    Basic parameters for CCSD procedure
!
!       no - # of occupied orbitals
!       nv - # of virtual orbitals
!       nc - # of Cholesky vectors
!       nfr- # if frozen orbitals
!
!       mhkey - FTN/BLAS switch
!       maxiter - Maximal number of iterations
!       restkey - key for restart [1/0]
!       generkey - key for generation of integrals from L vectors [1/0]
!       intkey - key for generation of (vv|vv) and (vv|vo) integrals [1/0]
!       W34DistKey - 1 - only required W34 files are stored on each node
!                    0 - all W34 files are stored on each node
!       JoinLkey - key to control the cummulation of L L(ml)->L(m)
!                  0 - no cummulation of L files are needed
!                  1 - cummulation immidiately afret L(ml) are produced
!                      (before I, -> I and W are generated using full L)
!                  2 - cummulation after I
!                      (I are generated per partes using L(ml) + gadgop
!                       and W are generated using full L)
!                  3 - cummulation after W
!                      (also W34 are generated per partes using L(ml) +g)
!
!       printkey -  1 - Only basic info is printed (default)
!                   2 - 1 + timings
!                  10 - debug info
!
!1.2    position of Permanent arrays
!
!       Fock matrix             - PosFoo,PosFvv,PosFvo
!       Orbital energies        - PosOE
!       T1 amplitudes, old, new - PosT1o,PosT1n
!       N2 Intermediates for T1 - PosHoo,PosHvv,PosHvo
!       N2 Intermediates for T2 - PosGoo,PosGvv
!       O2OO Intermediate A     - PosA - @@ na zamyslenie, ci nie cez worky
!       O2OO Intermediate Aex   - PosAex - @@ na zamyslenie, ci nie cez worky
!       position of Free space  - PosFree
!
!1.3    Checkeroo part
!       T1c, T2c, OEo, OEv, Q0, Q1, Q21, Q22, Q3, Q4, L0k, L1k, L2k,
!       Jc, Kc, Ac, Bc, Hooc, Hvvc, Hvoc, Gooc, Gvvc
!
!2      for routines evaluating o3v3 contributions
!
!2.1    Dimensions of groups
!       DimGrpv
!
!2.2    Dimensions of groups (a,b)
!       DimGrpaR
!
!3      for routines evaluating o2v4 contributions
!
!3.1    Groups borders (which subgroup is low and upper limit of given group)
!       GrpaLow, GrpbeLow, GrpaUp, GrpbeUp
!
!3.2    Dimensions of groups
!       DimGrpa, DimGrpbe
!
!3.3    Dimensions of subgroups
!       DimSGrpa, DimSGrpbe
!
!4      Timing
!       TWall, TCpu     - pomocne
!       TWall0, TCpu0   - zaciatocne hodnoty
!       TWall_l, TCpu_l - lokalne hodnoty
!
!5      Names of all used files
!
!5.1    Names of L0-L2 files
!       L0Name, L1Name, L2Name
!
!5.2    Names of T2 files
!       T2Name
!
!5.3    Names of I0-3 files
!       I0Name, I1Name, I2Name, I3Name
!
!5.4    Names of Tmp1,2 files
!       Tmp1Name, Tmp2Name
!
!5.5    Names of Tmp3 file
!       Tmp3Name
!
!6      for parallel computing
!
!6.1    BetaID(myRank,beGrp) - indicates if contributions that need to be
!       realized only once per each beta' will be realized on this node (myRank)
!        =1, contributions will be realized on this node
!        =0, not to be realized on this node
!
!6.2.1  BeAID(myRank,beGrp,bGrp) - indicates which be' a' combinations
!       are to be evaluated on this node (myRank)
!        =1, this be' a' combination is to be realized on this node
!        =0, not to be realized on this node
!
!6.2.2  ABID(myRank,aGrp,bGrp) - indicates which a'b' values are to be evaluated
!       on this node (myRank)
!        =1, this a'b' combination is to be realized on this node
!        =0, not to be realized on this node
!
!6.3    Xyes(p,q) - indicates if in Tmp2(p,q) is a record of X from
!       o3v3chol part on this node
!        =1, there is a record from o3v3chol part on this node
!        =0, file is leaving empty from o3v3chol  on this node
!
!6.4    XYyes(p,q) - indicates if in Tmp2(p,q) is a record of XY from
!       o3v3t2 part on this node
!        =1, there is a record from o3v3t2 part on this node
!        =0, there is no record from o3v3t2  on this node
!            (but there still can be record from o3v3chol on this node)
!
!6.5    T2o2v4yes(p,q) - indicates if in Tmp3(p,q) is a record of T2n from
!       o2v4 part on this node
!        =1, there is a record from o2v4 part on this node
!        =0, there is no record from o2v4  on this node
!
!6.6    Local number of Cholesky vectors on each node
!       NChLoc
!
!6.7    Files, containing information which W3 and W4 files are used on this node
!       InqW3, InqW4

module chcc_global

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: maxGrp = 32, maxSGrp = 64, &
                                maxNod = 512

integer(kind=iwp) :: ABID(0:maxNod-1,maxGrp,maxGrp), BeAID(0:maxNod-1,maxGrp,maxGrp), BetaID(0:maxNod-1,maxGrp), DimGrpa(maxGrp), &
                     DimGrpaR(maxGrp), DimGrpbe(maxGrp), DimGrpv(maxGrp), DimSGrpa(maxSGrp), DimSGrpbe(maxSGrp), generkey, &
                     GrpaLow(maxGrp), GrpaUp(maxGrp), GrpbeLow(maxGrp), GrpbeUp(maxGrp), intkey, JoinLkey, maxiter, mhkey, nc, &
                     NChLoc(0:maxnod-1), nfr, no, nv, PosA, PosAex, PosFoo, PosFree, PosFvo, PosFvv, PosGoo, PosGvv, PosHoo, &
                     PosHvo, PosHvv, PosOE, PosT1n, PosT1o, printkey, restkey, T2o2v4yes(maxSGrp,maxSGrp), W34DistKey, &
                     Xyes(maxGrp,maxGrp), XYyes(maxGrp,maxGrp)
real(kind=wp) :: conv, TCpu, TCpu0, TCpu_l, TWall, TWall0, TWall_l
character(len=6) :: I0Name, I1Name(maxGrp), I2Name(maxGrp,maxGrp), I3Name(maxGrp,maxGrp), L0Name, L1Name(maxGrp), &
                    L2Name(maxGrp,maxGrp), T2Name(maxGrp,maxGrp), Tmp1Name(maxGrp,maxGrp), Tmp2Name(maxGrp,maxGrp), &
                    Tmp3Name(maxSGrp,maxSGrp)
logical(kind=iwp) :: InqW3(maxSGrp*(maxSGrp+1)/2,maxSGrp), InqW4(maxSGrp*(maxSGrp+1)/2,maxSGrp*(maxSGrp+1)/2)
real(kind=wp), allocatable :: Ac(:,:,:,:), Bc(:,:,:,:), Gooc(:,:), Gvvc(:,:), Hooc(:,:), Hvoc(:,:), Hvvc(:,:), Jc(:,:,:,:), &
                              Kc(:,:,:,:), L0k(:,:,:), L1k(:,:,:), L2k(:,:,:), OEo(:), OEv(:), Q0(:,:,:,:), Q1(:,:,:,:), &
                              Q21(:,:,:,:), Q22(:,:,:,:), Q3(:,:,:,:), Q4(:,:,:,:), T1c(:,:), T2c(:,:,:,:)

public :: ABID, Ac, Bc, BeAID, BetaID, conv, deallocate_arrays, DimGrpa, DimGrpaR, DimGrpbe, DimGrpv, DimSGrpa, DimSGrpbe, &
          generkey, Gooc, GrpaLow, GrpaUp, GrpbeLow, GrpbeUp, Gvvc, Hooc, Hvoc, Hvvc, I0Name, I1Name, I2Name, I3Name, InqW3, &
          InqW4, intkey, Jc, JoinLkey, Kc, L0k, L0Name, L1k, L1Name, L2k, L2Name, maxGrp, maxiter, maxSGrp, mhkey, nc, NChLoc, &
          nfr, no, nv, OEo, OEv, PosA, PosAex, PosFoo, PosFree, PosFvo, PosFvv, PosGoo, PosGvv, PosHoo, PosHvo, PosHvv, PosOE, &
          PosT1n, PosT1o, printkey, Q0, Q1, Q21, Q22, Q3, Q4, restkey, T1c, T2c, T2Name, T2o2v4yes, TCpu, TCpu0, TCpu_l, Tmp1Name, &
          Tmp2Name, Tmp3Name, TWall, TWall0, TWall_l, W34DistKey, Xyes, XYyes

contains

subroutine deallocate_arrays()

  use stdalloc, only: mma_deallocate

  if (allocated(T1c)) call mma_deallocate(T1c)
  if (allocated(T2c)) call mma_deallocate(T2c)
  if (allocated(OEo)) call mma_deallocate(OEo)
  if (allocated(OEv)) call mma_deallocate(OEv)
  if (allocated(Q0)) call mma_deallocate(Q0)
  if (allocated(Q1)) call mma_deallocate(Q1)
  if (allocated(Q21)) call mma_deallocate(Q21)
  if (allocated(Q22)) call mma_deallocate(Q22)
  if (allocated(Q3)) call mma_deallocate(Q3)
  if (allocated(Q4)) call mma_deallocate(Q4)
  if (allocated(L0k)) call mma_deallocate(L0k)
  if (allocated(L1k)) call mma_deallocate(L1k)
  if (allocated(L2k)) call mma_deallocate(L2k)
  if (allocated(Jc)) call mma_deallocate(Jc)
  if (allocated(Kc)) call mma_deallocate(Kc)
  if (allocated(Ac)) call mma_deallocate(Ac)
  if (allocated(Bc)) call mma_deallocate(Bc)
  if (allocated(Hooc)) call mma_deallocate(Hooc)
  if (allocated(Hvvc)) call mma_deallocate(Hvvc)
  if (allocated(Hvoc)) call mma_deallocate(Hvoc)
  if (allocated(Gooc)) call mma_deallocate(Gooc)
  if (allocated(Gvvc)) call mma_deallocate(Gvvc)

end subroutine deallocate_arrays

end module chcc_global
