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
! Copyright (C) 2020, Bruno Tenorio                                    *
!***********************************************************************

subroutine MKRTDM2(IFSBTAB1,IFSBTAB2,ISSTAB,MAPORB,DET1,DET2,IF21,IF12,NRT2M,RT2M,SPIN,OrbTab)
! Given two CI expansions, using a biorthonormal set of SD's,
! calculate the 2-particle transition density matrix
! in the biorthonormal active orbital basis.
! It will build the contribution from high spin (J->beta,L->beta)
! and low spin (J->beta,L->alpha).
! The spin coupling matrix elements have the following index-code:
!SPIN=1 means  K2V (AAB+BBB)
!SPIN=-1 means SDA (AAA+BBA)
!SPIN=2 means: bbb
!SPIN=3 means: aaa
!SPIN=4 means: aab
!SPIN=5 means: bba
!SPIN=6 means: aba
!SPIN=7 means: bab
! Notice, SPIN here has nothing to do with the spin quantum number. It
! is just a printing code.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IFSBTAB1(*), IFSBTAB2(*), ISSTAB(*), MAPORB(*), NRT2M, SPIN, OrbTab(*)
real(kind=wp), intent(in) :: DET1(*), DET2(*)
logical(kind=iwp), intent(in) :: IF21, IF12
real(kind=wp), intent(out) :: RT2M(NRT2M)
integer(kind=iwp) :: IAJALA, IAJALB, IAJBLA, IBJALB, IBJBLA, IBJBLB, IJLTABS, IORB, IORBA, IORBB, ITABS, JLTABS, JORB, JORBA, &
                     JORBB, JTABS, LORB, LORBA, LORBB, LTABS, NASHT, NASORB, NSRT2M
real(kind=wp) :: GAAA, GAAB, GABA, GBAB, GBBA, GBBB, GVAL
real(kind=wp), allocatable :: SRT2M(:)

! Pick out nr of active orbitals from orbital table:

NASORB = OrbTab(4)
NASHT = NASORB/2
!NASGEM = nTri_Elem(NASORB-1)
NSRT2M = NASORB**3
call mma_allocate(SRT2M,nSRT2M,Label='SRT2M')
SRT2M(:) = Zero
call SRTDM2(ORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,DET1,DET2,IF21,IF12,SRT2M)

! Mapping from active spin-orbital to active orbital in external order.
! Note that these differ, not just because of the existence of two
! spin-orbitals for each orbital, but also because the active orbitals
! (external order) are grouped by symmetry and then RAS space, but the
! spin orbitals are grouped by subpartition.
GVAL = Zero
IAJALA = 0 ! dummy initialize
IAJALB = 0 ! dummy initialize
IAJBLA = 0 ! dummy initialize
IBJALB = 0 ! dummy initialize
IBJBLA = 0 ! dummy initialize
IBJBLB = 0 ! dummy initialize

! For high spin density it will keep only beta,beta,beta.
! Notice that beta,beta,beta when J=L is zero.
! For low spin density it'll keep only alpha,beta,alpha.

do IORB=1,NASHT
  IORBA = 2*IORB-1
  IORBB = 2*IORB
  ITABS = MAPORB(IORBA)
  do JORB=1,NASHT
    JORBA = 2*JORB-1
    JORBB = 2*JORB
    JTABS = MAPORB(JORBA)
    do LORB=1,NASHT
      LORBA = 2*LORB-1
      LORBB = 2*LORB
      LTABS = MAPORB(LORBA)
      JLTABS = LTABS+NASHT*(JTABS-1)
      if (JORB > LORB) then ! When J>L
        IAJALB = IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJBLA = IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        IBJBLB = IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
        IBJALB = IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJALA = IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
        IBJBLA = IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        select case (SPIN)
          case (1)
            ! K^1/2,1/2= (AAB+BBB)
            GAAB = SRT2M(IAJALB)
            GBBB = SRT2M(IBJBLB)
            GVAL = (GAAB+GBBB)
          case (-1)
            ! SDA. K^1/2,-1/2
            GBBA = SRT2M(IBJBLA)
            GAAA = SRT2M(IAJALA)
            GVAL = (GAAA+GBBA)
          case (2)
            GBBB = SRT2M(IBJBLB)
            GVAL = GBBB
          case (3)
            GAAA = SRT2M(IAJALA)
            GVAL = GAAA
          case (4)
            GAAB = SRT2M(IAJALB)
            GVAL = GAAB
          case (5)
            GBBA = SRT2M(IBJBLA)
            GVAL = GBBA
          case (6)
            GABA = SRT2M(IAJBLA)
            GVAL = GABA
          case (7)
            GBAB = SRT2M(IBJALB)
            GVAL = GBAB
        end select
      else if (JORB == LORB) then
        IAJALB = IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJBLA = IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        IBJBLB = IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
        IBJALB = IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJALA = IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
        IBJBLA = IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        select case (SPIN)
          case (1)
            ! K^1/2,1/2
            GAAB = SRT2M(IAJALB)
            GBBB = SRT2M(IBJBLB)
            GVAL = (GAAB+GBBB)
          case (-1)
            ! SDA. K^1/2,-1/2
            GBBA = SRT2M(IBJBLA)
            GAAA = SRT2M(IAJALA)
            GVAL = (GAAA+GBBA)
          case (2)
            GBBB = SRT2M(IBJBLB)
            GVAL = GBBB
          case (3)
            GAAA = SRT2M(IAJALA)
            GVAL = GAAA
          case (4)
            GAAB = SRT2M(IAJALB)
            GVAL = GAAB
          case (5)
            GBBA = SRT2M(IBJBLA)
            GVAL = GBBA
          case (6)
            GABA = SRT2M(IAJBLA)
            GVAL = GABA
          case (7)
            GBAB = SRT2M(IBJALB)
            GVAL = GBAB
        end select
      else if (JORB < LORB) then ! When J<L
        IAJALB = IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJBLA = IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        IBJBLB = IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
        IBJALB = IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJALA = IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
        IBJBLA = IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        select case (SPIN)
          case (1)
            ! K^1/2,1/2
            GAAB = SRT2M(IAJALB)
            GBBB = SRT2M(IBJBLB)
            GVAL = (GAAB+GBBB)
          case (-1)
            ! SDA. K^1/2,-1/2
            GBBA = SRT2M(IBJBLA)
            GAAA = SRT2M(IAJALA)
            GVAL = (GAAA+GBBA)
          case (2)
            GBBB = SRT2M(IBJBLB)
            GVAL = GBBB
          case (3)
            GAAA = SRT2M(IAJALA)
            GVAL = GAAA
          case (4)
            GAAB = SRT2M(IAJALB)
            GVAL = GAAB
          case (5)
            GBBA = SRT2M(IBJBLA)
            GVAL = GBBA
          case (6)
            GABA = SRT2M(IAJBLA)
            GVAL = GABA
          case (7)
            GBAB = SRT2M(IBJALB)
            GVAL = GBAB
        end select
      end if
      IJLTABS = ITABS+NASHT*(JLTABS-1)
      RT2M(IJLTABS) = GVAL
    end do
  end do
end do

call mma_deallocate(SRT2M)

end subroutine MKRTDM2
