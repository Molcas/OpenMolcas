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

implicit none
integer IFSBTAB1(*), IFSBTAB2(*)
integer ISSTAB(*), MAPORB(*), NRT2M
real*8 DET1(*), DET2(*)
logical IF21, IF12
integer NSRT2M
real*8 RT2M(NRT2M)
integer SPIN, OrbTab(*)
integer NASHT, NASORB
real*8 GVAL, GAAA, GAAB, GABA, GBAB, GBBA, GBBB
integer IAJBLA, IBJALB
integer IAJALA, IAJALB, IBJBLA, IBJBLB
integer LORB, JORB, IORB
integer JORBA, JORBB, LORBA, LORBB, IORBA, IORBB
integer ITABS, JTABS, LTABS, JLTABS, IJLTABS
real*8, allocatable :: SRT2M(:)

! Pick out nr of active orbitals from orbital table:

NASORB = OrbTab(4)
NASHT = NASORB/2
!NASGEM = (NASORB*(NASORB-1))/2
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
        if (SPIN == 1) then ! K^1/2,1/2= (AAB+BBB)
          GAAB = SRT2M(IAJALB)
          GBBB = SRT2M(IBJBLB)
          GVAL = (GAAB+GBBB)
        else if (SPIN == -1) then ! SDA. K^1/2,-1/2
          GBBA = SRT2M(IBJBLA)
          GAAA = SRT2M(IAJALA)
          GVAL = (GAAA+GBBA)
        else if (SPIN == 2) then
          GBBB = SRT2M(IBJBLB)
          GVAL = GBBB
        else if (SPIN == 3) then
          GAAA = SRT2M(IAJALA)
          GVAL = GAAA
        else if (SPIN == 4) then
          GAAB = SRT2M(IAJALB)
          GVAL = GAAB
        else if (SPIN == 5) then
          GBBA = SRT2M(IBJBLA)
          GVAL = GBBA
        else if (SPIN == 6) then
          GABA = SRT2M(IAJBLA)
          GVAL = GABA
        else if (SPIN == 7) then
          GBAB = SRT2M(IBJALB)
          GVAL = GBAB
        end if
      else if (JORB == LORB) then
        IAJALB = IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJBLA = IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        IBJBLB = IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
        IBJALB = IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJALA = IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
        IBJBLA = IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        if (SPIN == 1) then ! K^1/2,1/2
          GAAB = SRT2M(IAJALB)
          GBBB = SRT2M(IBJBLB)
          GVAL = (GAAB+GBBB)
        else if (SPIN == -1) then ! SDA. K^1/2,-1/2
          GBBA = SRT2M(IBJBLA)
          GAAA = SRT2M(IAJALA)
          GVAL = (GAAA+GBBA)
        else if (SPIN == 2) then
          GBBB = SRT2M(IBJBLB)
          GVAL = GBBB
        else if (SPIN == 3) then
          GAAA = SRT2M(IAJALA)
          GVAL = GAAA
        else if (SPIN == 4) then
          GAAB = SRT2M(IAJALB)
          GVAL = GAAB
        else if (SPIN == 5) then
          GBBA = SRT2M(IBJBLA)
          GVAL = GBBA
        else if (SPIN == 6) then
          GABA = SRT2M(IAJBLA)
          GVAL = GABA
        else if (SPIN == 7) then
          GBAB = SRT2M(IBJALB)
          GVAL = GBAB
        end if
      else if (JORB < LORB) then ! When J<L
        IAJALB = IORBA+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJBLA = IORBA+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        IBJBLB = IORBB+NASORB*(NASORB*(JORBB-1)+LORBB-1)
        IBJALB = IORBB+NASORB*(NASORB*(JORBA-1)+LORBB-1)
        IAJALA = IORBA+NASORB*(NASORB*(JORBA-1)+LORBA-1)
        IBJBLA = IORBB+NASORB*(NASORB*(JORBB-1)+LORBA-1)
        if (SPIN == 1) then ! K^1/2,1/2
          GAAB = SRT2M(IAJALB)
          GBBB = SRT2M(IBJBLB)
          GVAL = (GAAB+GBBB)
        else if (SPIN == -1) then ! SDA. K^1/2,-1/2
          GBBA = SRT2M(IBJBLA)
          GAAA = SRT2M(IAJALA)
          GVAL = (GAAA+GBBA)
        else if (SPIN == 2) then
          GBBB = SRT2M(IBJBLB)
          GVAL = GBBB
        else if (SPIN == 3) then
          GAAA = SRT2M(IAJALA)
          GVAL = GAAA
        else if (SPIN == 4) then
          GAAB = SRT2M(IAJALB)
          GVAL = GAAB
        else if (SPIN == 5) then
          GBBA = SRT2M(IBJBLA)
          GVAL = GBBA
        else if (SPIN == 6) then
          GABA = SRT2M(IAJBLA)
          GVAL = GABA
        else if (SPIN == 7) then
          GBAB = SRT2M(IBJALB)
          GVAL = GBAB
        end if
      end if
      IJLTABS = ITABS+NASHT*(JLTABS-1)
      RT2M(IJLTABS) = GVAL
    end do
  end do
end do

call mma_deallocate(SRT2M)

end subroutine MKRTDM2
