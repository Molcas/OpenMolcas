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

#include "macros.fh"

subroutine MKTDM2(LSYM1,MPLET1,MSPROJ1,IFSBTAB1,LSYM2,MPLET2,MSPROJ2,IFSBTAB2,ISSTAB,MAPORB,DET1,DET2,NTDM2,TDM2,ISTATE,JSTATE, &
                  OrbTab)

#ifdef _DMRG_
use rasscf_global, only: doDMRG
use qcmaquis_info
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: MUL
use Constants, only: Zero, One, Two, Half
#if defined (_DEBUGPRINT_) || defined (_DMRG_)
use Definitions, only: u6
#endif

implicit none
integer LSYM1, MPLET1, MSPROJ1
integer IFSBTAB1(*)
integer LSYM2, MPLET2, MSPROJ2
integer IFSBTAB2(*)
integer ISSTAB(*), MAPORB(*), NTDM2
real*8 DET1(*), DET2(*)
real*8 TDM2(NTDM2)
integer ISTATE, JSTATE
integer OrbTab(*)
integer IORB, ITABS, IUABS, JORB
integer NASHT, NASORB
real*8 SGNJL, SGNIK
real*8 GVAL, GAAAA, GABBA, GBAAB, GBBBB, GABAB, GBABA
integer ISYOP, MS2OP
integer IAAAA, IABAB, IABBA, IAKA, IAKB, IBAAB, IBABA, IBBBB, IBIA
integer IBKA, IBKB, IJ, IJIJ, IORBA, IORBB, ITU, ITUVX
integer IVABS, IVX, IXABS, JALA, JALB, JBJA, JBLA, JBLB
integer JORBA, JORBB, KORB, KORBA, KORBB, LORB, LORBA, LORBB
integer NASGEM, NSPD2
real*8, allocatable :: SPD2(:)

! Given two CI expansions, using a biorthonormal set of SD''s,
! calculate the spin-summed 2-particle transition density matrix
! in the biorthonormal active orbital basis.

! Pick out nr of active orbitals from orbital table:
NASORB = ORBTAB(4)
NASHT = NASORB/2
NASGEM = (NASORB*(NASORB-1))/2
NSPD2 = NASGEM**2
call mma_allocate(SPD2,nSPD2,Label='SPD2')
SPD2(:) = Zero
ISYOP = MUL(LSYM1,LSYM2)
MS2OP = MSPROJ1-MSPROJ2

#ifdef _DMRG_
if (.not. doDMRG) then
#endif
  call SPIND2(ISYOP,MS2OP,ORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,DET1,DET2,SPD2)
#ifdef _DMRG_
else

  write(u6,*) '2-TDM import with QCMaquis in MPSSI not implemented yet'
  call Quit_OnUserError()
end if
! Old interface import
!#define BLUBB
!if (doMPSSICheckpoints) then
!  call dmrg_interface_ctl( &
!    task='imp rdmY', &
!#   ifndef BLUBB
!    x2=spd2, &
!    mdim=nasgem, &
!#   else
!    x2=tdm2, &
!    mdim=ntdm2, &
!#   endif
!    checkpoint1=qcm_group_names(job1)%states(ist), &
!    checkpoint2=qcm_group_names(job2)%states(jst), &
!    msproj=msproj1, &
!    msprojL=msproj2, &
!    multiplet=MPLET1-1, & ! (we need 2*S)
!    multipletL=MPLET2-1, & ! (we need 2*S)
!    rdm1=.false., &
!    rdm2=.true. &
!  )
!else
!  call dmrg_interface_ctl( &
!    task='imp rdmY', &
!#   ifndef BLUBB
!    x2=spd2, &
!    mdim=nasgem, &
!#   else
!    x2=tdm2, &
!    mdim=ntdm2, &
!#   endif
!    state=iWork(lLROOT+ISTATE-1), &
!    stateL=iWork(lLROOT+JSTATE-1), &
!    msproj=msproj1, &
!    msprojL=msproj2, &
!    multiplet=MPLET1-1, & ! (we need 2*S)
!    multipletL=MPLET2-1, & ! (we need 2*S)
!    rdm1=.false., &
!    rdm2=.true. &
!  )
!end if
!#ifdef BLUBB
!goto 124
!#endif
#endif

SGNJL = One ! dummy initialize
SGNIK = One ! dummy initialize
IAKA = 0    ! dummy initialize
IAKB = 0    ! dummy initialize
IBIA = 0    ! dummy initialize
IBKA = 0    ! dummy initialize
IBKB = 0    ! dummy initialize
JALA = 0    ! dummy initialize
JALB = 0    ! dummy initialize
JBLA = 0    ! dummy initialize
JBLB = 0    ! dummy initialize
JBJA = 0    ! dummy initialize
do JORB=1,NASHT
  JORBA = 2*JORB-1
  JORBB = 2*JORB
  IUABS = MAPORB(JORBA)
  do IORB=1,NASHT
    IORBA = 2*IORB-1
    IORBB = 2*IORB
    ITABS = MAPORB(IORBA)
    ITU = ITABS+NASHT*(IUABS-1)
    do LORB=1,NASHT
      LORBA = 2*LORB-1
      LORBB = 2*LORB
      IXABS = MAPORB(LORBA)
      if (JORB > LORB) then
        SGNJL = One
        JALA = ((JORBA-1)*(JORBA-2))/2+LORBA
        JALB = ((JORBA-1)*(JORBA-2))/2+LORBB
        JBLA = ((JORBB-1)*(JORBB-2))/2+LORBA
        JBLB = ((JORBB-1)*(JORBB-2))/2+LORBB
      else if (JORB == LORB) then
        JBJA = ((JORBB-1)*(JORBB-2))/2+JORBA
      else
        SGNJL = -One
        JALA = ((LORBA-1)*(LORBA-2))/2+JORBA
        JALB = ((LORBB-1)*(LORBB-2))/2+JORBA
        JBLA = ((LORBA-1)*(LORBA-2))/2+JORBB
        JBLB = ((LORBB-1)*(LORBB-2))/2+JORBB
      end if
      do KORB=1,NASHT
        KORBA = 2*KORB-1
        KORBB = 2*KORB
        IVABS = MAPORB(KORBA)
        IVX = IVABS+NASHT*(IXABS-1)
        if (ITU < IVX) cycle
        if (IORB > KORB) then
          SGNIK = One
          IAKA = ((IORBA-1)*(IORBA-2))/2+KORBA
          IAKB = ((IORBA-1)*(IORBA-2))/2+KORBB
          IBKA = ((IORBB-1)*(IORBB-2))/2+KORBA
          IBKB = ((IORBB-1)*(IORBB-2))/2+KORBB
        else if (IORB == KORB) then
          IBIA = ((IORBB-1)*(IORBB-2))/2+IORBA
        else
          SGNIK = -One
          IAKA = ((KORBA-1)*(KORBA-2))/2+IORBA
          IAKB = ((KORBB-1)*(KORBB-2))/2+IORBA
          IBKA = ((KORBA-1)*(KORBA-2))/2+IORBB
          IBKB = ((KORBB-1)*(KORBB-2))/2+IORBB
        end if
        if (IORB /= KORB) then
          if (JORB /= LORB) then
            IAAAA = IAKA+NASGEM*(JALA-1)
            IABBA = IAKB+NASGEM*(JALB-1)
            IBAAB = IBKA+NASGEM*(JBLA-1)
            IBBBB = IBKB+NASGEM*(JBLB-1)
            GAAAA = SPD2(IAAAA)
            GABBA = SPD2(IABBA)
            GBAAB = SPD2(IBAAB)
            GBBBB = SPD2(IBBBB)
            GVAL = SGNIK*SGNJL*(GAAAA+GABBA+GBAAB+GBBBB)
          else
            IABAB = IAKB+NASGEM*(JBJA-1)
            IBAAB = IBKA+NASGEM*(JBJA-1)
            GABAB = SPD2(IABAB)
            GBAAB = SPD2(IBAAB)
            GVAL = SGNIK*(-GABAB+GBAAB)
          end if
        else
          if (JORB /= LORB) then
            IBABA = IBIA+NASGEM*(JALB-1)
            IBAAB = IBIA+NASGEM*(JBLA-1)
            GBABA = SPD2(IBABA)
            GBAAB = SPD2(IBAAB)
            GVAL = SGNJL*(-GBABA+GBAAB)
          else
            IBAAB = IBIA+NASGEM*(JBJA-1)
            GBAAB = SPD2(IBAAB)
            GVAL = Two*GBAAB
          end if
        end if
        ! Position determined by active orbital index in external order:
        ITUVX = (ITU*(ITU-1))/2+IVX
        TDM2(ITUVX) = GVAL
      end do
    end do
  end do
end do

!#ifdef BLUBB
!124  CONTINUE
!#endif

#ifdef _DEBUGPRINT_
write(u6,*) ' final 2-TDM'
do IJ=1,ntdm2
  write(u6,*) ' IJ, value = ',IJ,TDM2(IJ)
end do
#endif

call mma_deallocate(SPD2)
! DIAGONAL ELEMENTS HALF-SIZED (This is for proper contraction with TUVX):
IJIJ = 0
do IJ=1,NASHT**2
  IJIJ = IJIJ+IJ
  TDM2(IJIJ) = Half*TDM2(IJIJ)
end do

return

unused_var(MPLET1)
unused_var(MPLET2)
unused_var(ISTATE)
unused_var(JSTATE)

end subroutine MKTDM2
