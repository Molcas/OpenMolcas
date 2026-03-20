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

subroutine MKTDM1(LSYM1,MPLET1,MSPROJ1,IFSBTAB1,LSYM2,MPLET2,MSPROJ2,IFSBTAB2,ISSTAB,MAPORB,DET1,DET2,SIJ,NASHT,TDM1,TSDM1,WTDM1, &
                  ISTATE,JSTATE,job1,job2,ist,jst,OrbTab)
! Given two CI expansions, using a biorthonormal set of SD's,
! calculate the following quantities:
! (1) The overlap
! (2) The spin-summed 1-particle transition density matrix
! (3) The WE-reduced transition spin density matrix
! in the biorthonormal active orbital basis.

#ifdef _DMRG_
use rasscf_global, only: doDMRG
use qcmaquis_interface_cfg
use qcmaquis_info
use qcmaquis_interface_mpssi
use rassi_global_arrays, only: LROOT
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: MUL
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, u6

implicit none
integer LSYM1, MPLET1, MSPROJ1, LSYM2, MPLET2, MSPROJ2
#ifdef _DMRG_
real*8, allocatable :: TDMAA(:), TDMBB(:)
#endif
integer IFSBTAB1(*), IFSBTAB2(*), ISSTAB(*), MAPORB(*)
integer IORB, ISORB, ISYOP, ITABS, IUABS, JORB, JSORB
integer MS2OP, NASHT, NASORB, NSPD1
integer, intent(in) :: ISTATE, JSTATE, job1, job2, ist, jst
integer :: OrbTab(*)
real*8 DET1(*), DET2(*)
real*8 SIJ, TDM1(NASHT,NASHT), TSDM1(NASHT,NASHT), WTDM1(NASHT,NASHT)
real*8 S1, S2, SM, SM1, SM2, GAA, GAB, GBA, GBB
real*8 OVERLAP_RASSI, TMATEL, RED, FACT, CGCOEF, DCLEBS
real*8, allocatable :: SPD1(:)

! Pick out nr of active orbitals from orbital table:
NASORB = ORBTAB(4)

! Overlap:

SIJ = Zero

if ((MPLET1 == MPLET2) .and. (MSPROJ1 == MSPROJ2)) then

# ifdef _DMRG_
  if (.not. doDMRG) then
# endif

    SIJ = OVERLAP_RASSI(IFSBTAB1,IFSBTAB2,DET1,DET2)

# ifdef _DMRG_
  else
    sij = qcmaquis_mpssi_overlap(qcm_prefixes(job1),istate,qcm_prefixes(job2),jstate,.true.)

  end if ! DMRG or not

# endif
end if ! mmplet and msproj check

! General 1-particle transition density matrix:
NSPD1 = NASORB**2
#ifdef _DMRG_
if (.not. doDMRG) then
#endif
  call mma_allocate(SPD1,nSPD1,Label='SPD1')
  SPD1(:) = Zero
#ifdef _DMRG_
else
  ! For DMRG, we only need the AA and BB spin components
  ! Let's allocate two different arrays for them because that's
  ! easier for the new QCMaquis interface
  NSPD1 = NASHT**2
  call mma_allocate(TDMAA,nSPD1,Label='TDMAA')
  call mma_allocate(TDMBB,nSPD1,Label='TDMBB')
  TDMAA(:) = Zero
  TDMBB(:) = Zero
end if
#endif
ISYOP = MUL(LSYM1,LSYM2)
MS2OP = MSPROJ1-MSPROJ2

if (abs(MS2OP) <= 2) then
# ifdef _DMRG_
  if (.not. doDMRG) then
# endif
    !> spind constructs the 1-particle transition density matrix
    !> output in SPD1
    !> main input: DET1 and DET2
    call SPIND(ISYOP,MS2OP,ORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,DET1,DET2,SPD1)

# ifdef _DMRG_
  else
    if (isyop /= 1) then
      write(u6,*) 'MPS property density with spatial symm irrep > 1: FIXME!'
      call abend()
    end if
    ! calculate 1-TDMs: Must always be calculated with the higher multiplicity as <T|o|S>
    ! where T always has a higher multiplicity than S

    if (MPLET1 < MPLET2) then
      call qcmaquis_mpssi_get_onetdm_spin(qcm_prefixes(job2),LROOT(JSTATE),qcm_prefixes(job1),LROOT(ISTATE),TDMAA,TDMBB,NSPD1)
    else
      call qcmaquis_mpssi_get_onetdm_spin(qcm_prefixes(job1),LROOT(ISTATE),qcm_prefixes(job2),LROOT(JSTATE),TDMAA,TDMBB,NSPD1)
    end if
  end if
# endif
end if
! Create a scalar, and an WE-reduced spin, transition density matrix.
! The scalar matrix is simply the usual spin-summed density matrix.
! The WE-reduced matrix is the one defined through Wigner-Eckarts thm:
!  <A S1 M1 | T[K]_Q |B S2 M2>
!                   = FACT*CG(S2 M2 K Q;S1 M1)*<A S1 || T[K] ||B S2>
! i.e. with the usual Clebsch-Gordan factor, and a prefactor
!   FACT=(-1)**(MAX(S1,S2)-S1)/SQRT(2*S1+1)
S1 = Half*real(MPLET1-1,kind=wp)
S2 = Half*real(MPLET2-1,kind=wp)
SM1 = Half*real(MSPROJ1,kind=wp)
SM2 = Half*real(MSPROJ2,kind=wp)
do IORB=1,NASHT
  ISORB = 2*IORB-1
  do JORB=1,NASHT
    JSORB = 2*JORB-1
#   ifdef _DMRG_
    if (.not. doDMRG) then
#   endif
      GAA = SPD1(ISORB+NASORB*(JSORB-1))
      GAB = SPD1(ISORB+NASORB*(JSORB))
      GBA = SPD1(1+ISORB+NASORB*(JSORB-1))
      GBB = SPD1(1+ISORB+NASORB*(JSORB))
#   ifdef _DMRG_
    else
      GAA = TDMAA(JORB+NASHT*(IORB-1))
      GBB = TDMBB(JORB+NASHT*(IORB-1))
      ! transpose from row-major order,
      ! as it comes from C++ this way
      GAB = Zero
      GBA = Zero
    end if
#   endif

    ! Position determined by active orbital index in external order:
    ITABS = MAPORB(ISORB)
    IUABS = MAPORB(JSORB)

    !> scalar TDM
    TDM1(ITABS,IUABS) = GAA+GBB
    !> spin TDM
    TSDM1(ITABS,IUABS) = GAA-GBB

    ! Clebsch-Gordan coefficient:
    SM = SM1-SM2
    FACT = One/sqrt(real(MPLET1,kind=wp))
    if (MPLET1 == MPLET2-2) FACT = -FACT
    CGCOEF = FACT*DCLEBS(S2,One,S1,SM2,SM,SM1)
    ! Spin tensor component matrix element
    if (MSPROJ2 == MSPROJ1+2) then
      TMATEL = sqrt(Two)*GBA
    else if (MSPROJ2 == MSPROJ1-2) then
      TMATEL = -sqrt(Two)*GAB
    else if (MSPROJ2 == MSPROJ1) then
      TMATEL = Half*(GBB-GAA)
    end if
    ! Thus obtain reduced matrix element from Wigner-Eckart theorem:
    RED = Zero
    if (CGCOEF /= Zero) then
      RED = TMATEL/CGCOEF
    else
      if (abs(TMATEL) > 1.0e-12_wp) then
        call WarningMessage(1,'A possible bug was detected.')
        write(u6,*) ' WARNING: Non-zero matrix element computed'
        write(u6,*) ' which should be zero by spin symmetry!'
        write(u6,*) '              Spins S1, S2:',S1,S2
        write(u6,*) ' Spin projections SM1, SM2:',SM1,SM2
        write(u6,*) '    Operator has S=1.0, SM:',SM
        write(u6,*) ' Clebsch-Gordan:',CGCOEF
        write(u6,*) ' Size is TMATEL=',TMATEL
      end if
    end if

    !> W-reduced TDM
    WTDM1(ITABS,IUABS) = RED

  end do
end do

! Avoid unused argument warnings
unused_var(ISTATE)
unused_var(JSTATE)
unused_var(job1)
unused_var(job2)
unused_var(ist)
unused_var(jst)

#ifdef _DMRG_
if (.not. doDMRG) then
#endif
  call mma_deallocate(SPD1)

#ifdef _DMRG_
else
  call mma_deallocate(TDMAA)
  call mma_deallocate(TDMBB)
end if
#endif

end subroutine MKTDM1
