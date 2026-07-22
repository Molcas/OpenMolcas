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
! Copyright (C) 1988,1991,1992,1998, Per Ake Malmqvist                 *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine POLY3(mkF)
!  IBM TEST VERSION 0, 1988-06-23.
!  NEW VERSION 1991-02-23, FOR USE WITH RASSCF IN MOLCAS PACKAGE.
!  NEW VERSION 1992-12-05, FOR MOLCAS-3 VERSION.
!  NEW VERSION 1998-10-02
!  AUTHOR: PER-AAKE MALMQUIST

! THIS PROGRAM CALCULATES 1-EL, 2-EL, AND 3-EL
! DENSITY MATRICES FOR A CASSCF WAVE FUNCTION.
! IF THE INTEGER KEY IFF == 1, THEN
! IT ALSO PRODUCES THE CONTRACTIONS OF 1-EL -- 4-EL
! DENSITY MATRICES WITH THE FOCK OPERATOR USED IN
! THE CASSCF-MP2 PROGRAM. THE RESULTS ARE WRITTEN
! TO FILE IN SEVERAL FORMS, TO SUPPORT BOTH KERSTINS
! PRESENT PROGRAM AND ALSO SUCH NEW PROCEDURES WHICH
! MIGHT TAKE ADVANTAGE OF ALL INDEX PERMUTATION SYMMETRIES.
! THE RDSTAT AND THE GUGA ROUTINES USED IN THIS
! PROGRAM ASSUMES THE JOBIPH IS PRODUCED BY THE RASSCF PROGRAM.

use fciqmc_interface, only: DoFCIQMC
use PrintLevel, only: VERBOSE
use sguga_states, only: SGS, CIS
use caspt2_global, only: IDTCEX, iPrGlb, LUCIEX, LUSOLV
use general_data, only: NACTEL, STSym, nLev
use caspt2_module, only: CIThr, DoCumulant, EPSA, Eta, iSCF, jState, mState, NAshT, nConf, nG1, nG2, nG3, nG3Tot, nState
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_ || defined _DMRG_
use caspt2_module, only: DMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, byte

implicit none
logical(kind=iwp), intent(in) :: mkF

integer(kind=iwp) :: IDCI, ILEV, ILUID, IPARDIV, nCI, NG3MAX
integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: CI(:)
real(kind=wp), allocatable, target :: F1_H(:), F2_H(:), F3_H(:), G1(:), G2(:), G3(:)
real(kind=wp), pointer :: F1(:), F2(:), F3(:)
integer(kind=iwp), parameter :: istate=1
nLev = SGS(istate)%nLev

! Note that in case of FCIQMC nConf is set to 0.
nCI = CIS(istate)%NCSF(STSYM)

if (mkF) then
  ! ORBITAL ENERGIES IN CI-COUPLING ORDER:
  do ILEV=1,NLEV
    ETA(ILEV) = EPSA(SGS(istate)%L2ACT(ILEV))
  end do
end if

call mma_allocate(G1,NG1,LABEL='G1')
call mma_allocate(G2,NG2,LABEL='G2')

!-SVC20100831: recompute approximate max NG3 size needed
NG3MAX = iPARDIV(NG3TOT,NG2)

!-SVC20100831: allocate local G3 matrices
call mma_allocate(G3,NG3MAX,LABEL='G3')

call mma_allocate(idxG3,6,NG3MAX,label='idxG3')
idxG3(:,:) = 0

G1(1) = Zero
G2(1) = Zero
G3(1) = Zero

! ALLOCATE SPACE FOR CORRESPONDING COMBINATIONS WITH H0:
if (mkF) then
  call mma_allocate(F1_H,NG1,LABEL='F1_H')
  call mma_allocate(F2_H,NG2,LABEL='F2_H')
  call mma_allocate(F3_H,NG3MAX,LABEL='F3_H')
  F1 => F1_H
  F2 => F2_H
  F3 => F3_H
else
  ! This is just done such that in the case of mkF=.false. that
  ! F1, F2, and F3 refer to an actual array.
  F1 => G1
  F2 => G2
  F3 => G3
end if

! NG3 will change inside subroutine MKFG3 to the actual
! number of nonzero elements, that is why here we allocate
! with NG3MAX, but we only store (PT2_PUT) the first NG3
! elements of the G3 and F3
NG3 = NG3MAX

if (.not. DoFCIQMC) then
  call mma_allocate(CI,NCONF,Label='CI')

  if ((.not. DoCumulant) .and. (ISCF == 0)) then
    IDCI = IDTCEX(JSTATE)
    call DDAFILE(LUCIEX,2,CI,NCONF,IDCI)
    if (IPRGLB >= VERBOSE) then
      write(u6,*)
      if (NSTATE > 1) then
        write(u6,'(A,I4)') ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
      else
        write(u6,*) ' With new orbitals, the CI array is:'
      end if
      call PRWF_CP2(STSYM,NCONF,CI,CITHR)
    end if
  else
    CI(1) = One
  end if
end if

if ((ISCF /= 0) .and. (NACTEL /= 0)) then
  call SPECIAL(G1,G2,G3,F1,F2,F3,idxG3,nAshT,nG3)
else if (ISCF == 0) then
  !-SVC20100903: during mkfg3, NG3 is set to the actual value
# if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_ || defined _DMRG_
  if ((.not. DoCumulant) .and. (.not. DMRG)) then
# endif
    if (.not. allocated(CI)) call mma_allocate(CI,1,LABEL='CI')
    call MKFG3(mkF,CI,nCI,G1,F1,G2,F2,G3,F3,idxG3,nLev,nG1,nG2,nG3)
# if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_ || defined _DMRG_
  else
    call MKFG3DM(mkF,G1,F1,G2,F2,G3,F3,idxG3,nLev,nG3)
  end if
# endif
end if

if (allocated(CI)) call mma_deallocate(CI)

if (NLEV > 0) then
  call PT2_PUT(NG1,' GAMMA1',G1)
  call PT2_PUT(NG2,' GAMMA2',G2)
  call PT2_PUT(NG3,' GAMMA3',G3)
  iLUID = 0
  call I1DAFILE(LUSOLV,1,idxG3,6*NG3,iLUID)
  if (mkF) then
    call PT2_PUT(NG1,' DELTA1',F1)
    call PT2_PUT(NG2,' DELTA2',F2)
    call PT2_PUT(NG3,' DELTA3',F3)
  end if
end if

if (NLEV > 0) then
  call mma_deallocate(G1)
  call mma_deallocate(G2)
  call mma_deallocate(G3)
  call mma_deallocate(idxG3)
  if (mkF) then
    call mma_deallocate(F1_H)
    call mma_deallocate(F2_H)
    call mma_deallocate(F3_H)
  end if
  nullify(F1,F2,F3)
end if

end subroutine POLY3
