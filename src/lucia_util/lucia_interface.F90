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

module LUCIA_INTERFACE

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
private

public :: Lucia_Util

contains

!***********************************************************************
!  Lucia_Util
!
!> @brief
!>   Wrapper for using LUCIA utils in MOLCAS.
!> @author Jesper Wisborg Krogh
!>
!> @details
!> By using the LUCIA utils through this wrapper it is guaranteed
!> that all common blocks used have a common parent routine.
!>
!> @param[in] Module Identifier
!***********************************************************************
subroutine Lucia_Util(ModLab,iSym,iDisk,LU,Array,RVec,CI_VECTOR,SIGMA_VECTOR)

  use lucia_data, only: IREFSM, MXNTTS

  implicit none
  character(len=*), intent(in) :: ModLab
  integer(kind=iwp), intent(in), optional :: iSym, LU
  integer(kind=iwp), intent(inout), optional :: iDisk
  real(kind=wp), intent(in), optional :: Array(:), RVEC(:)
  real(kind=wp), intent(_IN_), optional :: CI_Vector(:)
  real(kind=wp), intent(out), optional :: SIGMA_Vector(:)
  character(len=72) :: Module_
  integer(kind=iwp), allocatable :: lVec(:)
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: COUNTER = 0

  COUNTER = COUNTER+1
  write(u6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)') '[',COUNTER,']','ENTRY LUCIA_UTIL','(',ModLab,')'
# endif

  ! Make sure the ModLab variable is in upper case.

  Module_ = ModLab
  call UpCase(Module_)

  ! Call the appropriate routines according to ModLab

  if (Module_(1:4) == 'DIAG') then
    call Diag_Master()
  else if (Module_(1:9) == 'SIGMA_CVB') then
    ! iSym_LI is the symmetry to be used.
    call Sigma_Master_CVB(CI_VECTOR,SIGMA_VECTOR,iSym)
  else if (Module_(1:5) == 'SIGMA') then
    !write(u6,*) 'blubbbbbbhc'
    call Sigma_Master(CI_VECTOR,SIGMA_VECTOR)
  else if (Module_(1:5) == 'TRACI') then
    !write(u6,*) 'blubbbbbbtraci'
    ! iDisk is the initial disk address (for read/write of JOBIPH)
    ! Lu is the file unit for JOBIPH
    ! Array is the transformation matrix (not sorted as LUCIA needs it).
    call mma_allocate(lVec,MXNTTS,Label='lVec')
    call Traci_Master(iDisk,LU,Array,lVec)
    call mma_deallocate(lVec)
  else if (Module_(1:5) == 'DENSI') then
    if (present(RVEC)) then
      call Densi_Master(CI_VECTOR,RVEC=RVEC(:))
    else
      call Densi_Master(CI_VECTOR)
    end if
  else if (Module_(1:3) == 'INI') then
    call Lucia_Ini()
    call DetCtl_Gas()
  else if (Module_(1:5) == 'CLOSE') then
    call CSFDIM_FREE(IREFSM)
    call LUCIA2MOLCAS_FREE()
    call Lucia_Close()
  else
    write(u6,*) 'Unknown module requested in Lucia_Util.'
    write(u6,*) 'Module = ',ModLab
    write(u6,*) 'Known modules are:'
    write(u6,*) 'Diag, Sigma, Sigma_CVB, Densi, DetCtl, Ini'
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  write(u6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)') '[',COUNTER,']','EXIT LUCIA_UTIL','(',ModLab,')'
# endif

end subroutine Lucia_Util

subroutine densi_master(CIVec,RVec)
  ! Controls the calculation of the densities, when Lucia is called
  ! from Molcas Rasscf.

  use lucia_data, only: DSTmp, Dtmp, DTOC, IDISK, IREFSM, LCSBLK, kvec3_length, LUC, LUHC, LUSC1, LUSC34, MXNTTS, MXSOOB, &
                        NCSF_PER_SYM, NSD_PER_SYM, NTOOB, PAtmp, PSSIGN, Ptmp, RHO1, SDREO, Sigma_on_Disk, SRHO1, VEC3, XISPSM
  use Constants, only: Zero

  implicit none
  real(kind=wp), intent(in) :: CIVec(:)
  real(kind=wp), intent(in), optional :: RVec(:)
  integer(kind=iwp) :: LBLK, LBLOCK, NCSF, NSD
  real(kind=wp) :: dummy(1), EXPS2
  logical(kind=iwp) :: iPack, tdm
  integer(kind=iwp), allocatable :: lVec(:)
  real(kind=wp), allocatable :: SCR2(:), SCR4(:), VEC1(:), VEC2(:)
  real(kind=wp), allocatable, target :: SCR1(:), SCR3(:)

  ! Put CI-vector from RASSCF on luc

  ! if rvec is associated, it should be a pointer to a second CI vector
  ! and a one-particle transition density matrix will be computed
  tdm = present(rVec)

  NSD = NSD_PER_SYM(IREFSM)
  NCSF = NCSF_PER_SYM(IREFSM)
  call mma_allocate(SCR1,NSD,Label='SCR1')
  call mma_allocate(SCR2,NSD,Label='SCR2')

  SCR1(1:NCSF) = CIVEC(1:NCSF)

  call mma_allocate(lVec,MXNTTS,Label='lVec')
  if (tdm) then
    call mma_allocate(SCR3,NSD,Label='SCR3')
    call mma_allocate(SCR4,NSD,Label='SCR4')

    SCR3(1:NCSF) = rvec(1:NCSF)
    call CSDTVC(SCR3,SCR4,1,DTOC,SDREO,IREFSM,1)
    call CPCIVC1(SCR3,NSD,LUHC,MXNTTS,IREFSM,lVec)
    call mma_deallocate(SCR3)
    call mma_deallocate(SCR4)
  end if
  call CSDTVC(SCR1,SCR2,1,DTOC,SDREO,IREFSM,1)
  call CPCIVC1(SCR1,NSD,LUC,MXNTTS,IREFSM,lVec)
  call mma_deallocate(lVec)

  ! Determine length of arrays VEC1 and VEC2

  !if (ISIMSYM == 0) then
  LBLOCK = MXSOOB
  !else
  !  LBLOCK = MXSOOB_AS
  !end if
  LBLOCK = max(LBLOCK,LCSBLK)
  ! JESPER : Should reduce I/O
  !PAM06 LBLOCK = max(XISPSM(IREFSM,1),real(MXSOOB,kind=wp))
  LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
  if (PSSIGN /= Zero) LBLOCK = 2*int(XISPSM(IREFSM,1))

  ! Allocate arrays

  call mma_allocate(VEC1,LBLOCK,Label='VEC1')
  call mma_allocate(VEC3,kvec3_length,Label='VEC3')

  ! Copy Sigma-vector from disc to core

  call mma_allocate(VEC2,LBLOCK,Label='VEC2')
  if (Sigma_on_disk) then
    call mma_allocate(lVec,MXNTTS,Label='lVec')
    call cpsivc(lusc34,mxntts,vec2,lVec)
    call mma_deallocate(lVec)
  else
    vec2(:) = Zero
  end if

  ! Information needed on file handling

  LBLK = -1

  ! Copy vector on file LUC to LUSC1 and LUHC

  IDISK(LUC) = 0
  IDISK(LUSC1) = 0
  call COPVCD(LUC,LUSC1,VEC1,0,LBLK)
  if (.not. tdm) call COPVCD(LUSC1,LUHC,VEC1,1,LBLK)

  ! Calculate one- and two-body densities

  IPACK = .true.
  DUMMY = Zero
  if (tdm) then
    call densi2(1,Dtmp,dummy,dummy,dummy,vec1,vec2,lusc1,luhc,exps2,1,DStmp,IPACK)
  else
    call densi2(2,rho1,dummy,Ptmp,PAtmp,vec1,vec2,lusc1,luhc,exps2,1,srho1,IPACK)
  end if

  ! Explanation of calling parameters
  !
  ! 2     : DONE!!! - Calculate both one and two body densities.
  ! rho1  : DONE!!! - Output - include in module lucia_data
  ! rho2  : DONE!!! - Output - include in module lucia_data
  ! vec1  : DONE!!! - CI-vector
  ! vec2  : DONE!!! - Sigma-vector
  ! lusc1 : DONE!!! - file pointer
  ! luhc  : DONE!!! - file pointer
  ! exps2 : DONE!!! - Output - expectation value of S**2.
  ! 1     : DONE!!! - Calculate spin density
  ! srho1 : DONE!!! - Comming with module lucia_data

  if (.not. tdm) then
    ! Save densities in trigonal format for use in Molcas

    call TriPak(rho1,Dtmp,ntoob,ntoob)
    call TriPak(srho1,DStmp,ntoob,ntoob)
  end if

  call CSDTVC(scr1,scr2,2,dtoc,SDREO,iRefSm,1)

  call mma_deallocate(SCR1)
  call mma_deallocate(SCR2)
  call mma_deallocate(VEC1)
  call mma_deallocate(VEC2)
  call mma_deallocate(VEC3)

end subroutine densi_master

subroutine sigma_master(CIVEC,SIGMAVEC)
  ! Controls the calculation of the sigma vector, when Lucia is called
  ! from Molcas Rasscf.

  use lucia_data, only: CI_VEC, ECORE, ECORE_ORIG, INI_H0, INT1, INT1O, IREFSM, KVEC3_LENGTH, LUC, LUSC34, MXNTTS, NSD_PER_SYM, &
                        SIGMA_VEC, VEC3

  implicit none
  real(kind=wp), intent(_IN_) :: CIVEC(:)
  real(kind=wp), intent(out) :: SIGMAVEC(:)
  integer(kind=iwp) :: nSD
  integer(kind=iwp), allocatable :: lVec(:)

  nSD = NSD_PER_SYM(IREFSM)

  ! Put CI-vector from RASSCF on luc and get h0 from Molcas enviroment.

  if (INI_H0 == 0) ECORE = ECORE_ORIG
  INI_H0 = 0
  INT1(:) = INT1O(:)
  ECORE_ORIG = ECORE
  !if (IUSE_PH == 1) then
  !  call FI(INT1,ECORE_HEX,1)
  !  ECORE = ECORE+ECORE_HEX
  !end if
  call mma_allocate(lVec,MXNTTS,Label='lVec')
  call CPCIVC1(CIVEC,nSD,LUC,MXNTTS,IREFSM,lVec)
  call mma_deallocate(lVec)

  ! Calculate the sigma vector:
  ! LUC       : unit from which the CI vector is picked
  ! LUSC34    : unit to which the sigma vector is put
  ! CI_VEC    : scratch area to store the CI vector temporarily
  ! SIGMA_VEC : the computed sigmavector in core.

  call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
  ! Note that CI_VEC is used as a scratch array!
  call MV7(CI_VEC,SIGMA_VEC,LUC,LUSC34)
  call mma_deallocate(VEC3)

  ! Export lusc34 to RASSCF

  call mma_allocate(lVec,MXNTTS,Label='lVec')
  call CPCIVC2(SIGMAVEC,nSD,LUSC34,MXNTTS,IREFSM,lVec)
  call mma_deallocate(lVec)

end subroutine SIGMA_MASTER

subroutine SIGMA_MASTER_CVB(CIVEC,SIGMAVEC,IREFSM_CASVB)

  use CandS, only: ICSM, ISSM
  use lucia_data, only: CI_VEC, ECORE, ECORE_ORIG, INI_H0, INT1, INT1O, IREFSM, KVEC3_LENGTH, LUC, LUSC34, MXNTTS, NSD_PER_SYM, &
                        SIGMA_ON_DISK, VEC3

  implicit none
  integer(kind=iwp), intent(in) :: IREFSM_CASVB
  real(kind=wp), intent(_IN_) :: CIVEC(:)
  real(kind=wp), intent(out) :: SIGMAVEC(:)
  integer(kind=iwp) :: nSD
  integer(kind=iwp), allocatable :: lVec(:)

  ! Set ICSM and ISSM (from module CandS to the correct symmetry for this call

  ICSM = IREFSM_CASVB
  ISSM = IREFSM_CASVB

  !nSD = NSD_PER_SYM(IREFSM)
  nSD = NSD_PER_SYM(IREFSM_CASVB)

  ! Get h0 from Molcas enviroment.

  if (INI_H0 == 0) ECORE = ECORE_ORIG
  INI_H0 = 0
  INT1(:) = INT1O(:)
  ECORE_ORIG = ECORE
  !if (IUSE_PH == 1) then
  !  call FI(INT1,ECORE_HEX,1)
  !  ECORE = ECORE+ECORE_HEX
  !end if

  ! Write CI-vector to disc

  call mma_allocate(lVec,MXNTTS,Label='lVec')
  call CPCIVC1(CIVEC,nSD,LUC,MXNTTS,ISSM,lVec)
  call mma_deallocate(lVec)

  ! Calculate the sigma vector:

  call DIAG_MASTER()
  call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
  ! Note that CI_VEC is used as a scratch array!
  call MV7(CI_VEC,SIGMAVec,LUC,LUSC34)
  call mma_deallocate(VEC3)

  ! Export lusc34 to RASSCF

  SIGMA_ON_DISK = .true.

  ! Set ICSM and ISSM (from CandS) back to IREFSM

  ICSM = IREFSM
  ISSM = IREFSM

end subroutine SIGMA_MASTER_CVB

subroutine cpcivc1(CIVec,nCIVEC,ifile,mxrec,isym,lrec)
  ! Copies the CI-vector from Molcas to Lucia (from core to disk unit ifile).

  use lucia_data, only: IDISK

  implicit none
  integer(kind=iwp), intent(in) :: nCIVEC, ifile, mxrec, isym
  real(kind=wp), intent(_IN_) :: CIVec(nCIVec)
  integer(kind=iwp), intent(out) :: lrec(mxrec)
  integer(kind=iwp) :: dum(1), nRec
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: iOff, iRec
# endif

  ! ==================
  ! Find nrec and lrec
  ! ==================

  call blkfo_min(isym,nrec,lrec)
  IDISK(IFILE) = 0

  ! =======================
  ! Write CI-vector to disc
  ! =======================

# ifdef _DEBUGPRINT_
  ioff = 1
  write(u6,*) 'CI-vector put to disk:'
  do IREC=1,NREC
    if (LREC(IREC) >= 0) then
      call wrtmat(CIVec(ioff),1,lrec(irec),1,lrec(irec))
      ioff = ioff+lrec(irec)
    end if
  end do
# endif
  call todscn(CIVec,nrec,lrec,-1,ifile)
  dum(1) = -1
  call itods(dum,1,-1,ifile)

end subroutine cpcivc1

subroutine cpcivc2(CIVec,nCIVEC,ifile,mxrec,isym,lrec)
  ! Copies the CI-vector from Lucia to Molcas (from disk unit ifile to core).

  use lucia_data, only: IDISK

  implicit none
  integer(kind=iwp), intent(in) :: nCIVEC, ifile, mxrec, isym
  real(kind=wp), intent(out) :: CIVec(nCIVec)
  integer(kind=iwp), intent(out) :: lrec(mxrec)
  integer(kind=iwp) :: nRec

  ! ==================
  ! Find nrec and lrec
  ! ==================

  call blkfo_min(isym,nrec,lrec)
  IDISK(IFILE) = 0

  ! ========================
  ! Read CI-vector from disc
  ! ========================
  call frmdscn(CIVec,nrec,-1,ifile)

end subroutine cpcivc2

subroutine cpsivc(ifile,mxrec,vec,lrec)
  ! Copies the Sigma-vector between Molcas Rasscf and Lucia enviroment

  use lucia_data, only: IDISK
  use general_data, only: STSYM

  implicit none
  integer(kind=iwp), intent(in) :: ifile, mxrec
  real(kind=wp), intent(inout) :: vec(mxrec)
  integer(kind=iwp), intent(out) :: lrec(mxrec)
  integer(kind=iwp) :: nrec

  ! ==================
  ! Find nrec and lrec
  ! ==================

  call blkfo_min(stsym,nrec,lrec)
  IDISK(IFILE) = 0

  ! ===========================
  ! Read Sigma-vector from disc
  ! ===========================

  call frmdscn(vec,nrec,-1,ifile)

end subroutine cpsivc

end module LUCIA_INTERFACE
