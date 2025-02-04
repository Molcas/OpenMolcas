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

private

public Lucia_Util

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
subroutine Lucia_Util(module,iSym,iDisk,LU,Array,RVec,CI_VECTOR,SIGMA_VECTOR)

  use stdalloc, only: mma_allocate, mma_deallocate
  use lucia_data, only: MXNTTS
  use Definitions, only: u6

  implicit none
  character(len=*) module
  integer, optional :: iSym
  integer, optional :: iDisk
  integer, optional :: LU
  real*8, optional :: Array(:)
  real*8, optional :: RVEC(:)
  real*8, optional :: CI_Vector(:)
  real*8, optional :: SIGMA_Vector(:)
  integer, parameter :: MxpLnc = 72
  character(len=MxpLnc) Module_
  integer, allocatable :: lVec(:)
# ifdef _DEBUGPRINT_
  integer, save :: COUNTER = 0

  COUNTER = COUNTER+1
  write(u6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)') '[',COUNTER,']','ENTRY LUCIA_UTIL','(',module,')'
# endif

  ! Make sure the Module variable is in upper case.

  Module_ = module
  call UppCas(Module_,MxpLnc)

  ! Call the appropriate routines according to Module

  if (Module_(1:4) == 'DIAG') then
    call Diag_Master()
  else if (Module_(1:9) == 'SIGMA_CVB') then
    ! iSym_LI is the symmetry to be used.
    call Sigma_Master_CVB(CI_VECTOR,SIGMA_VECTOR,size(CI_VECTOR),iSym)
  else if (Module_(1:5) == 'SIGMA') then
    !write(u6,*) 'blubbbbbbhc'
    call Sigma_Master(CI_VECTOR,SIGMA_VECTOR,size(CI_VECTOR))
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
      call Densi_Master(CI_VECTOR,size(CI_VECTOR),RVEC=RVEC(:))
    else
      call Densi_Master(CI_VECTOR,size(CI_VECTOR))
    end if
  else if (Module_(1:3) == 'INI') then
    call Lucia_Ini()
    call DetCtl_Gas()
  else if (Module_(1:5) == 'CLOSE') then
    call DetCtl_Free()
    call Lucia_Close()
  else
    write(u6,*) 'Unknown module requested in Lucia_Util.'
    write(u6,*) 'Module = ',module
    write(u6,*) 'Known modules are:'
    write(u6,*) 'Diag, Sigma, Sigma_CVB, Densi, DetCtl, Ini'
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  write(u6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)') '[',COUNTER,']','EXIT LUCIA_UTIL','(',module,')'
# endif

end subroutine Lucia_Util

subroutine densi_master(CIVec,nCIVec,RVec)

  use stdalloc, only: mma_allocate, mma_deallocate
  use GLBBAS, only: VEC3, DTOC, RHO1, SRHO1, SDREO
  use rasscf_lucia, only: kvec3_length, Sigma_on_Disk, PAtmp, Ptmp, DSTmp, Dtmp
  use lucia_data, only: NCSF_PER_SYM, NSD_PER_SYM
  use lucia_data, only: MXSOOB, MXNTTS, XISPSM
  use lucia_data, only: LUC, LUSC1, LUHC, LUSC34
  use lucia_data, only: LCSBLK
  use lucia_data, only: IREFSM, PSSIGN
  use lucia_data, only: IDISK
  use lucia_data, only: NTOOB
  use Constants, only: Zero

  ! Controls the calculation of the densities, when Lucia is called
  ! from Molcas Rasscf.

  implicit none
  integer nCIVec
  real*8 CIVec(nCIVEC)
  real*8, optional :: RVec(:)
  logical iPack, tdm
  real*8 dummy(1)
  real*8, allocatable :: VEC1(:), VEC2(:)
  integer, allocatable :: lVec(:)
  real*8, allocatable, target :: SCR1(:), SCR3(:)
  real*8, allocatable :: SCR2(:), SCR4(:)
  integer NSD, NCSF, LBLOCK, LBLK
  real*8 EXPS2

  ! Put CI-vector from RASSCF on luc

  ! if rvec is associated, it should be a pointer to a second CI vector
  ! and a one-particle transition density matrix will be computed
  tdm = present(rVec)

  NSD = NSD_PER_SYM(IREFSM)
  NCSF = NCSF_PER_SYM(IREFSM)
  call mma_allocate(SCR1,NSD,Label='SCR1')
  call mma_allocate(SCR2,NSD,Label='SCR2')

  call COPVEC(CIVEC,SCR1,NCSF)

  call mma_allocate(lVec,MXNTTS,Label='lVec')
  if (tdm) then
    call mma_allocate(SCR3,NSD,Label='SCR3')
    call mma_allocate(SCR4,NSD,Label='SCR4')

    call COPVEC(rvec,SCR3,NCSF)
    call CSDTVC(SCR3,SCR4,1,DTOC,SDREO,IREFSM,1)
    call CPCIVC(SCR3,NSD,LUHC,MXNTTS,IREFSM,1,lVec)
    call mma_deallocate(SCR3)
    call mma_deallocate(SCR4)
  end if
  call CSDTVC(SCR1,SCR2,1,DTOC,SDREO,IREFSM,1)
  call CPCIVC(SCR1,NSD,LUC,MXNTTS,IREFSM,1,lVec)
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
    call densi2_lucia(1,Dtmp,dummy,dummy,dummy,vec1,vec2,lusc1,luhc,exps2,1,DStmp,IPACK)
  else
    call densi2_lucia(2,rho1,dummy,Ptmp,PAtmp,vec1,vec2,lusc1,luhc,exps2,1,srho1,IPACK)
  end if

  ! Explanation of calling parameters
  !
  ! 2     : DONE!!! - Calculate both one and two body densities.
  ! rho1  : DONE!!! - Output - include in module glbbas
  ! rho2  : DONE!!! - Output - include in module glbbas
  ! vec1  : DONE!!! - CI-vector
  ! vec2  : DONE!!! - Sigma-vector
  ! lusc1 : DONE!!! - file pointer
  ! luhc  : DONE!!! - file pointer
  ! exps2 : DONE!!! - Output - expectation value of S**2.
  ! 1     : DONE!!! - Calculate spin density
  ! srho1 : DONE!!! - Comming with module glbbas

  if (.not. tdm) then
    ! Save densities in trigonal format for use in Molcas

    call TriPak(rho1,Dtmp,1,ntoob,ntoob)
    call TriPak(srho1,DStmp,1,ntoob,ntoob)
  end if

  call CSDTVC(scr1,scr2,2,dtoc,SDREO,iRefSm,1)

  call mma_deallocate(SCR1)
  call mma_deallocate(SCR2)
  call mma_deallocate(VEC1)
  call mma_deallocate(VEC2)
  call mma_deallocate(VEC3)

end subroutine densi_master

subroutine sigma_master(CIVEC,SIGMAVEC,nCIVEC)

  use stdalloc, only: mma_allocate, mma_deallocate
  ! Note that CI_VEC is used as a scratch array!
  use GLBBAS, only: INT1, INT1O, VEC3, CI_SCR => CI_VEC, SIGMA_SCR => SIGMA_VEC
  use rasscf_lucia, only: INI_H0, KVEC3_LENGTH
  use lucia_data, only: NSD_PER_SYM
  use lucia_data, only: ECORE, ECORE_ORIG
  use lucia_data, only: MXNTTS
  use lucia_data, only: LUC, LUSC34
  use lucia_data, only: IREFSM

  ! Controls the calculation of the sigma vector, when Lucia is called
  ! from Molcas Rasscf.

  implicit none
  integer nCIVEC
  real*8 CIVEC(nCIVEC), SIGMAVEC(nCIVEC)

  integer, allocatable :: lVec(:)
  integer nSD

  nSD = NSD_PER_SYM(IREFSM)

  ! Put CI-vector from RASSCF on luc and get h0 from Molcas enviroment.

  if (INI_H0 == 0) ECORE = ECORE_ORIG
  INI_H0 = 0
  INT1(:) = INT1O(:)
  ECORE_ORIG = ECORE
  !if (IUSE_PH == 1) then
  !   call FI(INT1,ECORE_HEX,1)
  !   ECORE = ECORE+ECORE_HEX
  !end if
  call mma_allocate(lVec,MXNTTS,Label='lVec')
  call CPCIVC(CIVEC,nSD,LUC,MXNTTS,IREFSM,1,lVec)
  call mma_deallocate(lVec)

  ! Calculate the sigma vector:
  ! LUC      : unit from which the CI vector is picked
  ! LUSC34   : unit to which the sigma vector is put
  ! CI_SCR   : scratch area to store the CI vector temporarily
  ! SIGMAVec : the computed sigmavector in core.

  call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
  call MV7(CI_SCR,SIGMA_SCR,LUC,LUSC34)
  call mma_deallocate(VEC3)

  ! Export lusc34 to RASSCF

  call mma_allocate(lVec,MXNTTS,Label='lVec')
  call CPCIVC(SIGMAVEC,nSD,LUSC34,MXNTTS,IREFSM,2,lVec)
  call mma_deallocate(lVec)

end subroutine SIGMA_MASTER

subroutine SIGMA_MASTER_CVB(CIVEC,SIGMAVEC,nCIVEC,IREFSM_CASVB)

  ! Note that CI_VEC is used as a scratch array!
  use GLBBAS, only: INT1, INT1O, SCR => CI_VEC, VEC3
  use stdalloc, only: mma_allocate, mma_deallocate
  use rasscf_lucia, only: INI_H0, KVEC3_LENGTH, SIGMA_ON_DISK
  use CandS, only: ICSM, ISSM
  use lucia_data, only: NSD_PER_SYM
  use lucia_data, only: ECORE, ECORE_ORIG
  use lucia_data, only: MXNTTS
  use lucia_data, only: LUC, LUSC34
  use lucia_data, only: IREFSM

  implicit none
  integer nCIVEC, IREFSM_CASVB
  real*8 CIVEC(nCIVEC), SIGMAVEC(nCIVEC)
  integer, allocatable :: lVec(:)
  integer nSD

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
  !   call FI(INT1,ECORE_HEX,1)
  !   ECORE = ECORE+ECORE_HEX
  !end if

  ! Write CI-vector to disc

  call mma_allocate(lVec,MXNTTS,Label='lVec')
  call CPCIVC(CIVEC,nSD,LUC,MXNTTS,ISSM,1,lVec)
  call mma_deallocate(lVec)

  ! Calculate the sigma vector:

  call DIAG_MASTER()
  call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
  call MV7(SCR,SIGMAVec,LUC,LUSC34)
  call mma_deallocate(VEC3)

  ! Export lusc34 to RASSCF

  SIGMA_ON_DISK = .true.

  ! Set ICSM and ISSM (from CandS) back to IREFSM

  ICSM = IREFSM
  ISSM = IREFSM

end subroutine SIGMA_MASTER_CVB

subroutine cpcivc(CIVec,nCIVEC,ifile,mxrec,isym,iway,lrec)
  ! Copies the CI-vector between Molcas Rasscf and Lucia enviroment
  ! IWAY = 1: from Molcas to Lucia (from core to disk unit ifile).
  ! IWAY = 2: from Lucia to Molcas (from disk unit ifile to core).

  use lucia_data, only: IDISK
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  implicit none
  integer nCIVEC, ifile, mxrec, isym, iway
  integer lrec(mxrec)
  real*8 CIVec(nCIVec)
  integer nRec
# ifdef _DEBUGPRINT_
  integer iOff, iRec
# endif

  ! ==================
  ! Find nrec and lrec
  ! ==================

  call blkfo_min(isym,nrec,lrec)
  IDISK(IFILE) = 0

  ! =======================
  ! Write CI-vector to disc
  ! =======================

  if (iway == 1) then
#   ifdef _DEBUGPRINT_
    ioff = 1
    write(u6,*) 'CI-vector put to disk:'
    do IREC=1,NREC
      if (LREC(IREC) >= 0) then
        call wrtmat(CIVec(ioff),1,lrec(irec),1,lrec(irec))
        ioff = ioff+lrec(irec)
      end if
    end do
#   endif
    call todscn(CIVec,nrec,lrec,-1,ifile)
    call itods([-1],1,-1,ifile)

  else
    ! ========================
    ! Read CI-vector from disc
    ! ========================
    call frmdscn(CIVec,nrec,-1,ifile)
  end if

end subroutine cpcivc

subroutine cpsivc(ifile,mxrec,vec,lrec)

  ! Copies the Sigma-vector between Molcas Rasscf and Lucia enviroment

  use lucia_data, only: IDISK
  use general_data, only: STSYM

  implicit real*8(a-h,o-z)
  integer ifile, mxrec
  integer lrec(mxrec)
  real*8 vec(mxrec)

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
