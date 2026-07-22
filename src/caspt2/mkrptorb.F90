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
! Copyright (C) 2006, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKRPTORB(FIFA,NFIFA,TORB,NTORB,CMO,NCMO)
! Transform to orbitals that diagonalize the diagonal blocks of FIFA.
! Affected data sets are CMO, C EPS, EPSI, EPSA, and EPSE. Also, the
! CI arrays are transformed on file LUCIEX.
! Note: FIFA is unchanged  and is not valid for the new orbitals.
! It will be recomputed later.
! The transformation matrices are returned in TORB.

use Index_Functions, only: nTri_Elem
use caspt2_qmc_interface, only: DoFCIQMC, NonDiagonal
use caspt2_global, only: IDCIEX, IDTCEX, LUCIEX
use general_data, only: STSym
use caspt2_module, only: EPS, EPSA, EPSE, EPSI, iSCF, nBas, nConf, nDel, nFro, nIsh, nOMx, nOrb, nRas1, nRas2, nRas3, nSsh, &
                         nState, nSym
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
use caspt2_module, only: DoCumulant
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_DMRG_)
use caspt2_module, only: jState, nAshT
use Constants, only: Zero
#endif
#if defined(_DMRG_)
use, intrinsic :: iso_c_binding, only: c_bool, c_int
use qcmaquis_interface, only: qcmaquis_interface_rotate_rdms
use caspt2_module, only: DMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NFIFA, NTORB, NCMO
real(kind=wp), intent(in) :: FIFA(NFIFA)
real(kind=wp), intent(out) :: TORB(NTORB)
real(kind=wp), intent(inout) :: CMO(NCMO)
integer(kind=iwp) :: I, ICMOEND, ICMOSTA, IDR, IDW, IEPS, IEPSA, IEPSE, IEPSI, II, IOEND, IOSTA, IST, ISYM, ITOEND, ITOSTA, NB, &
                     NCMOSCT, NFES, NFOCK, NO, NSCT
real(kind=wp), allocatable :: CI(:), CMO2(:), FOCK(:)
#if defined(_ENABLE_BLOCK_DMRG_) || defined(_DMRG_)
integer(kind=iwp) :: NXMAT
real(kind=wp), allocatable :: XMAT(:)
#endif

! Allocate space for temporary square Fock matrix in each symmetry:
! NBMX=Max number of basis functions in any symmetry, in common in caspt2_module
NFOCK = NOMX**2
call mma_allocate(FOCK,NFOCK,LABEL='FOCK')
! Allocate space for new CMO coefficients:
call mma_allocate(CMO2,NCMO,LABEL='CMO2')

! In the loop over symmetries, NFES is the nr of Fock matrix
! elements processed in earlier symmetries.
NFES = 0
IEPS = 0
IEPSI = 0
IEPSA = 0
IEPSE = 0
! The transformation matrices for each symmetry
! will be collected into TORB and returned. This is
! necessary for later backtransformation to original
! MO basis.
! ITOSTA,ITOEND: Section of TORB for each subspace.
ITOEND = 0
! ICMOSTA,ICMOEND: Section of CMO for each subspace.
ICMOEND = 0
do ISYM=1,NSYM
  NO = NORB(ISYM)
  NB = NBAS(ISYM)
  ! Put Fock matrix in square format in FOCK
  if (NO > 0) call SQUARE(FIFA(NFES+1),FOCK,NO,1,NO)
  ! the zero-ing out of the off-diagonal blocks happens here
  ! diafck is not structured like rasscf/fckpt2
  ! Number of orbitals processed so far in this symmetry:
  IOEND = 0
  ! Frozen orbitals: Just copy frozen CMO coefficients.
  NSCT = NFRO(ISYM)
  if (NSCT > 0) then
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    CMO2(ICMOSTA:ICMOEND) = CMO(ICMOSTA:ICMOEND)
  end if
  ! Inactive block: Section length NSCT=NISH(ISYM)
  NSCT = NISH(ISYM)
  if (NSCT > 0) then
    IOSTA = IOEND+1
    IOEND = IOEND+NSCT
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    ITOSTA = ITOEND+1
    ITOEND = ITOEND+NSCT**2
    call DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),NB,CMO(ICMOSTA),CMO2(ICMOSTA))
    do I=1,NSCT
      II = IOSTA-1+I
      IEPS = IEPS+1
      EPS(IEPS) = FOCK(II+NO*(II-1))
      IEPSI = IEPSI+1
      EPSI(IEPSI) = EPS(IEPS)
    end do
  end if
  ! RAS1 block: Section length NSCT=NRAS1(ISYM)
  NSCT = NRAS1(ISYM)
  if (NSCT > 0) then
    IOSTA = IOEND+1
    IOEND = IOEND+NSCT
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    ITOSTA = ITOEND+1
    ITOEND = ITOEND+NSCT**2
    call DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),NB,CMO(ICMOSTA),CMO2(ICMOSTA))
    do I=1,NSCT
      II = IOSTA-1+I
      IEPS = IEPS+1
      EPS(IEPS) = FOCK(II+NO*(II-1))
      IEPSA = IEPSA+1
      EPSA(IEPSA) = EPS(IEPS)
    end do
  end if
  ! RAS2 block: Section length NSCT=NRAS2(ISYM)
  NSCT = NRAS2(ISYM)
  if (NSCT > 0) then
    IOSTA = IOEND+1
    IOEND = IOEND+NSCT
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    ITOSTA = ITOEND+1
    ITOEND = ITOEND+NSCT**2
    call DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),NB,CMO(ICMOSTA),CMO2(ICMOSTA))
    do I=1,NSCT
      II = IOSTA-1+I
      IEPS = IEPS+1
      EPS(IEPS) = FOCK(II+NO*(II-1))
      IEPSA = IEPSA+1
      EPSA(IEPSA) = EPS(IEPS)
    end do
  end if
  ! RAS3 block: Section length NSCT=NRAS3(ISYM)
  NSCT = NRAS3(ISYM)
  if (NSCT > 0) then
    IOSTA = IOEND+1
    IOEND = IOEND+NSCT
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    ITOSTA = ITOEND+1
    ITOEND = ITOEND+NSCT**2
    call DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),NB,CMO(ICMOSTA),CMO2(ICMOSTA))
    do I=1,NSCT
      II = IOSTA-1+I
      IEPS = IEPS+1
      EPS(IEPS) = FOCK(II+NO*(II-1))
      IEPSA = IEPSA+1
      EPSA(IEPSA) = EPS(IEPS)
    end do
  end if
  ! Secondary (virtual) block: Section length NSCT=NSSH(ISYM)
  NSCT = NSSH(ISYM)
  if (NSCT > 0) then
    IOSTA = IOEND+1
    IOEND = IOEND+NSCT
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    ITOSTA = ITOEND+1
    ITOEND = ITOEND+NSCT**2
    call DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),NB,CMO(ICMOSTA),CMO2(ICMOSTA))
    do I=1,NSCT
      II = IOSTA-1+I
      IEPS = IEPS+1
      EPS(IEPS) = FOCK(II+NO*(II-1))
      IEPSE = IEPSE+1
      EPSE(IEPSE) = EPS(IEPS)
    end do
  end if
  ! Deleted orbitals: Just copy deleted CMO coefficients.
  NSCT = NDEL(ISYM)
  if (NSCT > 0) then
    NCMOSCT = NBAS(ISYM)*NSCT
    ICMOSTA = ICMOEND+1
    ICMOEND = ICMOEND+NCMOSCT
    CMO2(ICMOSTA:ICMOEND) = CMO(ICMOSTA:ICMOEND)
  end if

  NFES = NFES+nTri_Elem(NO)

end do

! Actually, we will not use the old CMO array any more, so just
! overwrite it with the new ones and get rid of the allocated array.
CMO(:) = CMO2(:)
call mma_deallocate(CMO2)

! We will not use the Fock matrix either. It was just used temporarily
! for each turn of the symmetry loop. Skip it.
call mma_deallocate(FOCK)

! Finally, loop again over symmetries, transforming the CI:

if (ISCF == 0) then
# ifdef _DMRG_
  if (DMRG) then
    NXMAT = NASHT**2
    call mma_allocate(XMAT,NXMAT,LABEL='XMAT')
    XMAT(:) = Zero
    call MKXMAT(TORB,XMAT)

    call qcmaquis_interface_rotate_rdms(int(JSTATE-1,kind=c_int),int(JSTATE-1,kind=c_int),0_c_int,XMAT,.false._c_bool)
    do I=1,NSTATE
      if (JSTATE /= I) then
        write(u6,*) 'QCMaquis> Rotating tRDMs',JSTATE-1,I-1
        call qcmaquis_interface_rotate_rdms(int(JSTATE-1,kind=c_int),int(I-1,kind=c_int),0_c_int,XMAT,.false._c_bool)
      end if
    end do
    call mma_deallocate(XMAT)
  end if
# endif

  if (DoFCIQMC) then
    if (NonDiagonal) then
      write(u6,*) 'Transforming CASPT2 intermediates to pseudo-canonical orbitals.'
    else
      write(u6,*) 'FCIQMC-CASPT2 assumes pseudo-canonical orbitals.'
    end if
  else
#   if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
    if (.not. DoCumulant) then
#   endif
      call mma_allocate(CI,NCONF,Label='CI')
      do IST=1,NSTATE
        IDR = IDCIEX(IST)
        call DDAFILE(LUCIEX,2,CI,NCONF,IDR)

        call mkTraCI(nTORB,TORB,STSYM,nConf,CI)

        IDW = IDTCEX(IST)
        call DDAFILE(LUCIEX,1,CI,NCONF,IDW)
      end do
      call mma_deallocate(CI)
#   ifdef _ENABLE_BLOCK_DMRG_
    else
      ! Transforming 2,3-RDMs from Block DMRG (1-RDM is computed from 2-RDM)
      ! NN.14 : For the time, Block's dump files of RDMs are directly loaded,
      !         but those should be stored in JobIph file eventually.
      NXMAT = NASHT**2
      ! Workspace for transformation matrix
      call mma_allocate(XMAT,NXMAT,LABEL='XMAT')
      XMAT(:) = Zero
      call MKXMAT(TORB,XMAT)

      call block_tran2pdm(NASHT,XMAT,JSTATE,JSTATE)
      call block_tran3pdm(NASHT,XMAT,JSTATE,JSTATE)

      call mma_deallocate(XMAT)
    end if
#   elif _ENABLE_CHEMPS2_DMRG_
  else
    write(u6,*) 'CHEMPS2> MKRPTORB assumes PSEUDOCANONICAL orbitals!'
  end if
# endif
end if
end if

end subroutine MKRPTORB
