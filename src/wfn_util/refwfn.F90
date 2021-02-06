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

module refwfn

use Definitions, only: wp, iwp, u6

implicit none
private

logical(kind=iwp) :: refwfn_active = .false., refwfn_is_h5
character(len=128) :: refwfn_filename
integer(kind=iwp) :: refwfn_id, IADR15(30)
character(len=100) :: ProgName
character(len=100), external :: Get_ProgName

public :: refwfn_active, refwfn_is_h5, refwfn_filename, refwfn_id, IADR15
public :: refwfn_init, refwfn_close, refwfn_info, refwfn_data

contains

!***********************************************************************
subroutine refwfn_init(Filename)
!***********************************************************************

# ifdef _HDF5_
  use mh5, only: mh5_is_hdf5, mh5_open_file_r
# endif

  character(len=*), intent(in) :: Filename
  integer(kind=iwp) :: I, IAD15

  refwfn_is_h5 = .false.

  ProgName = Get_ProgName()
  if (refwfn_active) then
    write(u6,*) ' trying to activate refwfn twice, aborting!'
    call abend()
  else
    refwfn_active = .true.
  end if

  ! if not a standard filename, call fileorb??
  if (FileName /= 'JOBIPH') then
    call fileorb(Filename,refwfn_filename)
  else
    refwfn_filename = 'JOBIPH'
  end if

# ifdef _HDF5_
  if (mh5_is_hdf5(refwfn_filename)) then
    refwfn_is_h5 = .true.
    write(u6,'(1X,A)') 'wavefunction data from HDF5 file:'
    write(u6,'(3X,A)') trim(refwfn_filename)
    refwfn_id = mh5_open_file_r(refwfn_filename)
  else
# endif
    refwfn_is_h5 = .false.
    ! Assume reference wavefunction is stored as JobIph format
    refwfn_id = 15
    call DANAME(refwfn_id,refwfn_filename)
    ! Read table of contents into IADR15() array.
    ! There are two possible different layouts, 15 or 30 integers:
    IAD15 = 0
    call IDAFILE(refwfn_id,2,IADR15,15,IAD15)
    if (IADR15(15) == -1) then
      IAD15 = 0
      call IDAFILE(refwfn_id,2,IADR15,30,IAD15)
    else
      do I=16,30
        IADR15(I) = 0
      end do
      call WarningMessage(1,'Old JOBIPH file layout.')
    end if
# ifdef _HDF5_
  end if
# endif

end subroutine refwfn_init

!***********************************************************************
subroutine refwfn_close
!***********************************************************************

# ifdef _HDF5_
  use mh5, only: mh5_close_file

  if (refwfn_is_h5) then
    call mh5_close_file(refwfn_id)
  else
# endif
    call DaClos(refwfn_id)
# ifdef _HDF5_
  end if
# endif
  refwfn_active = .false.

end subroutine refwfn_close

!***********************************************************************
subroutine refwfn_info
!***********************************************************************
!SVC: initialize the reference wavefunction info

  use stdalloc, only: mma_allocate, mma_deallocate
# ifdef _DMRG_
  use qcmaquis_info, only: qcmaquis_info_init, qcm_group_names
# endif
# ifdef _HDF5_
  use mh5, only: mh5_fetch_attr, mh5_exists_dset, mh5_fetch_dset
# endif

# include "rasdim.fh"
# include "caspt2.fh"

# ifdef _HDF5_
  character(len=1), allocatable :: typestring(:)
# endif
  integer(kind=iwp) :: iSym, ref_nSym, ref_nBas(mxSym), IAD15
  real(kind=wp) :: Weight(mxRoot)

  if (.not. refwfn_active) then
    write(u6,*) ' refwfn not yet activated, aborting!'
    call abend()
  end if

# ifdef _HDF5_
  if (refwfn_is_h5) then
    ! general wavefunction attributes
    !call mh5_fetch_attr(refwfn_id, 'TITLE', Title)
    call mh5_fetch_attr(refwfn_id,'SPINMULT',iSpin)
    call mh5_fetch_attr(refwfn_id,'NSYM',ref_nSym)
    call mh5_fetch_attr(refwfn_id,'LSYM',stSym)
    call mh5_fetch_attr(refwfn_id,'NBAS',ref_nBas)

    call mh5_fetch_attr(refwfn_id,'NACTEL',nActEl)
    call mh5_fetch_attr(refwfn_id,'NHOLE1',nHole1)
    call mh5_fetch_attr(refwfn_id,'NELEC3',nEle3)
    call mh5_fetch_attr(refwfn_id,'NCONF',nConf)
    call mh5_fetch_attr(refwfn_id,'NSTATES',nRoots)
    call mh5_fetch_attr(refwfn_id,'NROOTS',lRoots)
    call mh5_fetch_attr(refwfn_id,'STATE_ROOTID',iRoot)
    call mh5_fetch_attr(refwfn_id,'STATE_WEIGHT',Weight)

    call mma_allocate(typestring,sum(ref_nbas(1:nsym)))
    call mh5_fetch_dset(refwfn_id,'MO_TYPEINDICES',typestring)
    call tpstr2orb(ref_nsym,ref_nbas,typestring,nfro,nish,nras1,nras2,nras3,nssh,ndel)
    nash = nras1+nras2+nras3
    call mma_deallocate(typestring)
    ! Leon 14/6/2017 -- do not read CI vectors if NEVPT2 is attempted
    ! because for now we only support DMRG-NEVPT2
    if (ProgName(1:6) == 'caspt2') then
      if (.not. mh5_exists_dset(refwfn_id,'CI_VECTORS')) then
        write(u6,'(1X,A)') 'The HDF5 file does not contain CI vectors,'
        write(u6,'(1X,A)') 'make sure it was created by rasscf/caspt2.'
        call AbEnd()
      end if
    end if
    if (.not. mh5_exists_dset(refwfn_id,'MO_VECTORS')) then
      write(u6,'(1X,A)') 'The HDF5 file does not contain MO vectors,'
      write(u6,'(1X,A)') 'make sure it was created by rasscf/caspt2/nevpt2.'
      call AbEnd()
    end if
    IFQCAN = 0
# ifdef _DMRG_
    if (mh5_exists_dset(refwfn_id,'QCMAQUIS_CHECKPOINT')) then
      call qcmaquis_info_init(1,nroots,-1)
      call mh5_fetch_dset(refwfn_id,'QCMAQUIS_CHECKPOINT',qcm_group_names(1)%states)
    end if
# endif
  else
# endif
    ! Sizes in the GSLIST is counted in INTEGERS.
    ! Note that the title field in the JOBIPH file is not used for anything
    ! in this program, it is just a dummy read.
    ! Another title field is read from input a little later, it is called
    ! TITLE2. That one is printed out in PRINP_CASPT2.
    IAD15 = IADR15(1)
    call WR_RASSCF_Info(refwfn_id,2,iAd15,NACTEL,ISPIN,REF_NSYM,STSYM,NFRO,NISH,NASH,NDEL,REF_NBAS,8,NAME,LENIN8*MXORB,NCONF, &
                        HEADER,144,TITLE,4*18*mxTit,POTNUC,LROOTS,NROOTS,IROOT,MXROOT,NRAS1,NRAS2,NRAS3,NHOLE1,NELE3,IFQCAN,Weight)
    nssh = ref_nbas-nfro-nish-nash-ndel
# ifdef _HDF5_
  end if
# endif
  if (nSym /= ref_nSym) then
    write(u6,*) ' Number of irreps of the reference wavefunction'
    write(u6,*) ' does not match the data on the RunFile, abort!'
    call AbEnd()
  else
    do iSym=1,nSym
      if (nBas(iSym) /= ref_nBas(iSym)) then
        write(u6,*) ' Number of basis functions of the reference'
        write(u6,*) ' wavefunction does not match the data on the'
        write(u6,*) ' RunFile, abort!'
        call AbEnd()
      end if
    end do
  end if

end subroutine refwfn_info

!***********************************************************************
subroutine refwfn_data
!***********************************************************************
!SVC: initialize the reference wavefunction data

  use stdalloc, only: mma_allocate, mma_deallocate
# ifdef _HDF5_
  use mh5, only: mh5_fetch_attr, mh5_fetch_dset, mh5_fetch_dset_array_real
# endif

# include "rasdim.fh"
# include "caspt2.fh"
# include "pt2_guga.fh"

  integer(kind=iwp) :: I, IAD15, II, IDISK, ID, IAD, NEJOB, IT, NMAYBE, ISNUM
  real(kind=wp) :: Root_Energies(mxRoot), AEMAX, E
  real(kind=wp), allocatable :: tmp(:), ejob(:,:)

  if (.not. refwfn_active) then
    write(u6,*) ' refwfn not yet activated, aborting!'
    call abend()
  end if

  !---  Read the MO coefficients from HDF5/JOBIPH and store on LUONEM
  NCMO = NBSQT
  call mma_allocate(tmp,NCMO,label='LCMORAS')
# ifdef _HDF5_
  if (refwfn_is_h5) then
    call mh5_fetch_dset_array_real(refwfn_id,'MO_VECTORS',tmp)
  else
# endif
    IAD15 = IADR15(9)
    if (IFQCAN == 0) IAD15 = IADR15(2)
    call DDAFILE(refwfn_id,2,tmp,NCMO,IAD15)
# ifdef _HDF5_
  end if
# endif
  IEOF1M = 0
  IDISK = IEOF1M
  IAD1M(1) = IDISK
  call DDAFILE(LUONEM,1,tmp,NCMO,IDISK)
  call mma_deallocate(tmp)
  IEOF1M = IDISK

  ! IDCIEX: Present EOF on LUCIEX.
  ID = IDCIEX
  ! Skip when using cumulant reconstruction of (3-,) 4-RDM
  !     Leon 14/6/2017 -- do not read CI vectors if NEVPT2 is attempted
  !     because for now we only support DMRG-NEVPT2
  if (ProgName(1:6) == 'caspt2') then
    if ((.not. DoCumulant) .and. (ISCF == 0)) then
      call mma_allocate(tmp,NCONF,label='LCI')
      do I=1,NSTATE
        ISNUM = MSTATE(I)
#       ifdef _HDF5_
        if (refwfn_is_h5) then
          !---  Read the CI coefficients from the HDF5 file
          call mh5_fetch_dset_array_real(refwfn_id,'CI_VECTORS',tmp,[nconf,1],[0,ISNUM-1])
        else
#       endif
          !---  Read the CI coefficients from the JOBIPH file
          IDISK = IADR15(4)
          do II=1,ISNUM-1
            call DDAFILE(refwfn_id,0,tmp,NCONF,IDISK)
          end do
          call DDAFILE(refwfn_id,2,tmp,NCONF,IDISK)
#       ifdef _HDF5_
        end if
#       endif
        ! Copy selected vectors to LUCI:
        call DDAFILE(LUCIEX,1,tmp,NCONF,ID)
      end do
      ! Disk address = present EOF on LUCIEX.
      ! IDTCEX = Disk address to transformed CI.
      if (ORBIN == 'TRANSFOR') then
        IDTCEX = ID
        ! Dummy writes:
        do II=1,NSTATE
          call DDAFILE(LUCIEX,0,tmp,NCONF,ID)
        end do
      else
        IDTCEX = IDCIEX
      end if
      call mma_deallocate(tmp)
    else
      ! If this is Closed-shell or Hi-spin SCF case
      ! Just in case...
      if (.not. DoCumulant .and. (NSTATE /= 1 .or. NCONF /= 1)) then
        write(u6,*) ' readin_caspt2: A Closed-shell or Hi-spin SCF'
        write(u6,*) ' but nr of states is: NSTATE=',NSTATE
        write(u6,*) ' and nr of CSFs is    NCONF= ',NCONF
        write(u6,*) ' Program error?? Must stop.'
        call ABEND()
      end if
      ! This should be solved elsewhere in the code...just for the now,
      ! make a write of a CI vector to LUCIEX, so other routines do not get
      ! their knickers into a twist:
      NCONF = 1
      call mma_allocate(tmp,NCONF,label='LCI')
      tmp(1) = 1.0_wp
      call DDAFILE(LUCIEX,1,tmp,NCONF,ID)
      call mma_deallocate(tmp)
    end if
  end if
  ! Now, the selected original CASCI expansions are on LUCIEX
  ! beginning from disk address 0.

  !SVC: read the L2ACT and LEVEL arrays
# ifdef _HDF5_
  if (refwfn_is_h5) then
    call mh5_fetch_attr(refwfn_id,'L2ACT',L2ACT)
    call mh5_fetch_attr(refwfn_id,'A2LEV',LEVEL)
  else
# endif
    IAD15 = IADR15(18)
    call IDAFILE(refwfn_id,2,L2ACT,mxAct,IAD15)
    call IDAFILE(refwfn_id,2,LEVEL,mxAct,IAD15)
# ifdef _HDF5_
  end if
# endif

# ifdef _HDF5_
  if (refwfn_is_h5) then
    call mh5_fetch_dset(refwfn_id,'ROOT_ENERGIES',ROOT_ENERGIES)
  else
# endif
    !PAM 2015: We will no longer recompute the RASSCF energies
    ! but just assume they can be obtained from the JOBIPH file.
    ! There is table with unknown length with all energies from
    ! all iterations (!) there.
    NEJOB = MXROOT*MXITER
    call mma_allocate(ejob,MXROOT,MXITER,label='EJOB')
    IAD = IADR15(6)
    call DDAFILE(refwfn_id,2,ejob,NEJOB,IAD)
    ! Note that there is no info on nr of iterations
    ! so we cannot know what energies to pick...
    ! Let us make a guess: The correct set of energy values in the
    ! table of energies/iteration is the last one with not all zeroes.
    NMAYBE = 0
    do IT=1,MXITER
      AEMAX = 0.0_wp
      do I=1,MXROOT
        E = ejob(I,IT)
        AEMAX = max(AEMAX,abs(E))
      end do
      if (abs(AEMAX) < 1.0e-12_wp) exit
      NMAYBE = IT
    end do
    if (NMAYBE == 0) then
      write(u6,*) ' PT2INI tried to read energies from the'
      write(u6,*) ' JOBIPH file, but could not find any.'
      call ABEND()
    end if
    ! And then put the energies into the Hamiltonian matrix,
    ! unless already filled in by the EFFE keyword
    do I=1,mxRoot
      Root_Energies(I) = ejob(I,NMAYBE)
    end do
    ! No more use for the array EJOB.
    call mma_deallocate(ejob)
# ifdef _HDF5_
  end if
# endif
  ! Leon: MSTATE(I) is initialised only in caspt2
  ! if it's initialised somewhere else, feel free to change the line below
  if (ProgName(1:6) == 'caspt2') then
    do I=1,NSTATE
      REFENE(I) = ROOT_ENERGIES(MSTATE(I))
    end do
  else
    ! since it isn't initialised, we just assume that we need the states in the
    ! same order as they're written on JobIph
    ! Apparently, nstates is also initialised in caspt2 only, but nroots is
    ! available from here, so we use nroots
    ! > stknecht: initialize nstate here nevertheless - can be used in nevpt2 initialization as well
    nstate = nroots
    refene(1:nroots) = ROOT_ENERGIES(1:nroots)
  end if

end subroutine refwfn_data

end module refwfn
