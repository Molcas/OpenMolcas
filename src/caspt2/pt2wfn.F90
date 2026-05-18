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
! Copyright (C) 2016, Steven Vancoillie                                *
!               2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

module pt2wfn

use Definitions, only: wp, iwp

implicit none
private

# ifdef _HDF5_
integer(kind=iwp) :: pt2wfn_cicoef, pt2wfn_dens, pt2wfn_energy, pt2wfn_heff, pt2wfn_id, pt2wfn_mocoef, pt2wfn_occnum, &
                     pt2wfn_orbene, pt2wfn_refene
#endif
logical(kind=iwp) :: pt2wfn_is_h5 = .false.

public :: pt2wfn_close, pt2wfn_data, pt2wfn_densstore, pt2wfn_estore, pt2wfn_init

contains

subroutine pt2wfn_init()
  ! SVC: Create a wavefunction file. If another .wfn file already
  ! exists, it will be overwritten.

  use Molcas, only: MxAct
  use refwfn, only: refwfn_active
# ifdef _HDF5_
  use refwfn, only: refwfn_is_h5
  use mh5, only: mh5_close_dset, mh5_create_dset_real, mh5_create_dset_str, mh5_create_file, mh5_init_attr, mh5_put_dset
  use sguga, only: L2ACT, LEVEL
  use caspt2_global, only: do_grad
  use caspt2_module, only: DMRG, IfMix, IfMSCOUP, IfProp, iSpin, lRoots, mState, nActEl, nBas, nBasT, nBSqT, nConf, nDel, nDet, &
                           nEle3, nFro, nHole1, nIsh, nOrb, nRas1, nRas1T, nRas2, nRas3, nRas3T, nSsh, nState, nSym, Root2State, &
                           STSym
  use stdalloc, only: mma_allocate, mma_deallocate
# endif

# ifdef _HDF5_
  integer(kind=iwp) :: dsetid, i, ndmat
  character, allocatable :: typestring(:)
# endif

  if (refwfn_active) then
    call WarningMessage(2,'Active reference wavefunction file, cannot create new PT2 wavefunction file, aborting!')
    call AbEnd()
  end if

# ifdef _HDF5_
  if (refwfn_is_h5) then
    pt2wfn_is_h5 = .true.
    ! create a new wavefunction file!
    pt2wfn_id = mh5_create_file('PT2WFN')
    ! no JOBMIX for HDF5 files
    ifmix = .false.

    ! set module type
    call mh5_init_attr(pt2wfn_id,'MOLCAS_MODULE','CASPT2')

    ! copy basic molecular information to the HDF5 file
    call run2h5_molinfo(pt2wfn_id)
    call one2h5_ovlmat(pt2wfn_id,nsym,nbas)
    call one2h5_fckint(pt2wfn_id,nsym,nbas)
    call one2h5_crtmom(pt2wfn_id,nsym,nbas)

    ! set wavefunction type
    if (NRAS1T+NRAS3T == 0) then
      call mh5_init_attr(pt2wfn_id,'CI_TYPE','CAS')
    else
      call mh5_init_attr(pt2wfn_id,'CI_TYPE','RAS')
    end if

    ! general wavefunction attributes
    call mh5_init_attr(pt2wfn_id,'SPINMULT',iSpin)
    call mh5_init_attr(pt2wfn_id,'LSYM',stSym)
    call mh5_init_attr(pt2wfn_id,'NACTEL',nActEl)
    call mh5_init_attr(pt2wfn_id,'NHOLE1',nHole1)
    call mh5_init_attr(pt2wfn_id,'NELEC3',nEle3)
    call mh5_init_attr(pt2wfn_id,'NCONF',nConf)
    call mh5_init_attr(pt2wfn_id,'NSTATES',NSTATE)
    call mh5_init_attr(pt2wfn_id,'NDET',NDET)

    call mh5_init_attr(pt2wfn_id,'L2ACT',1,[mxAct],L2ACT)
    call mh5_init_attr(pt2wfn_id,'A2LEV',1,[mxAct],LEVEL)

    ! molecular orbital type index
    call mma_allocate(typestring,nbast)
    call orb2tpstr(NSYM,NBAS,NFRO,NISH,NRAS1,NRAS2,NRAS3,NSSH,NDEL,typestring)
    dsetid = mh5_create_dset_str(pt2wfn_id,'MO_TYPEINDICES',1,[NBAST],1)
    call mh5_init_attr(dsetid,'DESCRIPTION', &
                       'Type index of the molecular orbitals arranged as blocks of size [NBAS(i)], i=1,#irreps')
    call mh5_put_dset(dsetid,typestring)
    call mma_deallocate(typestring)
    call mh5_close_dset(dsetid)

    ! roots
    call mh5_init_attr(pt2wfn_id,'STATE_ROOTID',1,[NSTATE],MSTATE)
    call mh5_init_attr(pt2wfn_id,'ROOT2STATE',1,[LROOTS],ROOT2STATE)

    ! reference energy (for each CI root)
    pt2wfn_refene = mh5_create_dset_real(pt2wfn_id,'STATE_REFWF_ENERGIES',1,[NSTATE])
    call mh5_init_attr(pt2wfn_refene,'DESCRIPTION','Reference energy for each state, arranged as array of [NSTATES]')

    ! PT2 energy (for each CI root)
    pt2wfn_energy = mh5_create_dset_real(pt2wfn_id,'STATE_PT2_ENERGIES',1,[NSTATE])
    call mh5_init_attr(pt2wfn_energy,'DESCRIPTION','PT2 energy for each state, arranged as array of [NSTATES]')

    ! molecular orbital coefficients
    pt2wfn_mocoef = mh5_create_dset_real(pt2wfn_id,'MO_VECTORS',1,[NBSQT])
    call mh5_init_attr(pt2wfn_mocoef,'DESCRIPTION', &
                       'Coefficients of the average orbitals, arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

    ! molecular orbital occupation numbers
    ! (most probably empty, but left for compatibility)
    pt2wfn_occnum = mh5_create_dset_real(pt2wfn_id,'MO_OCCUPATIONS',1,[NBAST])
    call mh5_init_attr(pt2wfn_occnum,'DESCRIPTION', &
                       'Occupation numbers of the average orbitals arranged as blocks of size [NBAS(i)], i=1,#irreps')

    ! molecular orbital energies
    ! (most probably empty, but left for compatibility)
    pt2wfn_orbene = mh5_create_dset_real(pt2wfn_id,'MO_ENERGIES',1,[NBAST])
    call mh5_init_attr(pt2wfn_orbene,'DESCRIPTION', &
                       'Orbital energies of the average orbitals arranged as blocks of size [NBAS(i)], i=1,#irreps')

    ! CI data for each root
    if (.not. DMRG) then
      pt2wfn_cicoef = mh5_create_dset_real(pt2wfn_id,'CI_VECTORS',2,[nConf,NSTATE])
      call mh5_init_attr(pt2wfn_cicoef,'DESCRIPTION', &
                         'Coefficients of configuration state functions in Split-GUGA ordering for each STATE, arranged as '// &
                         'matrix of size [NCONF,NSTATES]')
    end if

    ! effective Hamiltonian coefficients
    if (IFMSCOUP) then
      pt2wfn_heff = mh5_create_dset_real(pt2wfn_id,'H_EFF',2,[NSTATE,NSTATE])
      call mh5_init_attr(pt2wfn_heff,'DESCRIPTION', &
                         'Effective (X)MS-CASPT2 Hamiltonian, arranged as matrix of size [NSTATES,NSTATES]')
    end if

    ! density matrices
    if (IFPROP .or. do_grad) then
      ndmat = 0
      do i=1,NSYM
        ndmat = ndmat+(NORB(i)**2+NORB(i))/2
      end do

      pt2wfn_dens = mh5_create_dset_real(pt2wfn_id,'DENSITY_MATRIX',2,[ndmat,NSTATE])
      call mh5_init_attr(pt2wfn_dens,'DESCRIPTION', &
                         '1-body density matrix, arranged as blocks of size NDMAT=sum([NORB(i)*(NORB(i)+1)/2], i=1,#irreps), '// &
                         'where NORB excludes frozen and deleted orbitals, for each state: [NDMAT,NSTATES].')
    end if

  else
# endif
    pt2wfn_is_h5 = .false.
# ifdef _HDF5_
  end if
# endif

end subroutine pt2wfn_init

subroutine pt2wfn_data()

# ifdef _HDF5_
  use mh5, only: mh5_put_dset
  use caspt2_global, only: IDCIEX, LUCIEX, LUONEM, NCMO
  use caspt2_module, only: DMRG, iAd1m, nConf, nState
  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), allocatable :: BUF(:)
  integer(kind=iwp) :: IDISK, ISTATE

  if (pt2wfn_is_h5) then
    if (.not. DMRG) then
      call mma_allocate(BUF,NCONF)
      do ISTATE=1,NSTATE
        IDISK = IDCIEX(ISTATE)
        call DDAFILE(LUCIEX,2,BUF,NCONF,IDISK)
        call mh5_put_dset(pt2wfn_cicoef,BUF,[NCONF,1],[0,ISTATE-1])
      end do
      call mma_deallocate(BUF)
    end if

    call mma_allocate(BUF,NCMO)
    IDISK = IAD1M(1)
    call DDAFILE(LUONEM,2,BUF,NCMO,IDISK)
    call mh5_put_dset(pt2wfn_mocoef,BUF)
    call mma_deallocate(BUF)
  end if
# endif

end subroutine pt2wfn_data

subroutine pt2wfn_estore(Heff,nState)

# ifdef _HDF5_
  use mh5, only: mh5_put_dset
  use caspt2_module, only: Energy, IfMSCOUP, RefEne
# endif

  integer(kind=iwp), intent(in) :: nstate
  real(kind=wp), intent(in) :: Heff(nstate,nstate)

# ifdef _HDF5_
  if (pt2wfn_is_h5) then
    call mh5_put_dset(pt2wfn_energy,ENERGY)
    call mh5_put_dset(pt2wfn_refene,REFENE)
    if (IFMSCOUP) call mh5_put_dset(pt2wfn_heff,Heff)
  end if
# else
  ! Avoid unused argument warnings
  if (.false.) call Unused_real_array(Heff)
#endif

end subroutine pt2wfn_estore

subroutine pt2wfn_densstore(Dmat,nDmat)

# ifdef _HDF5_
  use mh5, only: mh5_put_dset
  use caspt2_module, only: jState
# endif

  integer(kind=iwp), intent(in) :: nDmat
  real(kind=wp), intent(in) :: Dmat(nDmat)

# ifdef _HDF5_
  if (pt2wfn_is_h5) call mh5_put_dset(pt2wfn_dens,Dmat,[nDmat,1],[0,JSTATE-1])
# else
  ! Avoid unused argument warnings
  if (.false.) call Unused_real_array(Dmat)
# endif

end subroutine pt2wfn_densstore

subroutine pt2wfn_close()

# ifdef _HDF5_
  use mh5, only: mh5_close_file

  if (pt2wfn_is_h5) call mh5_close_file(pt2wfn_id)
  pt2wfn_id = -1
#endif

end subroutine pt2wfn_close

end module pt2wfn
