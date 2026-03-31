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
! Copyright (C) 2017, Stefan Knecht                                    *
!***********************************************************************

module mspt2_eigenvectors

#ifdef _HDF5_
use mh5, only: mh5_put_dset
#endif
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

real(kind=wp), allocatable :: Heff_evc_pc(:,:,:), Heff_evc_sc(:,:,:)

public :: deinit_mspt2_eigenvectors, Heff_evc_pc, Heff_evc_sc, init_mspt2_eigenvectors, prpdata_mspt2_eigenvectors

contains

subroutine init_mspt2_eigenvectors(ijob,nstates,tag)

  use stdalloc, only: mma_allocate
  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: ijob, nstates, tag

  select case (tag)
    case (0)
    case (1)
      call mma_allocate(Heff_evc_pc,nstates,nstates,ijob,Label='Heff_evc_pc')
      HEff_evc_pc(:,:,:) = Zero
    case (2)
      call mma_allocate(Heff_evc_sc,nstates,nstates,ijob,Label='Heff_evc_sc')
      HEff_evc_sc(:,:,:) = Zero
    case default
      write(u6,*) 'unknown tag in init_mspt2_eigenvectors'
      call Abend()
  end select

end subroutine init_mspt2_eigenvectors

subroutine deinit_mspt2_eigenvectors()

  use stdalloc, only: mma_deallocate

  if (allocated(Heff_evc_pc)) call mma_deallocate(Heff_evc_pc)
  if (allocated(Heff_evc_sc)) call mma_deallocate(Heff_evc_sc)

end subroutine deinit_mspt2_eigenvectors

subroutine prpdata_mspt2_eigenvectors(rtdm,stdm,wetdm,prop,nprop,nstate,istate,jstate,ntdmzz,addr,iempty,lu,put_so_data,put_h5_data)

# ifdef _HDF5_
  use RASSIWfn, only: wfn_SFS_TDM, wfn_SFS_TSDM, wfn_SFS_WETDM
# endif

  integer(kind=iwp), intent(in) :: nprop, nstate, istate, jstate, ntdmzz, addr, iempty, lu
  real(kind=wp), intent(inout) :: rtdm(ntdmzz), stdm(ntdmzz), wetdm(ntdmzz), prop(nstate,nstate,nprop)
  logical(kind=iwp), intent(in) :: put_so_data, put_h5_data
  integer(kind=iwp) :: iaddr, iGo, iOpt

  !> calculate property matrix elements
  call proper(prop,istate,jstate,rtdm,wetdm)

  !> put data to file
  if (put_so_data) then
    iOpt = 1
    iGo = 7
    iaddr = addr
    call dens2file(rtdm,stdm,wetdm,ntdmzz,lu,iaddr,iempty,iOpt,iGo,iState,jState)
  end if

# ifdef _HDF5_
  if (put_h5_data .or. put_so_data) then
    call mh5_put_dset(wfn_sfs_tdm,rtdm,[NTDMZZ,1,1],[0,ISTATE-1,JSTATE-1])
    call mh5_put_dset(wfn_sfs_tsdm,stdm,[NTDMZZ,1,1],[0,ISTATE-1,JSTATE-1])
    if (put_so_data) call mh5_put_dset(wfn_sfs_wetdm,wetdm,[NTDMZZ,1,1],[0,ISTATE-1,JSTATE-1])
  end if
# else
# include "macros.fh"
  unused_var(put_h5_data)
# endif

end subroutine prpdata_mspt2_eigenvectors

end module mspt2_eigenvectors
