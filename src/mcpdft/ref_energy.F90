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

subroutine ref_energy(mcscf_energy,nroots)
  use definitions,only:iwp,wp,u6
  use constants,only:zero
  use stdalloc,only:mma_allocate,mma_deallocate
  use mcpdft_input,only:mcpdft_options
  use mspdft,only:heff
  implicit none

  integer(kind=iwp),intent(in) :: nroots
  real(kind=wp),dimension(nroots),intent(out) :: mcscf_energy

#include "rasdim.fh"
#include "general.fh"

  integer(kind=iwp),dimension(15) :: iadr19
  integer(kind=iwp) :: root
  integer(kind=iwp) :: jobiph_lu,iad19,disk,nmaybe,i,it
  real(kind=wp) :: e,aemax

  real(kind=wp),allocatable :: elist(:,:)

  iadr19(:) = 0
  iad19 = 0

  if(mcpdft_options%mspdft) then
    do root = 1,nroots
      mcscf_energy(root) = heff(root,root)
    enddo
  else
    if(mcpdft_options%is_hdf5_wfn) then
      write(u6,*) 'cannot load energy from hdf5 file...'
      call abend()
    else
      call DaName(jobiph_lu,mcpdft_options%wfn_file)
      call iDaFile(mcpdft_options%wfn_file,2,iadr19,15,iad19)
      disk = iadr19(6)

      call mma_allocate(elist,mxroot,mxiter,label='EList')
      call DDaFile(jobiph_lu,2,elist,mxroot*mxiter,disk)

      nmaybe = 0
      do it = 1,mxiter
        aemax = zero
        do i = 1,mxroot
          e = elist(i,it)
          aemax = max(aemax,abs(e))
        enddo
        if(abs(aemax) <= 1.0D-12) then
          goto 11
        endif
        nmaybe = it
      enddo
11    continue

      do root = 1,nroots
        mcscf_energy(root) = elist(root,nmaybe)
      enddo

      call mma_deallocate(elist)
      call DaClos(jobiph_lu)
    endif
  endif

endsubroutine
