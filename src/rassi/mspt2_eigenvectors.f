************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2017, Stefan Knecht                                    *
************************************************************************
      module mspt2_eigenvectors

        implicit none

        type mspt2evc
          real*8, allocatable :: pc(:,:)
          real*8, allocatable :: sc(:,:)
        end Type

        type(mspt2evc), public, allocatable :: Heff_evc(:)

        save

        contains

        subroutine init_mspt2_eigenvectors(ijob,nstates,tag)
          integer, intent(in) :: ijob, nstates, tag
          if(tag == 0)then
            allocate(Heff_evc(ijob))
          else if(tag == 1)then
            allocate(Heff_evc(ijob)%pc(nstates,nstates))
            Heff_evc(ijob)%pc = 0
          else if(tag == 2)then
            allocate(Heff_evc(ijob)%sc(nstates,nstates))
            Heff_evc(ijob)%sc = 0
          else
            write(6,*) 'unknown tag in init_mspt2_eigenvectors'
            call Abend()
          end if
        end subroutine init_mspt2_eigenvectors

        subroutine deinit_mspt2_eigenvectors()
          integer             :: i
          do i = 1, size(Heff_evc)
            if(allocated(Heff_evc(i)%pc)) deallocate(Heff_evc(i)%pc)
            if(allocated(Heff_evc(i)%sc)) deallocate(Heff_evc(i)%sc)
          end do
          deallocate(Heff_evc)
        end subroutine deinit_mspt2_eigenvectors

        subroutine prpdata_mspt2_eigenvectors(
     &                                         rtdm,
     &                                         stdm,
     &                                         wetdm,
     &                                         prop,
     &                                         nprop,
     &                                         nstate,
     &                                         istate,
     &                                         jstate,
     &                                         ntdmzz,
     &                                         addr,
     &                                         iempty,
     &                                         lu,
     &                                         put_so_data,
     &                                         put_h5_data
     &                                       )

        integer, intent(in)    :: nprop
        integer, intent(in)    :: nstate
        integer, intent(in)    :: istate
        integer, intent(in)    :: jstate
        integer, intent(in)    :: ntdmzz
        integer, intent(in)    :: addr
        integer, intent(in)    :: iempty
        integer, intent(in)    :: lu
        logical, intent(in)    :: put_so_data
        logical, intent(in)    :: put_h5_data
        real*8,  intent(in)    :: rtdm(ntdmzz)
        real*8,  intent(in)    :: stdm(ntdmzz)
        real*8,  intent(in)    :: wetdm(ntdmzz)
        real*8,  intent(inout) :: prop(nstate,nstate,nprop)
        integer iOpt, iaddr, iGo
#include "rassiwfn.fh"


          !> calculate property matrix elements
          call proper(prop,istate,jstate,rtdm,wetdm)

          !> put data to file
          if(put_so_data)then
            iOpt=1
            iGo=7
            iaddr=addr
            call dens2file(
     &                     rtdm,
     &                     stdm,
     &                     wetdm,
     &                     ntdmzz,
     &                     lu,
     &                     iaddr,
     &                     iempty,
     &                     iOpt,
     &                     iGo,
     &                     iState,
     &                     jState
     &                     )
          end if

          if(put_h5_data.or.put_so_data)then
#ifdef _HDF5_
            call mh5_put_dset_array_real(wfn_sfs_tdm,
     &                                   rtdm, [NTDMZZ,1,1],
     &                                   [0,ISTATE-1,JSTATE-1])
            call mh5_put_dset_array_real(wfn_sfs_tsdm,
     &                                   stdm, [NTDMZZ,1,1],
     &                                   [0,ISTATE-1,JSTATE-1])
            if(put_so_data)then
              call mh5_put_dset_array_real(wfn_sfs_wetdm,
     &                                     wetdm, [NTDMZZ,1,1],
     &                                     [0,ISTATE-1,JSTATE-1])
            end if
#endif

          end if
        end subroutine prpdata_mspt2_eigenvectors

      end module mspt2_eigenvectors
