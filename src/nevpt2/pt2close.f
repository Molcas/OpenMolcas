************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************

      !> Cleanup and close files
      subroutine pt2close

      use nevpt2_cfg
      use info_state_energy  ! energies
      use info_orbital_space ! orbital specifications read from JobIph
      use nevpt2wfn
      implicit none

#include "mxdm.fh"
#include "caspt2.fh"

      call nevpt2wfn_close
      call deinit_energies
      call finalize_inforb_molcas

      if(allocated(MultGroup%State))
     &  deallocate(MultGroup%State)
      if(allocated(MultGroup%h5_file_name))
     &  deallocate(MultGroup%h5_file_name)

      !> close file LUONEM
      call DaClos(LUONEM)

      end subroutine pt2close
