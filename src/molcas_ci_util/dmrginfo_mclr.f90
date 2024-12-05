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
Module DMRGInfo
#include "dmrginfo_mclr.fh"
save

#ifdef _NOT_USED_
!     !> These are used for the MCLR part
logical :: doMCLR
logical :: doDMRG

integer :: ndets_RGLR          ! number of Slater deternimants
integer :: ncsfs_RGLR          ! number of CSFs

integer :: nele_RGLR
integer :: MS2_RGLR

integer :: LRras2(1:20)        ! CI space when solving LR equation
integer :: RGras2(1:20)        ! CI space when solving LR equation

integer :: nstates_RGLR        ! number of states
#endif

End Module DMRGInfo

