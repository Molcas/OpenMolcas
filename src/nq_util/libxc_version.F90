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
! Copyright (C) 2022, Roland Lindh                                     *
!               2022, Susi Lehtola                                     *
!***********************************************************************

subroutine libxc_version()

use xc_f03_lib_m, only: xc_f03_reference, xc_f03_reference_doi, xc_f03_version
use Definitions, only: iwp, LibxcInt, u6

implicit none
integer(kind=LibxcInt) :: vmajor, vminor, vmicro
character(len=128) :: libxc_reference, libxc_reference_doi
logical(kind=iwp), external :: Reduce_Prt

if (Reduce_Prt()) return
! Get the data from libxc
call xc_f03_version(vmajor,vminor,vmicro)
call xc_f03_reference(libxc_reference)
call xc_f03_reference_doi(libxc_reference_doi)
! Print out the version
write(u6,'(6X,"Using Libxc version: ",I0,".",I0,".",I0)') vmajor,vminor,vmicro
! Print out the Libxc literature reference
write(u6,'(6X,"Please cite the following reference:")')
write(u6,'(6X,A," doi:",A)') trim(libxc_reference),trim(libxc_reference_doi)

end subroutine libxc_version
