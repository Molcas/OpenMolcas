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
Module Jobiph_J
use Definitions, only: wp, iwp
Implicit None
#include "Molcas.fh"
integer(kind=iwp) iadr15_j(15)
integer nactel_j,ispin_j,nsym_j,lsym_j,nfro_j(mxsym),nish_j(mxsym),nash_j(mxsym),ndel_j(mxsym),nbas_j(mxsym)
character(Len=lenin8) name_j(mxorb)
integer(kind=iwp) nconf_j
character(LEN=2) header_j(72)
character(LEN=72) title_j(18)
Real(kind=wp) potnuc_j
integer(kind=iwp) lroots_j,nroots_j,iroot_j(mxroot),nrs1_j(mxsym),nrs2_j(mxsym),nrs3_j(mxsym),nhole1_j,nelec3_j,ipt2_j
Real(kind=wp) weight_j(mxroot)

Private
Public:: iadr15_j,nactel_j,ispin_j,nsym_j,lsym_j,nfro_j,nish_j,nash_j,ndel_j,nbas_j,name_j,nconf_j,header_j, title_j, &
         potnuc_j,lroots_j,nroots_j,iroot_j,nrs1_j,nrs2_j,nrs3_j,nhole1_j,nelec3_j,ipt2_j,weight_j

End Module Jobiph_J
