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

subroutine SG_Setup_RASSCF()

use general_data, only: iSpin, Level, nActel, nConf, NLEV, nRas, nRasEl, nRsPrt, NSM, nSym, STSYM
use CI_States, only: CI_Setup_SG
use CI_Info, only: CI_St
#ifdef _DMRG_
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: dum1, dum2, dum3, Eterna_1, Eterna_2
integer(kind=iwp), parameter :: istate = 1

! Construct the Guga tables

call Timing(Eterna_1,dum1,dum2,dum3)

call CI_Setup_SG(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,Level,Level,NLEV,NSM,'FULL')

call SETSXCI()

NCONF = CI_St(istate)%nCSF(STSYM)

call Timing(Eterna_2,dum1,dum2,dum3)

end subroutine sg_setup_rasscf
