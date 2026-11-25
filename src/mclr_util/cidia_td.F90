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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine CIDIA_TD(iSym)

use ipPage, only: ipget, ipin, ipnout, W
use CI_Lib_Defs, only: CI_Diag_Cnf_Avg
use MCLR_procedures, only: MCLR_CI_Diag
use MCLR_Data, only: FIMO, INT2, ipDia, nexp, nexp_max
use CI_Lib_Info, only: CI_Lib_Space
use input_mclr, only: ERASSCF, nCSF, PotNuc, rIn_Ene, TimeDep
use Definitions, only: iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iSym
integer(kind=iwp) :: ipDIAI, nD, nP1, nP2, nQ, nSD

! This is just a interface to hide Jeppe from the rest of the world
! we don't want to let world see the work of the Danish
! (I hope he never reads that)
! Anyway concerning the CSF/SD stuff.
! If we work with spin dependent perturbations
! we never use CSF's (to complicated), instead we use
! SD in all parts of the program,
! otherwise we will switch to SD representation in this routine

! The CSF diagonal through the CI interface (MCLR_CI_Diag -> CI_Lib_Diag)
nsd = max(ncsf(isym),CI_Lib_Space(iState)%nDet(ISYM))
ipdiai = ipget(nsd)
call ipin(ipdiai)
call MCLR_CI_Diag(CI_Diag_Cnf_Avg,iSym,size(FIMO),FIMO,size(INT2),INT2,rin_ene+potnuc,ncsf(iSym),W(ipdiai)%A(1:ncsf(iSym)))
nD = nCSF(iSym)
! Calculate explicit part of hamiltonian

np2 = min(nd,nexp_max)
np1 = 0
nq = 0
if (np2 /= 0) then
  call ipnout(ipdiai)
  call ipin(ipdiai)
  call h0(W(ipdiai)%A,np1,nexp_max,nq,isym,nexp,TimeDep)
else
  nexp = 0
end if

call ipin(ipdiai)
W(ipdiai)%A(1:nD) = W(ipdiai)%A(1:nD)-ERASSCF(1)

ipdia = ipdiai

end subroutine CIDIA_TD
