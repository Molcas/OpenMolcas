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

subroutine CIDIA_sa(iSym,ralp,S)

use CI_Info, only: CI_St
use ipPage, only: ipget, ipin, W
use CI_Op_Defs, only: CI_Diag_Cnf_Avg
use MCLR_procedures, only: MCLR_CI_Diag
use MCLR_Data, only: FANCY_PRECONDITIONER, FIMO, INT2, ipCI, ipDia
use CI_Info, only: CI_St
use input_mclr, only: ERASSCF, nCSF, nRoots, PotNuc, rIn_Ene, Weight
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iSym
real(kind=wp), intent(inout) :: ralp(nRoots), S(*)
integer(kind=iwp) :: iP2, ipDIAI, J, nSD
real(kind=wp) :: rTmp

! This is just a interface to hide Jeppe from the rest of the world
! we don't want to let world see the work of the Danish
! (I hope he never reads that)
! Anyway concerning the CSF/SD stuff.
! If we work with spin dependent perturbations
! we never use CSF's (to complicated), instead we use
! SD in all parts of the program,
! otherwise we will switch to SD representation in this routine

! The CSF diagonal through the CI interface (MCLR_CI_Diag -> CI_Lib_Diag)
nsd = max(ncsf(isym),CI_St(iState)%nDet(ISYM))
ipdiai = ipget(nsd)
call ipin(ipdiai)
call MCLR_CI_Diag(CI_Diag_Cnf_Avg,iSym,size(FIMO),FIMO,size(INT2),INT2,rin_ene+potnuc,ncsf(iSym),W(ipdiai)%A(1:ncsf(iSym)))
! Calculate explicit part of hamiltonian

ipdia = ipdiai

if (FANCY_PRECONDITIONER) then
  call ipin(ipdia)
  call SA_PREC(S,W(ipdia)%A)
else
  call ipin(ipdiai)
  call ipin(ipCI)
  ! DMinvCI_SA always uses S, so it must be defined here as well. Without
  ! the "fancy" preconditioner each root is only projected against its own
  ! reference vector:  S(j,j,j) = 1/<j|(H0-E_j)^-1|j>, all other elements zero.
  S(1:nRoots**3) = Zero
  ip2 = 0
  do j=1,nroots
    rTmp = sum(W(ipCI)%A(ip2+1:ip2+ncsf(CI_St(iState)%STSym))**2/(W(ipdiai)%A(1:ncsf(CI_St(iState)%STSym))-ERASSCF(j)))
    ralp(j) = Weight(j)*rTmp
    if (rTmp /= Zero) S(j+(j-1)*nRoots+(j-1)*nRoots**2) = One/rTmp
    ip2 = ip2+ncsf(CI_St(iState)%STSym)
  end do
end if

end subroutine CIDIA_sa
