!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,2020, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine Seward_Init()
!***********************************************************************
!                                                                      *
!     Object: to set data which is stored in common blocks             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose      *
!             January 1990                                             *
!***********************************************************************

use EFP_Module, only: lEFP, nEFP_fragments
use k2_arrays, only: XMem
use k2_structure, only: k2_processed
use Basis_Info, only: Seward_Activated
use RICD_Info, only: Do_RI, iRI_Type
use rmat, only: bParm, dipol, Dipol1, Epsabs, Epsq, Epsrel, keyr, lgamma, nagint, qCoul, Quadpack, RMat_On, RmatR, testint
use DCR_mod, only: DCR_Init
use NAC, only: isCSF, isNAC
use NDDO, only: twoel_NDDO
use PrintLevel, only: nPrint, Show
use Constants, only: Zero, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iPL
character(len=180) :: Env
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
! Info

Seward_Activated = .false.

! LInfo

call GetEnvF('MOLCAS_NEW_DEFAULTS',Env)
call UpCase(Env)
if (Env == 'YES') then
  Do_RI = .true.
  iRI_Type = 4
end if

iPL = iPrintLevel(-1)
if (iPL == 2) then
  iPL = 5
else if (iPL == 3) then
  iPL = 6
else if (iPL == 4) then
  iPL = 7
else if (iPL == 5) then
  iPL = 49  ! 99 would be just too much
end if
nPrint(:) = iPL
if ((Reduce_Prt() .and. (iPL < 6)) .or. (iPL == 0)) then
  Show = .false.
else
  Show = .true.
end if

twoel_NDDO = .false.

Seward_Activated = .true.
XMem = .false.
k2_processed = .false.

call Set_Binom()
call Set_CanInd()

! Set some default value for RMAT type integration

! rmat.fh

RmatR = Ten
Epsabs = 1.0e-9_wp
Epsrel = 1.0e-14_wp
qCoul = Zero
Epsq = 1.0e-8_wp
bParm = Zero
dipol(1) = Zero
dipol(2) = Zero
dipol(3) = Zero
Dipol1 = Zero
keyr = 6
Quadpack = .true.
nagint = .false.
testint = .false.
RMat_On = .false.
lgamma = 9

call DCR_Init()

call Set_Basis_Mode('Valence')

! nac.fh

isNAC = .false.
isCSF = .false.

! EFP stuff

lEFP = .false.
nEFP_fragments = 0

end subroutine Seward_Init
