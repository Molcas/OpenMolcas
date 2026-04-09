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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Rd1Int()
!***********************************************************************
!                                                                      *
!     Read header and matrices from the one-electron integral file     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use rasscf_global, only: BName, header, PotNuc
use general_data, only: NBAS, NSYM
use Molcas, only: LenIn
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nBas_tot

!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
!---  read file header  -----------------------------------------------*
call Get_cArray('Seward Title',Header,144)
!---  read number of symm. species ------------------------------------*
call Get_iScalar('nSym',nSym)
!---  read number of basis functions per symmetry species -------------*
call Get_iArray('nBas',nBas,nSym)
!---  read nuclear potential ------------------------------------------*
! Get POTNUC from the runfile, where it was set by seward.
! (Do not trust reading it from JOBIPH).
call Get_dScalar('potNuc',PotNuc)
!---  read basis function labels --------------------------------------*
nBas_tot = sum(nBas(1:nSym))
call Get_cArray('Unique Basis Names',BName,(LenIn+8)*nBas_tot)
!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*

end subroutine Rd1Int
