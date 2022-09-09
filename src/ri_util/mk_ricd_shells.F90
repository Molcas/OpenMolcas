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
! Copyright (C) 2007,2008, Roland Lindh                                *
!***********************************************************************

subroutine Mk_RICD_Shells()
!***********************************************************************
!                                                                      *
!    Objective: To generate aCD auxiliary basis sets on-the-fly.       *
!                                                                      *
! Called from: RdCtl_Seward                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chem. Phys., Lund Univ., Sweden.  *
!                                                                      *
!             Final implementation for aCD and acCD auxiliary          *
!             basis sets developed while visiting N. Ferre' at the     *
!             Univ. of Provance (champus Univ. Paul Cezanne) in        *
!             Marseille, France, 20 March - 19 April.                  *
!                                                                      *
!             Modified to transform the auxiliary basis to a true      *
!             Cholesky basis set while on TACC 2008 conference in      *
!             Songjiang District, Shanghai, China, 23-27 Sept. 2008.   *
!                                                                      *
!***********************************************************************

use Real_Spherical, only: Sphere, Sphere_Free
use Basis_Info, only: dbsc, nCnttp
use Sizes_of_Seward, only: S
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iCnttp, jCnttp, mCnttp, nDiff
logical(kind=iwp) :: DoRys, W2L

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iPrint = 49
!iPrint = 99
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

call StatusLine('Gateway:',' Generating aCD or acCD auxiliary basis set')
!                                                                      *
!***********************************************************************
!                                                                      *
! Preamble: Compute kOffAO  and lOffAO

call Setup_OffAO()

! Set up transformation matrix from Cartesian to real spherical harmonics.

call Sphere(S%iAngMx)

! Setup of tables for coefficients for the Rys roots and weights.

nDiff = 0
if (S%iAngMx == 0) nDiff = 2
DoRys = .true.
call SetUp_RW(DoRys,nDiff)

mCnttp = nCnttp
!                                                                      *
!***********************************************************************
!                                                                      *
! Add the DUMMY SHELL!

call Mk_Dummy_Shell()
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Loop now over all unique valence basis sets and generate the
! corresponding aCD auxiliary basis sets. Note that there are two
! different types of aCD auxiliary basis sets, aCD and acCD.

do iCnttp=1,mCnttp
  if (dbsc(iCnttp)%Frag .or. (dbsc(iCnttp)%nVal == 0)) cycle
# ifdef _DEBUGPRINT_
  if (iPrint >= 99) write(u6,*) 'Generating auxiliary basis set for valence basis:',iCnttp
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Procrastinate the printing of the RICD basis set to library
  ! until the last unique valence basis set is processed.

  W2L = .true.
  do jCnttp=iCnttp+1,mCnttp
    if (dbsc(iCnttp)%Bsl_old == dbsc(jCnttp)%Bsl_old) then
      W2L = .false.
      exit
    end if
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! aCD and acCD section

  call Mk_aCD_acCD_Shells(iCnttp,W2L)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
!                                                                      *
!***********************************************************************
!                                                                      *
! Cleanup the mess!

call CloseR()
call Sphere_Free()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_RICD_Shells
