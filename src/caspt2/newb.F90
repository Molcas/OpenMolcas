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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!***********************************************************************
! This file contains boiler-plate code to allow developers to experiment
! with modifications to the zero-order hamiltonian. To enable them,
! add appropriate code, and use 'HZero = Custom' in the input.
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine NEWB()

use Index_Functions, only: nTri_Elem
use EQSOLV, only: iDBMat, iDSMat
use caspt2_global, only: LUSBT
use caspt2_module, only: nASup, nISup, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ICASE, IDB, IDS, ISYM, NAS, NB, NCOEF, NIS, NS
real(kind=wp), allocatable :: S(:), B(:)

! Modify B matrices, if requested.

do ICASE=1,11
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NCOEF = NAS*NIS
    if (NCOEF == 0) cycle
    NS = nTri_Elem(NAS)
    call mma_allocate(S,NS,LABEL='S')
    IDS = IDSMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,S,NS,IDS)
    NB = NS
    call mma_allocate(B,NB,LABEL='B')
    IDB = IDBMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,B,NB,IDB)
    ! Modify B matrix, using S matrix and some other data.
    call mma_deallocate(B)
    call mma_deallocate(S)
  end do
end do

end subroutine NEWB
