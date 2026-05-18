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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------
! 1998  PER-AAKE MALMQUIST
! DEPARTMENT OF THEORETICAL CHEMISTRY
! UNIVERSITY OF LUND
! SWEDEN
!--------------------------------------------
! NOTE: This new MKRHS code produces ONLY the
! contravariant components. 980928, P-A Malmqvist
!--------------------------------------------

subroutine MKRHS(IVEC)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV. The coupling matrix elements from the
! root state to the 1st order interacting space are computed, as
! combinations of MO integrals.
! This is the RHS vector in contravariant representation.

use definitions, only: iwp, wp, u6
use caspt2_global, only: iPrGlb
use caspt2_global, only: FIMO
use PrintLevel, only: VERBOSE
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NASHT, NOMX

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) NERI, NFIMO
real(kind=wp), allocatable, target :: ERI(:)
real(kind=wp), pointer :: ERI0(:), ERI1(:), ERI2(:), SCR(:)

if (IPRGLB >= VERBOSE) write(u6,'(1X,A)') ' Using conventional MKRHS algorithm'

! INTEGRAL BUFFERS:
NERI = NOMX**2
call mma_allocate(ERI,3*NERI,Label='ERI')
ERI0(1:2*NERI) => ERI(1:2*NERI)
ERI1(1:NERI) => ERI(1:NERI)
ERI2(1:NERI) => ERI(NERI+1:2*NERI)
SCR(1:NERI) => ERI(2*NERI+1:3*NERI)

if (NASHT > 0) then
  NFIMO = size(FIMO)
  call MKRHSA(IVEC,FIMO,NFIMO,ERI0,2*NERI,SCR,NERI)
  call MKRHSB(IVEC,ERI0,2*NERI,SCR,NERI)
  call MKRHSC(IVEC,FIMO,NFIMO,ERI0,2*NERI,SCR,NERI)
  call MKRHSD(IVEC,FIMO,NFIMO,ERI1,NERI,ERI2,NERI,SCR,NERI)
  call MKRHSE(IVEC,ERI1,NERI,ERI2,NERI,SCR,NERI)
  call MKRHSF(IVEC,ERI1,NERI,ERI2,NERI,SCR,NERI)
  call MKRHSG(IVEC,ERI1,NERI,ERI2,NERI,SCR,NERI)
end if
call MKRHSH(IVEC,ERI1,NERI,ERI2,NERI,SCR,NERI)

ERI0 => null()
ERI1 => null()
ERI2 => null()
SCR => null()
call mma_deallocate(ERI)

end subroutine MKRHS
