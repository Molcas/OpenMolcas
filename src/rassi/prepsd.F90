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

subroutine PREPSD(WFTP,SGS,CIS,LSYM,ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,NCONF,CI,DET,detocc,detcoeff,SPTRA)
! Purpose: Given a RASSCF wave function in Split-GUGA format
! and an orbital transformation matrix for the purpose of
! getting biorthonormal orbitals, prepare a wave function
! in the general SD format, using transformed orbitals.

use gugx, only: CIStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
character(len=8) :: WFTP
type(SGStruct) :: SGS
type(CIStruct) :: CIS
integer(kind=iwp) :: LSYM, ICNFTAB(*), ISPNTAB(*), ISSTAB(*), IFSBTAB(*), NCONF
real(kind=wp) :: CI(*), DET(*)
character(len=*), intent(out) :: detocc(NCONF)
real(kind=wp), intent(out) :: detcoeff(NCONF)
real(kind=wp), intent(in) :: SPTRA(*)
integer(kind=iwp) :: IMODE
real(kind=wp), allocatable :: CTMP(:)

if (WFTP == 'GENERAL') then
  ! Transform SGUGA to SymmG:
  call mma_allocate(CTMP,NCONF,Label='CTMP')
  IMODE = 1
  call SYG2SGU(IMODE,SGS,CIS,LSYM,ICNFTAB,ISPNTAB,CI,CTMP)
  ! Transform SymmG to Slater Dets:
  call SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,CTMP,DET,detocc,detcoeff,SPTRA)
  call mma_deallocate(CTMP)
else
  DET(1) = CI(1)
end if

end subroutine PREPSD
