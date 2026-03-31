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
! Copyright (C) 1987, Per Ake Malmqvist                                *
!***********************************************************************

subroutine MKSXY(CMO1,CMO2,SXY)
! PURPOSE: FORM THE OVERLAP MATRIX SXY FOR ORBITAL BASES CMO1, CMO2.
! CODED 1987-02-18, P-AA M.

use OneDat, only: sNoNuc, sNoOri
use Symmetry_Info, only: nIrrep
use rassi_data, only: NBASF, NCMO, NOSH, NSXY
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO1(NCMO), CMO2(NCMO), SXY(NSXY)
integer(kind=iwp) :: ICMO, ICMP, IOPT, IRC, ISXY, ISY, ISYLAB, LSZZ1, NB, NO, NPROD, NSSQ, NSZZ
character(len=8) :: LABEL
real(kind=wp), allocatable :: PROD(:), SSQ(:), SZZ(:)

! CALCULATE SIZE AND ALLOCATE A FIELD SZZ FOR OVERLAP MATRIX
! IN COMMON BASIS SET (TRIANGULAR), SSQ TEMPORARY STORAGE
! FOR EACH OF ITS SYMMETRY BLOCKS (SQUARE), AND PROD FOR
! INTERMEDIATE MATRIX PRODUCTS.
NSZZ = 0
NSSQ = 0
NPROD = 0
do ISY=1,nIrrep
  NO = NOSH(ISY)
  NB = NBASF(ISY)
  NSZZ = NSZZ+(NB*(NB+1))/2
  NSSQ = max(NSSQ,NB**2)
  NPROD = max(NPROD,NO*NB)
end do
call mma_allocate(SZZ,NSZZ,Label='SZZ')
call mma_allocate(SSQ,NSSQ,Label='SSQ')
call mma_allocate(PROD,NPROD,Label='PROD')
! READ OVERLAP MATRIX SZZ:
IRC = -1
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICMP = 1
ISYLAB = 1
LABEL = 'MLTPL  0'
call RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '      *** ERROR IN SUBROUTINE MKSXY ***'
  write(u6,*) '     OVERLAP INTEGRALS ARE NOT AVAILABLE'
  write(u6,*)
  call ABEND()
end if
!  LOOP OVER SYMMETRIES:
LSZZ1 = 1
ISXY = 1
ICMO = 1
do ISY=1,nIrrep
  NB = NBASF(ISY)
  if (NB == 0) cycle
  NO = NOSH(ISY)
  if (NO /= 0) then
    call SQUARE(SZZ(LSZZ1),SSQ,1,NB,NB)
    ! PROD:=SSQ*CMO2
    call DGEMM_('N','N',NB,NO,NB,One,SSQ,NB,CMO2(ICMO),NB,Zero,PROD,NB)
    ! SXY:=(CMO1(TRANSP))*PROD
    call DGEMM_('T','N',NO,NO,NB,One,CMO1(ICMO),NB,PROD,NB,Zero,SXY(ISXY),NO)
    ISXY = ISXY+NO**2
    ICMO = ICMO+NO*NB
  end if
  LSZZ1 = LSZZ1+(NB*(NB+1))/2
end do
call mma_deallocate(SZZ)
call mma_deallocate(SSQ)
call mma_deallocate(PROD)

end subroutine MKSXY
