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

subroutine SDCI_MRCI()

use mrci_global, only: CISEL, CSPCK, DMO, FIJKL, FOCK, INDX, INTSY, IREFCI, ISAB, ITRANS, JREFX, LN, MBUF, NBAST, NBTRI, NCONF, &
                       NREF, TDMO
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NHREF, NIJ, NIJKL
integer(kind=iwp), allocatable :: ICI(:)
real(kind=wp), allocatable :: AREF(:,:), CI(:), EREF(:), HREF(:), PLEN(:), SGM(:)

! PUT THE SUBROUTINE NAME ONTO THE ENTRY NAME STACK
! INPUT AND MEMORY ALLOCATION:
call READIN_MRCI()
! INTEGRAL SORTING AND DIAGONAL ELEMENTS:
! USE COUPLING COEFFS FROM UNIT 10, TRANSFORMED INTEGRALS FROM 13 AND 17
! PRODUCE FILES UNIT 14, 15 AND 16 WITH SORTED INTEGRALS.
! ALSO FOCK MATRIX TO UNIT 25 AND DIAGONAL ELEMENTS TO UNIT 27.
!PAM04 ALLOCATION OF FOCK MATRIX MOVED HERE FROM ALLOC.
call mma_allocate(FOCK,NBTRI,label='FOCK')
call DIAGCT()
! CREATE REFERENCE CI HAMILTONIAN:
NHREF = (NREF*(NREF+1))/2
call mma_allocate(HREF,NHREF,label='HREF')
NIJ = (LN*(LN+1))/2
NIJKL = (NIJ*(NIJ+1))/2
call mma_allocate(FIJKL,NIJKL,label='FIJKL')
call MKHREF(HREF,FOCK,FIJKL,JREFX)
! SOLVE REFERENCE CI EQUATIONS:
call mma_allocate(AREF,NREF,NREF,label='AREF')
call mma_allocate(EREF,NREF,label='EREF')
call mma_allocate(PLEN,NREF,label='PLEN')
call REFCI(HREF,AREF,EREF,CSPCK,CISEL,PLEN)
call mma_deallocate(HREF)
call mma_deallocate(PLEN)
if (IREFCI /= 1) then
  ! SOLVE MRCI OR ACPF EQUATIONS:
  ! FIRST, SET UP START CI ARRAYS, AND ALSO TRANSFORM DIAGONAL ELEMENTS:
  !------
  call mma_allocate(ICI,MBUF,label='ICI')
  call mma_allocate(CI,NCONF,label='CI')
  call mma_allocate(SGM,NCONF,label='SGM')
  call CSTART(AREF,EREF,CI,ICI)
  call MQCT(AREF,EREF,CI,SGM,ICI)
  call mma_deallocate(CI)
  call mma_deallocate(SGM)
  call mma_deallocate(ICI)
  ! DENSITY (AND MAYBE TRANSITION DENSITY) MATRICES IN AO BASIS:
  !PAM04 ALLOCATION OF DMO AND TDMO MOVED HERE FROM ALLOC:
  call mma_allocate(DMO,NBTRI,label='DMO')
  if (ITRANS == 1) call mma_allocate(TDMO,NBAST,NBAST,label='TDMO')
  !PAM04 End of addition
  call DENSCT(AREF)
  call mma_deallocate(AREF)
  call mma_deallocate(EREF)
  ! NATURAL ORBITALS AND PROPERTIES (AND MAYBE TRANSITION PROPS):
  call PROPCT()
  !PAM04 EXPLICIT DEALLOCATION ADDED:
  call mma_deallocate(DMO)
  if (ITRANS == 1) call mma_deallocate(TDMO)
  !PAM04 End of addition
end if
call mma_deallocate(FOCK)
call mma_deallocate(FIJKL)
call mma_deallocate(CISEL)
call mma_deallocate(JREFX)
call mma_deallocate(ISAB)
call mma_deallocate(INDX)
call mma_deallocate(INTSY)
call mma_deallocate(CSPCK)

return

end subroutine SDCI_MRCI
