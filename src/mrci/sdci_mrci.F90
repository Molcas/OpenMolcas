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

use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "mrci.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: NHREF, NIJ, NIJKL
!real(kind=wp) :: H(MAXMEM)
!integer(kind=iwp) :: iH(RtoI*MAXMEM)

! PUT THE SUBROUTINE NAME ONTO THE ENTRY NAME STACK
! INPUT AND MEMORY ALLOCATION:
!call READIN(HWork,iHWork)
call READIN_MRCI()
! INTEGRAL SORTING AND DIAGONAL ELEMENTS:
! USE COUPLING COEFFS FROM UNIT 10, TRANSFORMED INTEGRALS FROM 13 AND 17
! PRODUCE FILES UNIT 14, 15 AND 16 WITH SORTED INTEGRALS.
! ALSO FOCK MATRIX TO UNIT 25 AND DIAGONAL ELEMENTS TO UNIT 27.
!PAM04 ALLOCATION OF FOCK MATRIX MOVED HERE FROM ALLOC.
call GETMEM('FOCK','ALLO','REAL',LFOCK,NBTRI)
call DIAGCT()
! CREATE REFERENCE CI HAMILTONIAN:
!PAM04 call MKHREF(HWork(LHREF),Hwork(LFOCK),HWork(LFIJKL),HWork(LJREFX))
NHREF = (NREF*(NREF+1))/2
call GETMEM('HREF','ALLO','REAL',LHREF,NHREF)
NIJ = (LN*(LN+1))/2
NIJKL = (NIJ*(NIJ+1))/2
call GETMEM('FIJKL','ALLO','REAL',LFIJKL,NIJKL)
call MKHREF(Work(LHREF),Work(LFOCK),Work(LFIJKL),IWork(LJREFX))
! SOLVE REFERENCE CI EQUATIONS:
!PAM04 call REFCI(HWork(LHREF),HWork(LAREF),HWork(LEREF),HWork(LCSPCK),HWork(LCISEL),HWork(LPLEN))
call GETMEM('AREF','ALLO','REAL',LAREF,NREF**2)
call GETMEM('EREF','ALLO','REAL',LEREF,NREF)
call GETMEM('PLEN','ALLO','REAL',LPLEN,NREF)
call REFCI(Work(LHREF),Work(LAREF),Work(LEREF),IWork(LCSPCK),Work(LCISEL),Work(LPLEN))
call GETMEM('PLEN','FREE','REAL',LPLEN,NREF)
call GETMEM('HREF','FREE','REAL',LHREF,NHREF)
if (IREFCI == 1) then
  call GETMEM('FOCK','FREE','REAL',LFOCK,NBTRI)
  call GETMEM('FIJKL','FREE','REAL',LFIJKL,NIJKL)
else
  ! SOLVE MRCI OR ACPF EQUATIONS:
  ! FIRST, SET UP START CI ARRAYS, AND ALSO TRANSFORM DIAGONAL ELEMENTS:
  !------
  ! POW: Initialize HSMALL(1,1)
  HSMALL(1,1) = Zero
  !------
  call GETMEM('ICI','ALLO','INTE',LICI,MBUF)
  call GETMEM('CI','ALLO','REAL',LCI,NCONF)
  call GETMEM('SGM','ALLO','REAL',LSGM,NCONF)
  call CSTART(Work(LAREF),Work(LEREF),Work(LCI),IWork(LICI))
  call MQCT(WORK(LAREF),WORK(LEREF),Work(LCI),Work(LSGM),IWork(LICI))
  call GETMEM('SGM','FREE','REAL',LSGM,NCONF)
  call GETMEM('CI','FREE','REAL',LCI,NCONF)
  call GETMEM('ICI','FREE','INTE',LICI,MBUF)
  !PAM04 EXPLICIT DEALLOCATION OF FOCK MATRIX
  call GETMEM('FOCK','FREE','REAL',LFOCK,NBTRI)
  ! DENSITY (AND MAYBE TRANSITION DENSITY) MATRICES IN AO BASIS:
  !PAM04 ALLOCATION OF DMO AND TDMO MOVED HERE FROM ALLOC:
  call GETMEM('DMO','ALLO','REAL',LDMO,NBTRI)
  if (ITRANS == 1) call GETMEM('TDMO','ALLO','REAL',LTDMO,NBAST**2)
  !PAM04 End of addition
  call DENSCT(WORK(LAREF))
  call GETMEM('AREF','FREE','REAL',LAREF,NREF**2)
  call GETMEM('EREF','FREE','REAL',LEREF,NREF)
  ! NATURAL ORBITALS AND PROPERTIES (AND MAYBE TRANSITION PROPS):
  call PROPCT()
  !PAM04 EXPLICIT DEALLOCATION ADDED:
  call GETMEM('DMO','FREE','REAL',LDMO,NBTRI)
  if (ITRANS == 1) call GETMEM('TDMO','FREE','REAL',LTDMO,NBAST**2)
  !PAM04 End of addition
end if
call GETMEM('FIJKL','FREE','REAL',LFIJKL,NIJKL)
call GETMEM('CISEL','FREE','REAL',LCISEL,NSEL*NREF)
call GETMEM('JREFX','FREE','INTE',LJREFX,NCVAL)
call GETMEM('ISAB','FREE','INTE',LISAB,NVIRT**2)
call GETMEM('INDX','FREE','INTE',LINDX,NIWLK)
call GETMEM('INTSY','FREE','INTE',LINTSY,NINTSY)
call GETMEM('CSPCK','FREE','INTE',LCSPCK,NCSPCK)

return

end subroutine SDCI_MRCI
