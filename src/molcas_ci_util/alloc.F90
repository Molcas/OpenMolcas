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

subroutine ALLOC()
! RASSCF: allocation of core memory
!
! Called from inpctl
!
! No subroutine calls
!
! ********** IBM-3090 Release 88 10 11 **********

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use rasscf_global, only: ISTORD, ISTORP, nFint
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
#include "rasdim.fh"
#include "general.fh"
integer(kind=iwp) :: IORD, IORP, NAP, NAQ, NAR, NAS, NOP, NRS, NSP, NSPQ, NSPQR, NSQ, NSR, NSS

#ifdef _DEBUGPRINT_
write(u6,*) ' Entering ALLOC'
#endif

! Compute space needed for transformed two-electron integrals

ISTORD(1) = 0
ISTORP(1) = 0
IORD = 0
IORP = 0
do NSP=1,NSYM
  NOP = NORB(NSP)
  NAP = NASH(NSP)
  do NSQ=1,NSYM
    NAQ = NASH(NSQ)
    NSPQ = Mul(NSP,NSQ)
    do NSR=1,NSYM
      NSPQR = Mul(NSPQ,NSR)
      NAR = NASH(NSR)
      do NSS=1,NSR
        if (NSPQR /= NSS) cycle
        NAS = NASH(NSS)
        NRS = NAR*NAS
        if (NSS == NSR) NRS = nTri_Elem(NAR)
        IORD = IORD+NOP*NAQ*NRS
        IORP = IORP+NAP*NAQ*NRS
      end do
    end do
  end do
  ISTORD(NSP+1) = IORD
  ISTORP(NSP+1) = IORP
end do
NFINT = ISTORD(NSYM+1)

#ifdef _DEBUGPRINT_
write(u6,'(1X,A,5X,9I5)') 'ISTORD-vector:',ISTORD(1:NSYM+1)
#endif

end subroutine ALLOC
