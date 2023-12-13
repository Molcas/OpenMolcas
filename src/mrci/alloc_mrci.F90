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

subroutine ALLOC_MRCI()

use mrci_global, only: IFIRST, IPASS, IROW, ISMAX, JSC, KBUFF1, LN, MBUF, MEMWRK, MXVEC, NBITM1, NBITM2, NBITM3, NBMN, NCHN1, &
                       NCHN2, NCHN3, NREF, NRROOT, NSECT, NSYM, NVIR, NVIRT, NVMAX, NVSQ
use guga_util_global, only: IAD10
use Definitions, only: wp, iwp, RtoI

implicit none
#include "warnings.h"
integer(kind=iwp) :: I, ILIM, MBUF1, MBUF2, MEMB, MEMX, NARR, NBSIZ1, NBSIZ2, NBSIZ3, NBUFBI, NHREF, NIJ, NIJKL, NOT2, NOVLY1, &
                     NPER, NPLEN, NVT

ILIM = 4
if (IFIRST /= 0) ILIM = 2
NVSQ = 0
NVMAX = 0
do I=1,NSYM
  NVMAX = max(NVMAX,NVIR(I))
  NVSQ = NVSQ+NVIR(I)**2
end do
NVT = (NVIRT*(NVIRT+1))/2
if (NVIRT == 0) then
  call SysAbendMsg('alloc_mrci.f:','no virtual orbitals in the basis',' ')
end if
!PAM04 MOVE ALLOCATION OF FOCK MATRIX TO SDCI.
! FOCK MATRIX IN SORT,IIJJ,IJIJ,FIJ,AI.
!-----------------------------------------------
!PAM97 ! BUFFER FOR MOTRA INTEGRALS, TIBUF, IN SORT, SORTA, SORTB.
!PAM06: Evidently some miscounting margin needed:
MEMX = int(0.9_wp*real(MEMWRK,kind=wp))
! ALLOCATION FOR SORTA.
NCHN1 = LN*NVIRT+1
if (IFIRST /= 0) NCHN1 = 1
NBSIZ1 = MEMX/NCHN1-1
NBUFBI = KBUFF1
NBSIZ1 = min(NBSIZ1,MEMX-2*ISMAX-NBUFBI-1)
NBSIZ1 = max(NBSIZ1,256)
!PAM96 NBSIZ1=Size counted in real*8 words.
!PAM96 Must contain NBITM1 real*8 + NBITM1 integ + 2 integ:
NBITM1 = (RTOI*NBSIZ1-2)/(RTOI+1)
NBITM1 = min(NBITM1,NVSQ)
NBITM1 = ((NBITM1+2)/RTOI)*RTOI-2
NBSIZ1 = ((RTOI+1)*NBITM1+2+(RTOI-1))/RTOI
! BIAC
! DYNAMIC ALLOCATION FOR SORTING ABCD
NBITM2 = 1
if (IFIRST == 0) then
  IPASS = 0
  do
    IPASS = IPASS+1
    NCHN2 = (NVT-1)/IPASS+1
    NBSIZ2 = (MEMX-2*ISMAX-KBUFF1)/NCHN2
    if (RTOI*NBSIZ2 > ((RTOI+1)*NVSQ+2)) exit
    if (IPASS == 5) exit
    if (NBSIZ2 >= 1024) exit
  end do
  NBITM2 = (RTOI*NBSIZ2-2)/(RTOI+1)
  NBITM2 = min(NBITM2,NVSQ)
  NBITM2 = ((NBITM2+2)/RTOI)*RTOI-2
  NBSIZ2 = ((RTOI+1)*NBITM2+2+(RTOI-1))/RTOI
end if
! DYNAMIC ALLOCATION FOR SORTING AIBJ, AND CREATING HDIAG:
NOT2 = IROW(LN+1)
NCHN3 = 3*NOT2
NBSIZ3 = (MEMWRK-1)/NCHN3
NBSIZ3 = max(NBSIZ3,256)
NBITM3 = (RTOI*NBSIZ3-2)/(RTOI+1)
NBITM3 = min(NBITM3,NVSQ)
NBITM3 = ((NBITM3+2)/RTOI)*RTOI-2
NBSIZ3 = ((RTOI+1)*NBITM3+2+(RTOI-1))/RTOI
! CALCULATE HOW MUCH SCRATCH WILL BE NEEDED FOR PERS PART:
! FIRST, SET ASIDE WHATS NEEDED FOR SIGMA GENERATION:
NIJ = (LN*(LN+1))/2
NIJKL = NIJ*(NIJ+1)/2
NBMN = IAD10(1)
if (IFIRST /= 0) NBMN = 0
NPER = 5*NVSQ+NBSIZ3+2*NVMAX**2
NPER = max(NPER,NIJKL)
NPER = max(NPER,2*NBMN+2*ISMAX+KBUFF1)
NPER = max(NPER,NBSIZ3+2*NVMAX**2+2*NVSQ)
! OVERLAY CI,((HREF,PLEN) & (SGM,PERS PART))
NHREF = (NREF*(NREF+1))/2
NPLEN = NREF
NOVLY1 = JSC(ILIM)+max(JSC(ILIM)+NPER,NHREF+NPLEN)
! THIS IS TO BE OVERLAYED WITH (CBUF,...,LSCR) IN MQCT. TWO ALT:
NARR = 11*NRROOT**2
MEMB = MEMWRK
MBUF1 = MEMB-NOVLY1-NARR
MBUF2 = (MEMB-NARR-(3*NRROOT+2*MXVEC)*NSECT)/(3*MXVEC+2)
MBUF = min(MBUF1,MBUF2,20249)
MBUF = max(MBUF,1259)
! ALLOCATION OF PERS PART ENDS HERE ------------------------------

end subroutine ALLOC_MRCI
