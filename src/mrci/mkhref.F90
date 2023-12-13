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

subroutine MKHREF(HREF,FC,FIJKL,JREFX)

use mrci_global, only: IAD25S, INDSRT, IRC, IROW, IVVER, LASTAD, LN, Lu_25, Lu_70, LUSYMB, NBTRI, NCVAL, NREF, NSRTMX, POTNUC, &
                       VALSRT
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: HREF(*), FC(*), FIJKL(*)
integer(kind=iwp), intent(in) :: JREFX(NCVAL)
integer(kind=iwp) :: i, IADD10, IADD25, IADR, IBUF, IC1, IC2, ICHK, IIN, IIR, IK, ILEN, IND, INDI, IP, IR, IVEC, IVL, JP, KP, &
                     LENGTH, LP, NA, NAT, NB, NHREF, NI, NIJ, NIJKL, NK, NKL
real(kind=wp) :: FINI

NHREF = (NREF*(NREF+1))/2
HREF(1:NHREF) = Zero
ICHK = 0
IK = 0
FINI = Zero
IADD25 = 0
call dDAFILE(Lu_25,2,FC,NBTRI,IADD25)
IADD10 = IAD10(8)
do
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN < 0) exit
  do IIN=1,ILEN
    IND = ICOP1(IIN)
    if (ICHK /= 0) then
      ICHK = 0
      INDI = IND
      NI = ibits(INDI,0,10)
      NK = ibits(INDI,10,10)
      IK = IROW(NK)+NI
      FINI = FC(IK)
    else if (IND == 0) then
      ICHK = 1
    else
      IVL = ibits(IND,0,6)
      if (IVL /= IVVER) cycle
      IC2 = ibits(IND,6,13)
      NA = JREFX(IC2)
      if (NA == 0) cycle
      IC1 = ibits(IND,19,13)
      NB = JREFX(IC1)
      if (NB == 0) cycle
      if (NA < NB) then
        NAT = NA
        NA = NB
        NB = NAT
      end if
      IVEC = (NA*(NA-1))/2+NB
      HREF(IVEC) = HREF(IVEC)+COP(IIN)*FINI
    end if
  end do
end do
ICHK = 0
NIJ = IROW(LN+1)
NIJKL = NIJ*(NIJ+1)/2
FIJKL(1:NIJKL) = Zero
FINI = Zero
IADR = LASTAD(1)
do
  call dDAFILE(Lu_70,2,VALSRT,NSRTMX,IADR)
  call iDAFILE(Lu_70,2,INDSRT,NSRTMX+2,IADR)
  LENGTH = INDSRT(NSRTMX+1)
  IADR = INDSRT(NSRTMX+2)
  do i=1,length
    FIJKL(INDSRT(i)) = VALSRT(i)
  end do
  if (IADR == -1) exit
end do
IADD10 = IAD10(5)
do
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN < 0) exit
  do IIN=1,ILEN
    IND = ICOP1(IIN)
    if (ICHK /= 0) then
      ICHK = 0
      INDI = IND
      IP = ibits(INDI,0,8)
      JP = ibits(INDI,8,8)
      KP = ibits(INDI,16,8)
      LP = ibits(INDI,24,8)
      NIJ = IROW(IP)+JP
      NKL = IROW(KP)+LP
      IND = NIJ*(NIJ-1)/2+NKL
      FINI = FIJKL(IND)
    else if (IND == 0) then
      ICHK = 1
    else
      IVL = ibits(IND,0,6)
      if (IVL /= 0) cycle
      IC2 = ibits(IND,6,13)
      NA = JREFX(IC2)
      if (NA == 0) cycle
      IC1 = ibits(IND,19,13)
      NB = JREFX(IC1)
      if (NB == 0) cycle
      if (NA < NB) then
        NAT = NA
        NA = NB
        NB = NAT
      end if
      IVEC = (NA*(NA-1))/2+NB
      HREF(IVEC) = HREF(IVEC)+COP(IIN)*FINI
    end if
  end do
end do
IADD25 = IAD25S
IBUF = nCOP
do I=1,IRC(1)
  IBUF = IBUF+1
  if (IBUF > nCOP) then
    call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
    IBUF = 1
  end if
  IR = JREFX(I)
  if (IR > 0) then
    IIR = (IR*(IR+1))/2
    HREF(IIR) = HREF(IIR)+COP(IBUF)+POTNUC
  end if
end do

return

end subroutine MKHREF
