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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine NPSET(JSY,INDEX,C,TPQ,ENP,T,S,W,EPP,ICASE)

implicit real*8(A-H,O-Z)
dimension JSY(*), index(*), C(*), TPQ(*), ENP(*), T(*), S(*)
dimension W(*), EPP(*), ICASE(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
! Statement function
JSYM(L) = JSUNP_CPF(JSY,L)

if (IDENS == 1) GO TO 65
if (ITPUL /= 1) GO TO 60
IAD = 0
IADDP(1) = 0
call dDAFILE(Lu_CI,1,C,NCONF,IAD)
IADDP(2) = IAD

! VALENCE

60 IQ = IRC(1)
do I=1,IQ
  T(I) = C(I)*C(I)
end do

! SINGLES

IQ = IRC(2)-IRC(1)
IND = IRC(1)
IN = IRC(1)
do I=1,IQ
  IND = IND+1
  NS1 = JSYM(IN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = index(IN+I)+1
  T(IND) = DDOT_(INUM,C(IST),1,C(IST),1)
end do

! DOUBLES

IQ = IRC(4)-IRC(2)
IN = IRC(2)
do I=1,IQ
  IND = IND+1
  NS1 = JSYM(IN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = index(IN+I)+1
  T(IND) = DDOT_(INUM,C(IST),1,C(IST),1)
end do
IP = IRC(4)
do I=1,IP
  call TPQSET(ICASE,TPQ,I)
  ENP(I) = DDOT_(IP,TPQ,1,T,1)
  ENP(I) = ENP(I)+D1
end do
IP = IRC(4)
if (IPRINT > 5) write(6,12) (ENP(I),I=1,IP)
12 format(6X,'ENP ',5F14.8)

! VALENCE

65 IQ = IRC(1)
do I=1,IQ
  if (IDENS == 0) EMPI = D1/sqrt(ENP(I))
  if (IDENS == 1) EMPI = sqrt(ENP(I))
  C(I) = C(I)*EMPI
end do

! SINGLES

IQ = IRC(2)-IRC(1)
IN = IRC(1)
do I=1,IQ
  NS1 = JSYM(IN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = index(IN+I)+1
  if (IDENS == 0) EMPI = D1/sqrt(ENP(IN+I))
  if (IDENS == 1) EMPI = sqrt(ENP(IN+I))
  call VSMUL(C(IST),1,EMPI,C(IST),1,INUM)
end do

! DOUBLES

IQ = IRC(4)-IRC(2)
IN = IRC(2)
do I=1,IQ
  NS1 = JSYM(IN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = index(IN+I)+1
  if (IDENS == 0) EMPI = D1/sqrt(ENP(IN+I))
  if (IDENS == 1) EMPI = sqrt(ENP(IN+I))
  call VSMUL(C(IST),1,EMPI,C(IST),1,INUM)
end do
if (IPRINT >= 15) write(6,13) (C(I),I=1,NCONF)
13 format(6X,'C(NP)',5F10.6)
if (IDENS == 1) return

call SETZ(EPP,IRC(4))
call SETZ(S,JSC(4))
if ((ICPF /= 1) .and. (ISDCI /= 1) .and. (INCPF /= 1)) call SETZ(W,JSC(4))

return

end subroutine NPSET
