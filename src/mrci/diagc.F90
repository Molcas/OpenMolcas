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

subroutine DIAGC(INTSYM,C,S)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: INTSYM(*)
real(kind=wp) :: C(*), S(*)
#include "mrci.fh"
integer(kind=iwp) :: IADD25, IIC, ILIM, IND, INDA, IRL, NA, NA1, NA2, NB, NB1, NB2, NSA, NSS
integer(kind=iwp), external :: JSUNP
!Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

IADD25 = IAD25S
call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
IIC = 0
IND = 0
ILIM = 4
if (IFIRST /= 0) ILIM = 2
IRL = IRC(ILIM)
do INDA=1,IRL
  NSS = MUL(JSYM(INDA),LSYM)
  if (INDA > IRC(1)) GO TO 120
  IIC = IIC+1
  IND = IND+1
  S(IND) = S(IND)+COP(IIC)*C(IND)
  if (IIC < nCOP) GO TO 100
  call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
  IIC = 0
  GO TO 100
120 if (INDA > IRC(2)) GO TO 130
  NA1 = NVIRP(NSS)+1
  NA2 = NVIRP(NSS)+NVIR(NSS)
  if (NA2 < NA1) GO TO 100
  do NA=NA1,NA2
    IIC = IIC+1
    IND = IND+1
    S(IND) = S(IND)+COP(IIC)*C(IND)
    if (IIC < nCOP) GO TO 121
    call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
    IIC = 0
121 continue
  end do
  GO TO 100
130 do NA=1,NVIRT
    NSA = MUL(NSS,NSM(LN+NA))
    NB1 = NVIRP(NSA)+1
    NB2 = NVIRP(NSA)+NVIR(NSA)
    if (NB2 > NA) NB2 = NA
    if (NB2 < NB1) GO TO 141
    do NB=NB1,NB2
      IIC = IIC+1
      IND = IND+1
      S(IND) = S(IND)+COP(IIC)*C(IND)
      if (IIC < nCOP) GO TO 142
      call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
      IIC = 0
142   continue
    end do
141 continue
  end do
100 continue
end do

return

end subroutine DIAGC
