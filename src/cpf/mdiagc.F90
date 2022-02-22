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

subroutine MDIAGC(JSY,C,S,W,THET,ENP,NII)

use Constants, only: One, Half, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: JSY(*), NII
real(kind=wp) :: C(*), S(*), W(*), THET(NII,NII), ENP(*)
#include "cpfmcpf.fh"
#include "files_cpf.fh"
integer(kind=iwp) :: IADD25, IIC, ILIM, IND, INDA, IRL, NA, NA1, NA2, NB, NB1, NB2, NSA, NSS
real(kind=wp) :: ENPQ, FACS, FACW
integer(kind=iwp), external :: JSUNP_CPF

IADD25 = IAD25S
call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
IIC = 0
IND = 0
ILIM = 4
if (IFIRST /= 0) ILIM = 2
IRL = IRC(ILIM)
do INDA=1,IRL
  NSS = MUL(JSUNP_CPF(JSY,INDA),LSYM)
  ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
  FACS = sqrt(ENP(INDA))*sqrt(ENP(INDA))/ENPQ
  FACW = (FACS*(Two-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
  if (INDA <= IRC(1)) then
    IIC = IIC+1
    IND = IND+1
    S(IND) = S(IND)+FACS*COP(IIC)*C(IND)
    W(IND) = W(IND)+FACW*COP(IIC)*C(IND)
    if (IIC >= nCOP) then
      call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
      IIC = 0
    end if
  else if (INDA <= IRC(2)) then
    NA1 = NSYS(NSS)+1
    NA2 = NSYS(NSS+1)
    do NA=NA1,NA2
      IIC = IIC+1
      IND = IND+1
      S(IND) = S(IND)+FACS*COP(IIC)*C(IND)
      W(IND) = W(IND)+FACW*COP(IIC)*C(IND)
      if (IIC >= nCOP) then
        call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
        IIC = 0
      end if
    end do
  else
    do NA=1,NVIRT
      NSA = MUL(NSS,NSM(LN+NA))
      NB1 = NSYS(NSA)+1
      NB2 = NSYS(NSA+1)
      if (NB2 > NA) NB2 = NA
      do NB=NB1,NB2
        IIC = IIC+1
        IND = IND+1
        S(IND) = S(IND)+FACS*COP(IIC)*C(IND)
        W(IND) = W(IND)+FACW*COP(IIC)*C(IND)
        if (IIC >= nCOP) then
          call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
          IIC = 0
        end if
      end do
    end do
  end if
end do

return

end subroutine MDIAGC
