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

use mrci_global, only: IAD25S, IFIRST, IRC, LN, LSYM, Lu_25, NSM, NVIR, NVIRP, NVIRT
use guga_util_global, only: COP, nCOP
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*)
real(kind=wp), intent(in) :: C(*)
real(kind=wp), intent(inout) :: S(*)
integer(kind=iwp) :: IADD25, IIC, ILIM, IND, INDA, IRL, NA, NA1, NA2, NB, NB1, NB2, NSA, NSS
integer(kind=iwp), external :: JSUNP

IADD25 = IAD25S
call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
IIC = 0
IND = 0
ILIM = 4
if (IFIRST /= 0) ILIM = 2
IRL = IRC(ILIM)
do INDA=1,IRL
  NSS = MUL(JSUNP(INTSYM,INDA),LSYM)
  if (INDA <= IRC(1)) then
    IIC = IIC+1
    IND = IND+1
    S(IND) = S(IND)+COP(IIC)*C(IND)
    if (IIC >= nCOP) then
      call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
      IIC = 0
    end if
  else if (INDA <= IRC(2)) then
    NA1 = NVIRP(NSS)+1
    NA2 = NVIRP(NSS)+NVIR(NSS)
    do NA=NA1,NA2
      IIC = IIC+1
      IND = IND+1
      S(IND) = S(IND)+COP(IIC)*C(IND)
      if (IIC >= nCOP) then
        call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
        IIC = 0
      end if
    end do
  else
    do NA=1,NVIRT
      NSA = MUL(NSS,NSM(LN+NA))
      NB1 = NVIRP(NSA)+1
      NB2 = NVIRP(NSA)+NVIR(NSA)
      if (NB2 > NA) NB2 = NA
      do NB=NB1,NB2
        IIC = IIC+1
        IND = IND+1
        S(IND) = S(IND)+COP(IIC)*C(IND)
        if (IIC >= nCOP) then
          call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
          IIC = 0
        end if
      end do
    end do
  end if
end do

return

end subroutine DIAGC
