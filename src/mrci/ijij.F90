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

subroutine IJIJ(INTSYM,HDIAG,FIJIJ)

use mrci_global, only: IAD25S, IDVER, IRC, IROW, IVVER, LN, LSYM, Lu_25, Lu_27, LUSYMB, NSM, NVIR, NVIRT, NVPAIR, NVIRP, POTNUC
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*)
real(kind=wp), intent(_OUT_) :: HDIAG(*)
real(kind=wp), intent(in) :: FIJIJ(*)
integer(kind=iwp) :: IAD27, IADD10, IADD25, ICHK, ICOUP, ICOUPS, IFS, II, IIJ, IIJ1, IIJ2, IJJ, INB, IND, INDI, INS, IOUT, IREF0, &
                     ITYP, IVL, IVSAVE, J, JJ, KK, LENGTH, NA, NA1, NA2, NB, NB1, NB2, NSA, NSS
real(kind=wp) :: HCOUT(nCOP), TERM
integer(kind=iwp), external :: JSUNP

!------
! POW: Unnecessary but warning stopping initializations
INB = -1234567
!------
IADD25 = IAD25S
IAD27 = 0
IREF0 = 1
call dDAFILE(Lu_27,2,HDIAG,IRC(1),IAD27)

!write(u6,*) ' Hdiag'
!write(u6,*) (Hdiag(i),i=1,IRC(1))

IFS = 0
IVL = 0
IVSAVE = 0
ICOUPS = 0
ICOUP = 0
NSS = 1
IOUT = 0
ICHK = 0
IADD10 = IAD10(3)
TERM = Zero
do
  ! READ A COP BUFFER:
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  LENGTH = ICOP1(nCOP+1)
  if (LENGTH < 0) exit
  ! LOOP OVER THE COP BUFFER:
  do II=1,LENGTH
    IND = ICOP1(II)
    if (ICHK == 0) then
      if (IND == 0) then
        ICHK = 1
      else
        !ITYP = mod(IND,2)
        !IJJ = mod(IND/2,2**11)
        ITYP = ibits(IND,0,1)
        IJJ = ibits(IND,1,11)
        if (ITYP == 0) TERM = COP(II)*FIJIJ(IJJ)
        if (IVL == IVVER) then
          INB = ICOUP
          HDIAG(INB) = HDIAG(INB)+TERM
        else if (IVL == IDVER) then
          INB = 0
          NA1 = NVIRP(NSS)+1
          NA2 = NVIRP(NSS)+NVIR(NSS)
          do NA=NA1,NA2
            INB = INB+1
            if (ITYP == 1) then
              IIJ = IROW(LN+NA)+IJJ
              TERM = COP(II)*FIJIJ(IIJ)
            end if
            HDIAG(INB) = HDIAG(INB)+TERM
          end do
        else
          INB = 0
          do NA=1,NVIRT
            NSA = MUL(NSS,NSM(LN+NA))
            NB1 = NVIRP(NSA)+1
            NB2 = NVIRP(NSA)+NVIR(NSA)
            if (NB2 > NA) NB2 = NA
            if (NB2 < NB1) cycle
            IIJ1 = IROW(LN+NA)+IJJ
            do NB=NB1,NB2
              INB = INB+1
              if (ITYP == 1) then
                IIJ2 = IROW(LN+NB)+IJJ
                TERM = COP(II)*(FIJIJ(IIJ1)+FIJIJ(IIJ2))
              end if
              HDIAG(INB) = HDIAG(INB)+TERM
            end do
          end do
        end if
      end if
    else
      ICHK = 0
      INDI = IND
      !ICOUP = mod(INDI,2**16)
      !IVL = mod(INDI/2**16,2**8)
      ICOUP = ibits(INDI,0,16)
      IVL = ibits(INDI,16,8)
      ICHK = 0
      INS = 1
      if (IVSAVE == IVVER) then
        INS = ICOUPS
        INB = ICOUPS
      end if
      if (INB /= 0) then
        do J=INS,INB
          IOUT = IOUT+1
          HCOUT(IOUT) = HDIAG(J)
          if (IOUT < nCOP) cycle
          if (IFS == 0) then
            POTNUC = HCOUT(IREF0)
            IFS = 1
          end if
          do KK=1,nCOP
            HCOUT(KK) = HCOUT(KK)-POTNUC
          end do
          call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
          IOUT = 0
        end do
      end if
      if (IVL /= IVVER) then
        JJ = IRC(IVL)+ICOUP
        NSS = MUL(JSUNP(INTSYM,JJ),LSYM)
        if (IVL == IDVER) then
          INB = NVIR(NSS)
        else
          INB = NVPAIR(NSS)
        end if
        if (INB > 0) call dDAFILE(Lu_27,2,HDIAG,INB,IAD27)
      end if
      IVSAVE = IVL
      ICOUPS = ICOUP
    end if
  end do
end do
! EMPTY LAST BUFFER
do J=1,INB
  IOUT = IOUT+1
  HCOUT(IOUT) = HDIAG(J)
  if (IOUT < nCOP) cycle
  if (IFS == 0) then
    POTNUC = HCOUT(IREF0)
    IFS = 1
  end if
  do KK=1,nCOP
    HCOUT(KK) = HCOUT(KK)-POTNUC
  end do
  call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
  IOUT = 0
end do
if (IFS == 0) then
  POTNUC = HCOUT(IREF0)
  IFS = 1
end if
do KK=1,IOUT
  HCOUT(KK) = HCOUT(KK)-POTNUC
end do
call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)

return

end subroutine IJIJ
