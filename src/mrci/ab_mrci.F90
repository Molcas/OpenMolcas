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

subroutine AB_MRCI(ICSPCK,INTSYM,INDX,C,S,FC,A,B,FK)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension ICSPCK(*), INTSYM(*), INDX(*), C(*), S(*), FC(*), A(*), B(*), FK(*)
dimension IPOA(9), IPOF(9)
!PAM97 external UNPACK
!PAM97 integer UNPACK
!Statement function
!PAM96 JSYM(L) = UNPACK(INTSYM((L+9)/10),3*mod(L-1,10)+1,3)+1
JSYM(L) = JSUNP(INTSYM,L)

call CSCALE(INDX,INTSYM,C,SQ2)
call CSCALE(INDX,INTSYM,S,SQ2INV)
NCLIM = 4
if (IFIRST /= 0) NCLIM = 2
! MOVE FOCK MATRIX TO FK IN SYMMETRY BLOCKS
call IPO(IPOF,NVIR,MUL,NSYM,1,-1)
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      NAB = IROW(LN+NA)+LN+NB
      if (NB > NA) NAB = IROW(LN+NB)+LN+NA
      FK(IAB) = FC(NAB)
      if (NA == NB) FK(IAB) = 0.0d00
    end do
  end do
end do
ITAIL = IRC(NCLIM)
do INDA=1,ITAIL
  if (INDA <= IRC(1)) GO TO 40
  MYSYM = JSYM(INDA)
  MYL = MUL(MYSYM,LSYM)
  INMY = INDX(INDA)+1
  if (INDA > IRC(2)) GO TO 25
  ! DOUBLET-DOUBLET INTERACTIONS
  if (NVIR(MYL) /= 0) then
    call FZERO(A,NVIR(MYL))
    call FMMM(FK(IPOF(MYL)+1),C(INMY),A,NVIR(MYL),1,NVIR(MYL))
    call DAXPY_(NVIR(MYL),1.0d00,A,1,S(INMY),1)
  end if
  GO TO 40
  ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
25 IFT = 1
  if (INDA > IRC(3)) IFT = 0
  call IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
  !PAM97 IN = 0
  !PAM97 TSUM = 0.0D00
  do IASYM=1,NSYM
    IAB = IPOF(IASYM+1)-IPOF(IASYM)
    if (IAB == 0) GO TO 70
    ICSYM = MUL(MYL,IASYM)
    if (NVIR(ICSYM) == 0) GO TO 70
    if (MYL /= 1) GO TO 30
    if (IFT == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
    if (IFT == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
    NAA = NVIR(IASYM)*NVIR(IASYM)
    call FZERO(B,NAA)
    call FMMM(FK(IPOF(IASYM)+1),A,B,NVIR(IASYM),NVIR(IASYM),NVIR(IASYM))
    call FZERO(A,NAA)
    call DAXPY_(NAA,1.0d00,B,1,A,1)
    if (IFT == 1) GO TO 230
    call SIADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
    call FZERO(A,NAA)
    GO TO 70
230 call TRADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
    call FZERO(A,NAA)
    GO TO 70
30  NAC = NVIR(IASYM)*NVIR(ICSYM)
    call FZERO(A,NAC)
    if (IASYM > ICSYM) GO TO 31
    call FMMM(FK(IPOF(IASYM)+1),C(INMY+IPOA(ICSYM)),A,NVIR(IASYM),NVIR(ICSYM),NVIR(IASYM))
    call DAXPY_(NAC,1.0d00,A,1,S(INMY+IPOA(ICSYM)),1)
    GO TO 70
31  call FMMM(C(INMY+IPOA(IASYM)),FK(IPOF(IASYM)+1),A,NVIR(ICSYM),NVIR(IASYM),NVIR(IASYM))
    call DAXPY_(NAC,1.0d00,A,1,S(INMY+IPOA(IASYM)),1)
    GO TO 70
70  continue
  end do
40 continue
end do
call CSCALE(INDX,INTSYM,C,SQ2INV)
call CSCALE(INDX,INTSYM,S,SQ2)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(ICSPCK)

end subroutine AB_MRCI
