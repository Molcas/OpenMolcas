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

subroutine ABTD(ICSPCK,INTSYM,INDX,C1,C2,TDMO,A1,A2,F)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension ICSPCK(*), INTSYM(*), INDX(*), C1(*), C2(*), TDMO(NBAST,NBAST), A1(*), A2(*), F(*)
dimension IPOA(9), IPOF(9)
dimension IOC(55)
!PAM97 external UNPACK
!PAM97 integer UNPACK
!Statement functions
!PAM97 JO(L) = UNPACK(CSPCK((L+29)/30),2*L-(2*L-1)/60*60,2)
JO(L) = ICUNP(ICSPCK,L)
!PAM96 JSYM(L) = UNPACK(INTSYM((L+9)/10),3*mod(L-1,10)+1,3)+1
JSYM(L) = JSUNP(INTSYM,L)

! CALCULATE A) TRANSITION DENSITY ELEMENTS OF TYPE TDMO(A,B)
!           B) DIAGONAL ELEMENTS TDMO(I,I) AND TDMO(A,A)
! SCRATCH SPACES: A1(),A2(), SIZE NEEDED IS NVMAX**2
!                  ALSO F(), SIZE NEEDED IS NVSQ
call CSCALE(INDX,INTSYM,C1,SQ2)
call CSCALE(INDX,INTSYM,C2,SQ2)
NCLIM = 4
if (IFIRST /= 0) NCLIM = 2
! MOVE TRANSITION DENSITY MATRIX TO F IN SYMMETRY BLOCKS
call IPO(IPOF,NVIR,MUL,NSYM,1,-1)
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      F(IAB) = TDMO(LN+NA,LN+NB)
    end do
  end do
end do
II1 = 0
ITAIL = IRC(NCLIM)
do INDA=1,ITAIL
  do I=1,LN
    II1 = II1+1
    IOC(I) = (1+JO(II1))/2
  end do
  if (INDA <= IRC(1)) then
    TSUM = C1(INDA)*C2(INDA)
    GO TO 106
  end if
  MYSYM = JSYM(INDA)
  MYL = MUL(MYSYM,LSYM)
  INMY = INDX(INDA)+1
  if (INDA > IRC(2)) GO TO 25
  ! DOUBLET-DOUBLET INTERACTIONS
  if (NVIR(MYL) == 0) GO TO 40
  IPF = IPOF(MYL)+1
  NVIRA = NVIR(MYL)
  call DGER(NVIRA,NVIRA,1.0d00,C1(INMY),1,C2(INMY),1,F(IPF),NVIRA)
  LNA = LN+NVIRP(MYL)
  TSUM = 0.0d00
  do I=1,NVIRA
    TERM = C1(INMY-1+I)*C2(INMY-1+I)
    IA = LNA+I
    TDMO(IA,IA) = TDMO(IA,IA)+TERM
    TSUM = TSUM+TERM
  end do
  GO TO 106
  ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
25 IFT = 1
  if (INDA > IRC(3)) IFT = 0
  call IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
  TSUM = 0.0d00
  do IASYM=1,NSYM
    IAB = IPOF(IASYM+1)-IPOF(IASYM)
    if (IAB == 0) GO TO 70
    ICSYM = MUL(MYL,IASYM)
    NVIRA = NVIR(IASYM)
    NVIRC = NVIR(ICSYM)
    if (NVIRC == 0) GO TO 70
    if (MYL /= 1) then
      if (IASYM > ICSYM) then
        call MTRANS(C1(INMY+IPOA(IASYM)),1,A1,1,NVIRA,NVIRC)
        call MTRANS(C2(INMY+IPOA(IASYM)),1,A2,1,NVIRA,NVIRC)
      else
        NAC = NVIRA*NVIRC
        if (IFT == 0) then
          call DCOPY_(NAC,C1(INMY+IPOA(ICSYM)),1,A1,1)
          call DCOPY_(NAC,C2(INMY+IPOA(ICSYM)),1,A2,1)
        else
          call VNEG(C1(INMY+IPOA(ICSYM)),1,A1,1,NAC)
          call VNEG(C2(INMY+IPOA(ICSYM)),1,A2,1,NAC)
        end if
      end if
    else
      if (IFT == 0) then
        call SQUAR(C1(INMY+IPOA(IASYM)),A1,NVIRA)
        call SQUAR(C2(INMY+IPOA(IASYM)),A2,NVIRA)
      else
        call SQUARM(C1(INMY+IPOA(IASYM)),A1,NVIRA)
        call SQUARM(C2(INMY+IPOA(IASYM)),A2,NVIRA)
      end if
    end if
    IPF = IPOF(IASYM)+1
    call DGEMM_('N','T',NVIRA,NVIRA,NVIRC,1.0d00,A1,NVIRA,A2,NVIRA,1.0d00,F(IPF),NVIRA)
    INN = 1
    LNC = LN+NVIRP(ICSYM)
    do I=1,NVIRC
      TERM = DDOT_(NVIRA,A1(INN),1,A2(INN),1)
      TSUM = TSUM+TERM
      IC = LNC+I
      TDMO(IC,IC) = TDMO(IC,IC)+TERM
      INN = INN+NVIRA
    end do
70  continue
  end do
  TSUM = TSUM/2
106 continue
  do I=1,LN
    TDMO(I,I) = TDMO(I,I)+IOC(I)*TSUM
  end do
40 continue
end do
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      if (NA == NB) goto 420
      TDMO(LN+NA,LN+NB) = F(IAB)
420   continue
    end do
  end do
end do
call CSCALE(INDX,INTSYM,C1,SQ2INV)
call CSCALE(INDX,INTSYM,C2,SQ2INV)

return

end subroutine ABTD
