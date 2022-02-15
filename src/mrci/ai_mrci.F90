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

!subroutine AI(INTSYM,INDX,C,S,FC,BUFIN,IBUFIN,A,B,FK,DBK,KTYP)
subroutine AI_MRCI(INTSYM,INDX,C,S,FC,A,B,FK,DBK,KTYP)

use mrci_global, only: IRC, IREST, IROW, ITER, LASTAD, LN, LSYM, Lu_60, LUSYMB, NBITM3, NBTRI, NSM, NSYM, NVIR, NVIRP, NVIRT, SQ2, &
                       SQ2INV
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: INTSYM(*), INDX(*), KTYP
real(kind=wp) :: C(*), S(*), FC(*), A(*), B(*), FK(*), DBK(*) !, BUFIN(*), IBUFIN(*)
#include "cop.fh"
integer(kind=iwp) :: i, IADR, ICHK, ICP1, ICP2, IFT, II, IJ, IJOLD, ILEN, IND, INDA, INDB, INDI, INMY, INNY, IOUT, IPOB(9), ITYP, &
                     J, LENGTH, MYEXTS, MYINTS, NA, NAK, NI, NJ, NK, NKM, NOTT, NOVST, NSA, NSI, NSIJ, NSJ, NSK, NVIRA, NVM, NVT, &
                     NYEXTS, NYINTS
real(kind=wp) :: COPI, SGN, TERM
integer(kind=iwp), allocatable :: iBuf(:)
real(kind=wp), allocatable :: Buf(:)
integer(kind=iwp), external :: JSUNP
real(kind=r8), external :: DDOT_
!Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

! KTYP=0,  (A/I)   INTEGRALS
! KTYP=1,  (AI/JK) INTEGRALS

call mma_allocate(Buf,NBITM3,label='BUF')
call mma_allocate(iBuf,NBITM3+2,label='IBUF')

call CSCALE(INDX,INTSYM,C,SQ2)
call CSCALE(INDX,INTSYM,S,SQ2INV)
NVT = IROW(NVIRT+1)
ICHK = 0
IJOLD = 0
NK = 0
NSA = 1
NOTT = LN*(LN+1)
NOVST = LN*NVIRT+1+NVT
!PAM97 New portable code:
!PAM04 NBCMX3 = (RTOI*NBSIZ3-2)/(RTOI+1)
!PAM04 IBOFF3 = RTOI*NBCMX3
!PAM04 IBBC3 = IBOFF3+NBCMX3+1
!PAM04 IBDA3 = IBBC3+1

if (KTYP == 0) IADD10 = IAD10(9)
if (KTYP == 1) IADD10 = IAD10(7)
! READ A COP BUFFER
do
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN < 0) exit
  ! LOOP THROUGH THE COP BUFFER:
  do II=1,ILEN
    IND = ICOP1(II)
    if (ICHK /= 0) then
      ! BEGIN A RATHER LONG IF-BLOCK.
      ! ICHK FLAG IS SET. THIS SIGNALS THAT PREVIOUS IND WAS 0, WHICH IS
      ! USED TO INDICATE CHANGE TO A NEW BLOCK OF COUPLING COEFFICIENTS.
      ! RESET ICHK FLAG.
      ICHK = 0
      if (KTYP == 0) then
        ! AI CASE. SAVE INTERNAL ORBITAL INDEX IN NK:
        NK = IND
        IJOLD = NK
        NSK = NSM(NK)
        NSA = NSK
      else
        ! AIJK CASE. UNPACK INTERNAL ORBITAL INDICES INTO NI,NJ,NK:
        INDI = IND
        !NI = mod(INDI,2**10)
        !NJ = mod(INDI/2**10,2**10)
        !NK = mod(INDI/2**20,2**10)
        NI = ibits(INDI,0,10)
        NJ = ibits(INDI,10,10)
        NK = ibits(INDI,20,10)
        NSI = NSM(NI)
        NSJ = NSM(NJ)
        NSK = NSM(NK)
        NSIJ = MUL(NSI,NSJ)
        NSA = MUL(NSIJ,NSK)
        IJ = IROW(NI)+NJ
        if (IJ /= IJOLD) then
          ! NEW INTERNAL PAIR IJ. LOAD A NEW SET OF INTEGRALS INTO FC:
          IJOLD = IJ
          IADR = LASTAD(NOVST+NOTT+IJ)
          call FZERO(FC,NBTRI)

          do
            !PAM04 call dDAFILE(Lu_60,2,IBUFIN,NBSIZ3,IADR)
            call iDAFILE(Lu_60,2,iBuf,NBITM3+2,IADR)
            call dDAFILE(Lu_60,2,Buf,NBITM3,IADR)
            LENGTH = iBuf(NBITM3+1)
            IADR = iBuf(NBITM3+2)
            !PAM04 LENGTH = IBUFIN(IBBC3)
            !PAM04 IADR = IBUFIN(IBDA3)
            !call SCATTER(LENGTH,FC,IBUFIN(IBOFF3+1),BUFIN)
            do i=1,length
              !PAM04 fc(IBUFIN(IBOFF3+i)) = bufin(i)
              fc(iBuf(i)) = Buf(i)
            end do
            if (IADR == -1) exit
          end do
        end if
      end if
      ! FOR THIS PARTICULAR K, TRANSFER FC(NK,NA) TO ARRAY FK:
      NVIRA = NVIR(NSA)
      do I=1,NVIRA
        NA = NVIRP(NSA)+I
        NAK = IROW(LN+NA)+NK
        FK(I) = FC(NAK)
      end do
      ! END OF THE LONG IF-BLOCK.
    else if (IND == 0) then
      ! IND=0 SIGNALS SWITCH TO A NEW SET OF INTEGRALS.
      ICHK = 1
    else if (NVIRA /= 0) then
      ! WE ARE PROCESSING A COUPLING COEFFICIENT AS USUAL.
      !ITYP = mod(IND,2**6)
      !ICP2 = mod(IND/2**6,2**13)
      !ICP1 = mod(IND/2**19,2**13)
      ITYP = ibits(IND,0,6)
      ICP2 = ibits(IND,6,13)
      ICP1 = ibits(IND,19,13)
      if (ITYP <= 1) then
        ! ITYP=1. VALENCE-SINGLES CASE.
        INDA = ICP1
        INDB = IRC(1)+ICP2
        INNY = INDX(INDB)+1
        COPI = COP(II)*C(INDA)
        call DAXPY_(NVIRA,COPI,FK,1,S(INNY),1)
        TERM = DDOT_(NVIRA,FK,1,C(INNY),1)
        S(INDA) = S(INDA)+COP(II)*TERM
      else if ((ITER /= 1) .or. (IREST /= 0)) then
        INDA = IRC(1)+ICP1
        INDB = IRC(ITYP)+ICP2
        INMY = INDX(INDA)+1
        INNY = INDX(INDB)+1
        MYINTS = JSYM(INDA)
        NYINTS = MUL(MYINTS,NSA)
        MYEXTS = MUL(MYINTS,LSYM)
        NYEXTS = MUL(NYINTS,LSYM)
        IFT = 0
        if (ITYP == 2) IFT = 1
        call IPO(IPOB,NVIR,MUL,NSYM,NYEXTS,IFT)
        NVM = NVIR(MYEXTS)
        call FZERO(DBK,NVIRA)
        call DAXPY_(NVIRA,COP(II),FK,1,DBK,1)
        if (NYEXTS == 1) then
          if (IFT == 0) call SQUAR(C(INNY+IPOB(MYEXTS)),A,NVM)
          if (IFT == 1) call SQUARM(C(INNY+IPOB(MYEXTS)),A,NVM)
          call FZERO(B,NVM)
          call FMMM(DBK,A,B,1,NVM,NVIRA)
          call DAXPY_(NVM,One,B,1,S(INMY),1)
          SGN = One
          if (IFT == 1) SGN = -One
          IOUT = INNY+IPOB(MYEXTS)-1
          do I=1,NVM
            do J=1,I
              IOUT = IOUT+1
              TERM = DBK(I)*C(INMY+J-1)+SGN*DBK(J)*C(INMY+I-1)
              S(IOUT) = S(IOUT)+TERM
            end do
            if (IFT == 1) cycle
            TERM = DBK(I)*C(INMY+I-1)
            S(IOUT) = S(IOUT)-TERM
          end do
        else
          NKM = NVIRA*NVM
          call FZERO(B,NVM)
          if (NSA <= MYEXTS) then
            if (IFT == 1) call VNEG(DBK,1,DBK,1,NVIRA)
            call FMMM(DBK,C(INNY+IPOB(MYEXTS)),B,1,NVM,NVIRA)
            call DAXPY_(NVM,One,B,1,S(INMY),1)
            call FZERO(B,NKM)
            call FMMM(DBK,C(INMY),B,NVIRA,NVM,1)
            call DAXPY_(NKM,One,B,1,S(INNY+IPOB(MYEXTS)),1)
          else
            call FMMM(C(INNY+IPOB(NSA)),DBK,B,NVM,1,NVIRA)
            call DAXPY_(NVM,One,B,1,S(INMY),1)
            call FZERO(B,NKM)
            call FMMM(C(INMY),DBK,B,NVM,NVIRA,1)
            call DAXPY_(NKM,One,B,1,S(INNY+IPOB(NSA)),1)
          end if
        end if
      end if
    end if
  end do
end do
call CSCALE(INDX,INTSYM,C,SQ2INV)
call CSCALE(INDX,INTSYM,S,SQ2)
call mma_deallocate(Buf)
call mma_deallocate(iBuf)

return

end subroutine AI_MRCI
