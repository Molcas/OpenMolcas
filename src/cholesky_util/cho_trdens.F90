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
subroutine CHO_TRDENS(irc,DLT,Salpha,istate,jstate,iType,DoExch,labB)

!C
!C Author:  F. Aquilante and F. Segatta, Bologna ca. Oct 2015
!C Modified by A. Kaiser, 2022
!C     Coulomb term from TDMAT  (level 2 BLAS) from Cholesky vectors
!C
!C --- V(J) = sum_gd  TDMAT(gd) * Lgd,J
!C
!C*********************************************************************

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use Data_Structures, only: DSBA_Type, SBA_Type, Allocate_DT, &
                           Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit real(kind=wp) (A-H,O-Z)
integer(kind=iwp), intent(inout) :: irc
type(SBA_Type), target :: Ya(1)
type(DSBA_Type), intent(in) :: DLT, Salpha(1)

#include "chotime.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "debug.fh"
integer(kind=iwp) :: i, iBatch, iLoc, IVEC2, iVrs, JNUM, JRED, JRED1, JRED2, JSYM, JVEC, LREAD, LWork, MUSED, nBatch, nDen, nRS, &
                     NUMV, nVec, nVrs, iSym, iCase, LuT_, LuT2, LuT1, LuT, dimX
integer(kind=iwp), intent(in) :: iType, istate, jstate
integer(kind=iwp), external :: isFreeUnit
real(kind=wp) :: TCC1, TCC2, tcoul(2), TCR1, TCR2, TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TWC1, &
                 TWC2, TWR1, TWR2
logical, intent(in) :: DoExch, labB
logical :: add
character(len=50) :: CFmt
real(kind=wp), allocatable :: Drs(:), Frs(:), Lrs(:), VJ(:)
character(len=*), parameter :: SECNAM = 'CHO_TRDENS'
character(len=13) :: filnam
character(len=20) :: filnam1, filnam2
logical :: DoRead

#ifdef _DEBUGPRINT_
Debug= .true.
#else
Debug= .false.
#endif
irc=0

IREDC = -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For Coulomb only, the vectors symmetry is restricted to 1
JSYM = 1
if (NumCho(JSYM) < 1) return
write(filnam,'(A,I1,A,I3.3,A,I3.3)') 'WK_C',iType,'_',ISTATE,'_',JSTATE
LuT_ = 10
LuT = isFreeUnit(LuT_)
call molcas_open(LuT,filnam)
write(LuT,*) NumCho(JSYM)

! C create and open file WKX for DoExch
if (DoExch) then
  write(filnam1,'(A,I1,I2.2,I2.2)')'X',iType,ISTATE,JSTATE
  LuT2 = 11
  LuT1 = isFreeUnit(LuT_)
  call DANAME(LuT1,filnam1)
  LuT2 = isFreeUnit(LuT_)
  write(filnam2,'(A)') 'nRedX'
  call DANAME(LuT2,filnam2)

  iSym = 1
  nRS = nDimRS(JSYM,1)
  dimX = nBas(iSym)*nBas(iSym)*NumCho(JSYM)
  call GetMem('idimX_ptr','Allo','Inte',idimX_ptr,1)
  Work(idimX_ptr) = dimX
  iAddr = 0
  call dDaFile(LuT2,1,Work(idimX_ptr),1,iAddr)
  call DACLOS(LuT2)
  call GetMem('idimX_ptr','Free','Inte',idimX_ptr,1)
end if

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

do i=1,2          ! 1 --> CPU   2 --> Wall
  tread(i) = Zero ! time for reading the vectors
  tcoul(i) = Zero ! time for computing Coulomb
end do

iLoc = 3 ! use scratch location in reduced index arrays

JRED1 = InfVec(1,2,jSym)            ! red set of the 1st vec
JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec

do JRED=JRED1,JRED2

  ! Memory management section -----------------------------
  call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

  if (nVrs == 0) cycle

  if (nVrs < 0) then
    write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs < 0. STOP!!'
    call abend()
  end if

  call Cho_X_SetRed(irc,iLoc,JRED) ! set index arrays at iLoc
  if (irc /= 0) then
    write(u6,*) SECNAM//'cho_X_setred non-Zero return code. rc= ',irc
    call abend()
  end if

  nRS = nDimRS(JSYM,JRED)

  call mma_allocate(Drs,nRS,Label='Drs')
  call mma_allocate(Frs,nRS,Label='Frs')
  Drs(:) = Zero
  Frs(:) = Zero

  call mma_maxDBLE(LWork)
  nVec = min(LWORK/(nRS+1),nVrs)

  if (nVec < 1) then
    write(u6,*) SECNAM//': Insufficient memory for batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'min. mem. need= ',nRS+1
    irc = 33
    call Abend()
    nBatch = -9999 ! dummy assignment
  end if

  LREAD = nRS*nVec

  call mma_allocate(Lrs,LREAD,Label='Lrs')

  call mma_allocate(VJ,nVec,Label='VJ')
  add = .false.
  nDen = 1
  ! Transform the density to reduced storage

  call swap_full2rs(irc,iLoc,nRS,nDen,JSYM,[DLT],Drs,add)

  nBatch = (nVrs-1)/nVec+1

    iSym = 1

  do iBatch=1,nBatch

    if (iBatch == nBatch) then
      JNUM = nVrs-nVec*(nBatch-1)
    else
      JNUM = nVec
    end if

    JVEC = nVec*(iBatch-1)+iVrs
    IVEC2 = JVEC-1+JNUM

    call CWTIME(TCR1,TWR1)
    call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,JRED,MUSED)
    if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
      irc = 77
      return
    end if

    call CWTIME(TCR2,TWR2)
    tread(1) = tread(1)+(TCR2-TCR1)
    tread(2) = tread(2)+(TWR2-TWR1)

    ! ************ BEGIN COULOMB CONTRIBUTION  ****************
    !
    !-Computing the intermediate vector V(J)
    !
    ! Contraction with the density matrix
    ! -----------------------------------
    ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
    !==========================================================

    call CWTIME(TCC1,TWC1)

    call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs,1,Zero,VJ,1)

! Coulomb intermediate wrote to file
    write(LuT,*) (VJ(k),k=1,JNUM)

    call CWTIME(TCC2,TWC2)
    tcoul(1) = tcoul(1)+(TCC2-TCC1)
    tcoul(2) = tcoul(2)+(TWC2-TWC1)


! C --- Exchange term
! C ******************  L(ps,{#J}) * G(q,s)  ****************

    if (labB) then
      iSwap = 1 ! L(k,b,J) are returned
    else
      iSwap = 0 ! L(a,k,J) are returned
    end if

    if (doexch) then

      iSym=1
      iCase=0
      DoRead = .false.

      call Allocate_DT(Ya(1),nBas,nBas,nVec,JSYM,nSYm,iCase)

      kMOs = 1
      nMOs = 1

!********************************************************

      ! sending in Salpha and then print ddot of the result
      call Cho_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM, &
                         iSwap,IREDC,nMOs,kMOs,Salpha,&
                         Ya,DoRead)
      if (Debug) then
        write(u6,*) 'ddot Y * Y from Salpha'
        write(u6,*) ddot_(dimX,Ya(1)%A0,1,Ya(1)%A0,1)
      end if
      iAddr = 0
      call dDaFile(LuT1,1,Ya(1)%A0(1),dimX,iAddr)
      call Deallocate_DT(Ya(1))

    end if

  end do  !end batch loop
  close(LuT)
  if (DoExch) then
    call DACLOS(LuT1)
  end if

  ! free memory
  call mma_deallocate(VJ)
  call mma_deallocate(Lrs)
  call mma_deallocate(Frs)
  call mma_deallocate(Drs)

end do ! loop over red sets


call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky SCF timing from '//SECNAM
  write(u6,CFmt) '-----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in SCF-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(u6,'(6X,A)')
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
      call TRIPRT('Coulomb Fmat',' ',FLT(1)%SB(ISYM)%A1,NB)
    end if
  end do

end if
#endif

irc = 0

return

end subroutine CHO_TRDENS
