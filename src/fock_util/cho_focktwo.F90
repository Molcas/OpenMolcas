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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
!***********************************************************************
!  Author : F. Aquilante
!
!  Purpose:
!          For a given set of total symmetric AO-densities
!          (hermitian matrix representation) this routine
!          computes the Coulomb and Exchange contributions
!          to the Fock matrix from Cholesky vectors in full storage
!
!          Coulomb term (Drs is the lower triangular matrix
!                        of r=s symmetry block. The diagonals
!                        are scaled by a factor 1/2 already in
!                        the packed density DLT, given as input)
!
!    F(pq) = F(pq) + FactC * Sum_J Lpq,J * (Sum_rs Drs*Lrs,J)
!
!          Exchange term (Drs is the squared matrix
!                        of r=s symmetry block.)
!
!    F(pq) = F(pq) - FactX * Sum_rJ  Lpr,J * (Sum_s Drs*Lqs,J)
!
!    Integral-fashion formula:
!
!    F(pq) =  Sum_rs Drs * [(pq|rs)  - 0.5 (ps|rq)]
!
!----------------------------------------------------------------------
!  DoCoulomb(nDen) :  specifies the densities for which coulomb term
!                     has to be computed
!
!  DoExchange(nDen) : specifies the densities for which exchange term
!                     has to be computed
!
!  FactC(nDen) : factor for the coulomb of the corresponding density
!  FactX(nDen) : factor for the exchange of the corresponding density
!
!  ip{X}LT(nDen) : pointer to the array containing {X} in LT storage
!  ip{X}FS(nDen) : pointer to the array containing {X} in SQ storage
!    {X=D,F --- Density, Fock matrix}
!
!  pNocc(nDen) : pointer to the array of the Occupation numbers
!                 for the corresponding density
!
!  MinMem(nSym) : minimum amount of memory required to read
!                 a single Cholesky vector in full storage
!
!***********************************************************************

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Fock_util_global, only: Lunit
use Data_Structures, only: Deallocate_DT, DSBA_Type, Integer_Pointer, SBA_Type
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nDen, MinMem(nSym)
logical(kind=iwp), intent(in) :: DoCoulomb(nDen), DoExchange(nDen)
real(kind=wp), intent(in) :: FactC(nDen), FactX(nDen)
type(DSBA_Type), intent(in) :: DLT(nDen), DSQ(nDen)
type(DSBA_Type), intent(inout) :: FLT(nDen), FSQ(nDen)
type(Integer_Pointer), intent(in) :: pNocc(nDen)
#include "chotime.fh"
integer(kind=iwp) :: iBatch, iE, iS, iSkip(nSym), iSym, ISYMA, ISYMB, ISYMD, ISYMG, iSymp, ISYMQ, iSymr, iSymr_Occ, iSyms, iVec, &
                     jD, jDen, jE, jR, jS, jSR, jSym, JVEC, k, kRdMem, kSym, l, lu, LWORK, MaxSym, Naa, nb, nBatch, nd, ndim3, nk, &
                     Nmax, np, NumB, NumCho(nSym), NumV, nVec
real(kind=wp) :: TC1X1, TC1X2, TC2X1, TC2X2, TCC1, TCC2, tcoul(2), TCR1, TCR2, TCREO1, TCREO2, TCREO3, TCREO4, texch(2), TOTCPU, &
                 TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TW1X1, TW1X2, TW2X1, TW2X2, TWC1, TWC2, TWR1, TWR2, &
                 TWREO1, TWREO2, TWREO3, TWREO4
logical(kind=iwp) :: DoSomeC, DoSomeX
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
character(len=6) :: Fname
type(SBA_Type), target :: LqJs, Wab
real(kind=wp), pointer :: LrJs(:,:,:) => null(), VJ(:) => null(), XdJb(:) => null(), XpJs(:) => null()
character(len=*), parameter :: BaseNm = 'CHFV', SECNAM = 'CHO_FOCKTWO'
integer(kind=iwp), external :: isfreeunit

!*************************************************
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in SCF-debug
#endif

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero ! time read/rreorder vectors
tcoul(:) = zero ! time for computing Coulomb
texch(:) = zero ! time for computing Exchange

! Tests on the type of calculation
DoSomeC = .false.
DoSomeX = .false.
MaxSym = 0

jD = 0
do while ((jD < nDen) .and. (.not. DoSomeX))
  jD = jD+1
  DoSomeX = DoExchange(jD)
end do
jD = 0
do while ((jD < nDen) .and. (.not. DoSomeC))
  jD = jD+1
  DoSomeC = DoCoulomb(jD)
end do

if (DoSomeX) then
  MaxSym = nSym ! we want to do some exchange
else
  if (DoSomeC) then
    MaxSym = 1 ! we want to do only Coulomb
  else
    return ! we are lazy and won't do anything
  end if
end if

call get_iarray('NumCho',NumCho,nSym)

! *************** BIG LOOP OVER VECTORS SYMMETRY *****************
do jSym=1,MaxSym

  if (NumCho(jSym) < 1) cycle

  ! Open Files
  LUnit(:) = -1
  do ksym=1,nSym
    if (nBas(ksym) /= 0) then
      iSymp = Mul(ksym,jSym)
      if (iSymp >= ksym) then
        lu = 7
        lu = isfreeunit(lu)
        Lunit(iSymp) = lu
        write(Fname,'(A4,I1,I1)') BaseNm,iSymp,ksym
        call DANAME_MF_WA(Lunit(iSymp),Fname)
      end if
    end if
  end do

  ! SET UP THE READING
  ! ------------------
  call mma_maxDBLE(LWORK)

  if (MinMem(jSym) > 0) then
    nVec = min(LWORK/MinMem(jSym),NumCho(jSym))
  else
    ! ***QUIT*** bad initialization
    write(u6,*) 'Cho_FockTwo: bad initialization'
    rc = 99
    call Abend()
    nVec = -9999 ! dummy assignment - avoid compiler warnings
  end if

  !...this is ONLY for debug:
  !      nVec = min(nVec,1)

  if (nVec > 0) then
    nBatch = (NumCho(jSym)-1)/nVec+1
  else
    ! ***QUIT*** insufficient memory
    write(u6,*) 'Cho_FockTwo: Insufficient memory for batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'min. mem. need= ',MinMem(jSym)
    write(u6,*) 'NumCho= ',NumCho(jsym)
    write(u6,*) 'jsym= ',jsym
    rc = 205
    call Abend()
    nBatch = -9999 ! dummy assignment
  end if

  ! *************** BATCHING  *****************
  do iBatch=1,nBatch
    write(u6,*) 'iBatch=',iBatch

    if (iBatch == nBatch) then
      NumV = NumCho(jSym)-nVec*(nBatch-1)
    else
      NumV = nVec
    end if

    iVec = nVec*(iBatch-1)+1

    ! Allocate memory for reading the vectors
    ! We allocate the array as an element of a SBA_type.
    kRdMem = MinMem(jSym)*NumV
    call mma_allocate(Wab%A0,kRdMem,Label='Wab%A0')
    Wab%nsym = nSym
    Wab%isym = jSym
    Wab%icase = 7

    ! setup the skipping flags according to # of Occupied
    do k=1,nSym
      iSkip(k) = 0
      l = Mul(k,jsym) ! L(kl) returned if nOcc(k or l) /= 0
      if (k == l) then
        iSkip(k) = 666 ! always contribute to Coulomb
      else
        do jDen=1,nDen
          iSkip(k) = iSkip(k)+pNocc(jDen)%I1(k)+pNocc(jDen)%I1(l)
        end do
      end if
    end do

    ! Max dimension of a symmetry block
    Nmax = 0
    do iSym=1,nSym
      if ((NBAS(iSym) > Nmax) .and. (iSkip(iSym) /= 0)) then
        Nmax = NBAS(iSym)
      end if
    end do

    call CWTIME(TCR1,TWR1)

    iE = 0
    do ksym=1,nSym

      iSymp = Mul(ksym,jSym)

      iS = iE+1

      nk = nBas(kSym)
      np = nBas(iSymp)
      if (nk*np <= 0) cycle

      if ((iSymp == ksym) .and. (iSkip(iSymp) /= 0)) then
        NumB = nk*(nk+1)/2
        ! Special trick for the vector L11 ; used to store X(a,Jb)
        if ((ksym == 1) .and. (jSym == 1) .and. DoSomeX) then
          iE = iE+(Nmax**2)*NumV
          Wab%SB(iSymp)%A1(1:Nmax**2*NumV) => Wab%A0(iS:iE)
        else
          iE = iE+NumB*NumV
        end if
        Wab%SB(iSymp)%A2(1:NumB,1:NumV) => Wab%A0(iS:iS-1+NumB*NumV)
      else
        if ((iSymp > ksym) .and. (iSkip(iSymp) /= 0)) then
          NumB = nk*np
          iE = iE+NumB*NumV
          Wab%SB(iSymp)%A3(1:np,1:nk,1:NumV) => Wab%A0(iS:iE)
          Wab%SB(iSymp)%A2(1:NumB,1:NumV) => Wab%A0(iS:iE)
        end if
      end if

    end do ! ends the loop over symmetries

    do kSym=1,nSym
      iSymp = Mul(ksym,jSym)
      if (.not. associated(Wab%SB(iSymp)%A2)) cycle
      NumB = size(Wab%SB(iSymp)%A2,1)
      call RdChoVec(Wab%SB(iSymp)%A2,NumB,NumV,iVec,Lunit(iSymp))
    end do

    call CWTIME(TCR2,TWR2)
    tread(1) = tread(1)+(TCR2-TCR1)
    tread(2) = tread(2)+(TWR2-TWR1)

    ! first "free" position in Wab%A0
    iS = iE+1

#   ifdef _DEBUGPRINT_
    write(u6,*) 'Batch ',iBatch,' of  ',nBatch,': NumV = ',NumV
    write(u6,*) 'Total allocated:     ',kRdMem
    write(u6,*) 'iE:                  ',iE
    write(u6,*) 'iS:                ',iS
    write(u6,*) 'JSYM:                ',JSYM
#   endif

    ! Reading Done!

    ! ********************************************************************
    do jDen=1,nDen

      if ((jSym == 1) .and. DoCoulomb(jDen)) then

        ! *************** BEGIN COULOMB CONTRIBUTION  ****************
        !
        ! Computing the intermediate vector V(J)
        !
        ! Contraction with the density matrix
        ! -----------------------------------
        ! V{#J} <- V{#J} + sum_rr L(rr,{#J})*D(rr)
        !==========================================================

        call CWTIME(TCC1,TWC1)

        VJ(1:NumV) => Wab%A0(iS:iE+NumV)
        VJ(:) = Zero

        do iSymr=1,nSym
          if ((nBas(iSymr) /= 0) .and. (pNocc(jDen)%I1(iSymr) /= 0)) then

            Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

            call DGEMV_('T',Naa,NumV,ONE,Wab%SB(iSymr)%A2,Naa,DLT(jDen)%SB(ISYMR)%A1,1,ONE,VJ,1)

          end if
        end do

        ! Contraction with the 2nd vector
        ! -----------------------------------
        ! Fss{#J} <- Fss{#J} + FactC * sum_J L(ss,{#J})*V{#J}
        !==========================================================

        do iSyms=1,nSym
          if (nBas(iSyms) /= 0) then

            Naa = nBas(iSyms)*(nBas(iSyms)+1)/2

            call DGEMV_('N',Naa,NumV,FactC(jDen),Wab%SB(iSyms)%A2,Naa,VJ,1,ONE,FLT(jDen)%SB(ISYMS)%A1,1)

          end if
        end do

        VJ => null()

        call CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1)+(TCC2-TCC1)
        tcoul(2) = tcoul(2)+(TWC2-TWC1)

      end if ! jSym=1 & DoCoulomb

    end do ! loop over densities for the Coulomb term

    ! *************** END COULOMB CONTRIBUTION  ****************

    if ((jSym == 1) .and. DoSomeX) then
      ! WE INCLUDE THE DIAGONAL PART OF THE EXCHANGE HERE!!!!
      !
      ! ********** REORDER AND SQUARE ONE VECTOR AT THE TIME *****
      !     Reorder: L(rs,J) -> L(r,J,s).
      !     CHOVEC(nrs,numv) ---> CHOVEC(nr,numv,ns)

      do iSymr=1,nSym

        iSymr_Occ = 0
        do jDen=1,nDen
          iSymr_Occ = iSymr_Occ+pNocc(jDen)%I1(iSymr)
        end do
        if ((nBas(iSymr) /= 0) .and. (iSymr_Occ > 0)) then

          call CWTIME(TCREO1,TWREO1)

          LrJs(1:nBas(iSymr),1:NUMV,1:nBas(iSymr)) => Wab%A0(iS:iE+nBas(iSymr)**2*NUMV)

          do JVEC=1,NUMV

            do jR=1,NBAS(iSymr)
              do jS=jR,NBAS(iSymr)

                jSR = iTri(jS,jR)

                LrJs(jR,JVEC,jS) = Wab%SB(iSymr)%A2(jSR,JVEC)

                LrJs(jS,JVEC,jR) = Wab%SB(iSymr)%A2(jSR,JVEC)

              end do
            end do
          end do

          call CWTIME(TCREO2,TWREO2)
          tread(1) = tread(1)+(TCREO2-TCREO1)
          tread(2) = tread(2)+(TWREO2-TWREO1)

          do jDen=1,nDen

            if (DoExchange(jDen)) then

              if (pNocc(jDen)%I1(iSymr) /= 0) then

                call CWTIME(TC1X1,TW1X1)

                ndim3 = NBAS(ISYMR)*NUMV*NBAS(ISYMR)
                XpJs(1:nDim3) => Wab%A0(1:nDim3)

                ! Calculate intermediate:
                ! X(p,Js) = Sum(q) D(p,q) * L(q,Js).
                ! -----------------------------------

                call DGEMM_('N','N',NBAS(ISYMR),NUMV*NBAS(ISYMR),NBAS(ISYMR),ONE,DSQ(jDen)%SB(ISYMR)%A2,NBAS(ISYMR),LrJs, &
                            NBAS(ISYMR),ZERO,XpJs,NBAS(ISYMR))

                ! Calculate exchange contribution:
                ! F(p,q) = F(p,q) - FactX Sum(sJ) L(sJ,p) * X(sJ,q).
                ! --------------------------------------------------

                call DGEMM_('T','N',NBAS(ISYMR),NBAS(ISYMR),NBAS(ISYMR)*NUMV,-FactX(jDen),LrJs,NBAS(ISYMR)*NUMV,XpJs, &
                            NBAS(ISYMR)*NUMV,One,FSQ(jDen)%SB(ISYMR)%A2,NBAS(ISYMR))

                call CWTIME(TC1X2,TW1X2)
                texch(1) = texch(1)+(TC1X2-TC1X1)
                texch(2) = texch(2)+(TW1X2-TW1X1)

                XpJs => null()

              end if ! nOcc /= 0

            end if ! DoExchange(jDen)

          end do ! loop over the densities

        end if ! nbas /= 0 & nOcc /= 0

        LrJs => null()

      end do ! loop over Fock mat symmetries

    end if ! jSym=1 and DoSomeX

    ! ****************** END DIAGONAL EXCHANGE **********************

    ! ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************

    if (jSym /= 1) then

      ! Reordering of the Cholesky vectors (ALL AT ONCE)
      !
      !       Reorder: L(qs,J) -> L(qJ,s).
      !
      ! CHOVEC(nq,ns,numv) ---> CHOVEC(nq,numv,ns)
      ! ------------------------------------------
      call CWTIME(TCREO3,TWREO3)

      jE = iE
      do ISYMS=1,NSYM

        ISYMQ = Mul(ISYMS,JSYM)

        if (nBas(iSyms)*nBas(iSymq) <= 0) cycle

        if ((ISYMQ > ISYMS) .and. (iSkip(iSymq) /= 0)) then

          jS = jE+1
          jE = jE+nBas(iSymq)*NumV*nBas(iSyms)
          LqJs%SB(ISYMQ)%A3(1:nBas(iSymq),1:NumV,1:nBas(iSyms)) => Wab%A0(jS:jE)

          do JVEC=1,NUMV

            do jS=1,NBAS(iSyms)

              LqJs%SB(ISYMQ)%A3(:,JVEC,jS) = Wab%SB(ISYMQ)%A3(:,jS,JVEC)

            end do

          end do ! loop over the vectors

        end if ! sym(Q) > sym(S)

      end do ! loop over symmetry blocks

      call CWTIME(TCREO4,TWREO4)
      tread(1) = tread(1)+(TCREO4-TCREO3)
      tread(2) = tread(2)+(TWREO4-TWREO3)

      ! -------------- REORDERING DONE !! --------

      do jDen=1,nDen

        if (DoExchange(jDen)) then
          ! COMPUTE EXCHANGE FOR OFF-DIAGONAL VECTORS

          call CWTIME(TC2X1,TW2X1)

          do ISYMG=1,NSYM

            ISYMB = Mul(ISYMG,JSYM)

            if (nBas(iSymb)*nBas(iSymg) /= 0) then

              ISYMD = ISYMG ! only total symmetric Density
              ISYMA = ISYMB ! block diagonal Fock Matrix

              ! mem used to store the intermediate
              nd = NBAS(ISYMD)
              nb = NBAS(ISYMB)
              XdJb(1:nd*NUMV*nb) => Wab%A0(1:nd*NUMV*nb)

              if ((ISYMG > ISYMB) .and. (iSkip(iSymg) /= 0)) then

                ! ---------------------------
                ! F(a,b) = - D(g,d) * (ad|gb)
                ! ---------------------------
                if (pNocc(jDen)%I1(iSymg) /= 0) then

                  ! Calculate intermediate:
                  ! X(d,Jb) = Sum(g) D(d,g) * L(g,Jb).
                  ! ----------------------------------

                  call DGEMM_('N','N',NBAS(ISYMD),NUMV*NBAS(ISYMB),NBAS(ISYMG),ONE,DSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMD), &
                              LqJs%SB(ISYMG)%A3,NBAS(ISYMG),ZERO,XdJb,NBAS(ISYMD))

                  ! F(a,b) = F(a,b) - Sum(dJ) L(dJ,a) * X(dJ,b).
                  ! -------------------------------------------

                  call DGEMM_('T','N',NBAS(ISYMA),NBAS(ISYMB),NBAS(ISYMD)*NUMV,-FactX(jDen),LqJs%SB(ISYMG)%A3,NBAS(ISYMD)*NUMV, &
                              XdJb,NBAS(ISYMD)*NUMV,ONE,FSQ(jDen)%SB(ISYMB)%A2,NBAS(ISYMA))

                end if
                ! ---------------------------
                ! F(g,d) = - D(a,b) * (ad|gb)
                ! ---------------------------
                if (pNocc(jDen)%I1(iSyma) /= 0) then

                  ! Calculate intermediate:
                  ! X(gJ,b) = Sum(a) L(gJ,a)* D(a,b).
                  ! ----------------------------------

                  call DGEMM_('N','N',NBAS(ISYMG)*NUMV,NBAS(ISYMB),NBAS(ISYMA),ONE,LqJs%SB(ISYMG)%A3,NBAS(ISYMG)*NUMV, &
                              DSQ(jDen)%SB(ISYMB)%A2,NBAS(ISYMA),ZERO,XdJb,NBAS(ISYMG)*NUMV)

                  ! F(g,d) = F(g,d) - Sum(Jb) X(g,Jb) * L(d,Jb).
                  ! -------------------------------------------

                  call DGEMM_('N','T',NBAS(ISYMG),NBAS(ISYMD),NUMV*NBAS(ISYMB),-FactX(jDen),XdJb,NBAS(ISYMG),LqJs%SB(ISYMG)%A3, &
                              NBAS(ISYMD),ONE,FSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMG))

                end if

              end if

              XdJb => null()

            end if
            LqJs%SB(ISYMG)%A3 => null()

          end do ! loop over iorbital symmetries

          call CWTIME(TC2X2,TW2X2)
          texch(1) = texch(1)+(TC2X2-TC2X1)
          texch(2) = texch(2)+(TW2X2-TW2X1)

          ! ************ END "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ***********

        end if ! DoExchange section

      end do ! end of the loop over densities

    end if ! jSym /= 1

    ! Free the memory
    call Deallocate_DT(Wab)

  end do ! end of the batch procedure

  ! Close Files
  do ksym=1,nSym
    if (Lunit(ksym) /= -1) then
      call DACLOS(Lunit(ksym))
      Lunit(ksym) = -1
    end if
  end do

end do ! loop over JSYM

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky SCF timing from '//SECNAM
  write(u6,CFmt) '------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ/REORDER VECTORS                      ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'EXCHANGE                                  ',texch(1),texch(2)

  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in SCF-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM CHO_FOCKTWO.'
  write(u6,'(6X,A)') '***** EXCHANGE MATRIX ***** '
  do jDen=1,nDen
    write(u6,'(6X,A,I2)') 'DENSITY TYPE: ',jDen
    write(u6,'(6X,A,I2)') 'DoExchange: ',DoExchange(jDen)
    write(u6,*)
    if (DoExchange(jDen)) then
      do ISYM=1,NSYM
        NB = NBAS(ISYM)
        if (NB > 0) then
          write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
          call cho_output(FSQ(jDen)%SB(ISYM)%A2,1,NB,1,NB,NB,NB,1,u6)
        end if
      end do
    end if
  end do

end if
#endif

rc = 0

return

end subroutine CHO_FOCKTWO
