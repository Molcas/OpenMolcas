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

subroutine CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
!***********************************************************************
!  Author : F. Aquilante
!
!  Purpose:
!          For a given set of total symmetric AO-densities
!          (hermitian matrix representation) and MO coefficients
!          this routine
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
!          Exchange term (Crk is the squared matrix of MO coeff.
!                             of r=s symmetry block)
!
!    F(pq) = F(pq) - FactX * Sum_Jk [(Sum_r Crk * Lpr,J) * (Sum_s Csk * Lqs,J)]
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
!  ip{X}SQ(nDen) : pointer to the array containing {X} in SQ storage
!    {X=D,F,M --- Density, Fock matrix, MOs coeff.}
!
!  MinMem(nSym) : minimum amount of memory required to read
!                 a single Cholesky vector in full storage
!
!***********************************************************************

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Fock_util_global, only: Deco, DensityCheck, Lunit
use Data_Structures, only: Deallocate_DT, DSBA_Type, Integer_Pointer, SBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nDen, lOff1, MinMem(nSym)
logical(kind=iwp), intent(in) :: DoCoulomb(nDen), DoExchange(nDen)
real(kind=wp), intent(in) :: FactC(nDen), FactX(nDen)
type(DSBA_Type), intent(in) :: DLT(nDen), DSQ(nDen), MSQ(nDen)
type(DSBA_Type), intent(inout) :: FLT(nDen), FSQ(nDen)
type(Integer_Pointer), intent(in) :: pNocc(nDen)
#include "chotime.fh"
integer(kind=iwp) :: iBatch, iE, irc, iS, iSkip(nSym), iSym, ISYMA, ISYMB, ISYMD, ISYMG, iSymp, ISYMQ, iSymr, iSymr_Occ, iSyms, &
                     iVec, jB, jD, jDen, jE, jjB, jjS, jR, jS, jSR, jSym, JVEC, k, kOcc(nSym), kRdMem, kSym, l, LKV, LVK, LWORK, &
                     MaxSym, Naa, nb, nBatch, NBL, ng, nk, nMax, np, nq, nr, ns, NumB, NumCho(nSym), NumV, nVec
real(kind=wp) :: TC1X1, TC1X2, TC2X1, TC2X2, TCC1, TCC2, tcoul(2), TCR1, TCR2, TCREO1, TCREO2, TCREO3, TCREO4, texch(2), TOTCPU, &
                 TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TW1X1, TW1X2, TW2X1, TW2X2, TWC1, TWC2, TWR1, TWR2, &
                 TWREO1, TWREO2, TWREO3, TWREO4, xf
logical(kind=iwp) :: DoSomeC, DoSomeX, Square
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
character(len=6) :: Fname
type(SBA_Type), target :: Wab, LqJs
real(kind=wp), allocatable :: Dchk(:)
real(kind=wp), pointer :: LrJs(:,:,:) => null(), XdJb(:) => null(), XgJk(:) => null(), XkJs(:) => null(), VJ(:) => null()
real(kind=wp), parameter :: Thr = 1.0e-12_wp
character(len=*), parameter :: BaseNm = 'CHFV', SECNAM = 'CHO_FTWO_MO'

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in SCF-debug
#endif

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero ! time read/rreorder vectors
tcoul(:) = zero ! time for computing Coulomb
texch(:) = zero ! time for computing Exchange

if (DensityCheck) then
  Square = .true.
  xf = Two
  if (DECO) xf = One
  do jDen=1,nDen
    do jSym=1,nSym
      if ((nBas(jSym) /= 0) .and. (pNocc(jDen)%I1(jSym) /= 0)) then
        call mma_allocate(Dchk,nBas(jSym)**2,Label='Dchk')
        call Cho_X_Test(DSQ(jDen)%SB(jSym)%A2,nBas(jSym),Square,MSQ(jDen)%SB(jSym)%A2,pNocc(jDen)%I1(jSym),xf,Dchk,nBas(jSym)**2, &
                        Thr,irc)
        call mma_deallocate(Dchk)
        if (irc == 0) then
          write(u6,*) '*** DENSITY CHECK : OK! *** SYMM= ',jSym
        else
          write(u6,*) '*** DENSITY CHECK FAILED IN SYMMETRY: ',jSym,'FOR THE DENSITY NR. ',jDen
        end if
      end if
    end do
  end do
end if

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

! ********* Data from the Runfile
call get_iarray('NumCho',NumCho,nSym)

! *************** BIG LOOP OVER VECTORS SYMMETRY *****************
do jSym=1,MaxSym

  if (NumCho(jSym) < 1) cycle

  ! Reading of the vectors is done in the full dimension
  Lunit(:) = -1

  ! Open Files
  do ksym=1,nSym
    if (nBas(ksym) /= 0) then
      iSymp = Mul(ksym,jSym)
      if (iSymp >= ksym) then
        Lunit(iSymp) = 7
        write(Fname,'(A4,I1,I1)') BaseNm,iSymp,ksym
        call DANAME_MF_WA(Lunit(iSymp),Fname)
      end if
    end if
  end do

  ! ------------------------------------------------------

  ! SET UP THE READING
  ! ------------------
  call mma_maxDBLE(LWORK)

  if (MinMem(jSym) > 0) then
    nVec = min(LWORK/MinMem(jSym),NumCho(jSym))
  else
    ! ***QUIT*** bad initialization
    write(u6,*) 'Cho_FTwo_MO: bad initialization'
    rc = 99
    call Abend()
    nVec = -9999 ! dummy assignment - avoid compiler warnings
  end if

  !...this is ONLY for debug:
  !nVec = min(nVec,1)

  if (nVec > 0) then
    nBatch = (NumCho(jSym)-1)/nVec+1
  else
    ! ***QUIT*** insufficient memory
    write(u6,*) 'Cho_FTwo_MO: Insufficient memory for batch'
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

    if (iBatch == nBatch) then
      NumV = NumCho(jSym)-nVec*(nBatch-1)
    else
      NumV = nVec
    end if

    iVec = nVec*(iBatch-1)+1

    ! Allocate memory for reading the vectors
    kRdMem = MinMem(jSym)*NumV
    call mma_allocate(Wab%A0,kRdMem,Label='Wab%A0')
    Wab%iSym = JSYM
    Wab%nSym = nSym
    Wab%iCase = 7

    !setup the skipping flags according to # of Occupied
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

      nk = nBas(kSym)
      np = nBas(iSymp)

      if (nk*np <= 0) cycle
      iS = iE+1

      if ((iSymp == ksym) .and. (iSkip(iSymp) /= 0)) then
        NumB = nk*(nk+1)/2
        ! Special trick for the vector L11 ; used to store X(a,Jb)
        if ((ksym == 1) .and. (jSym == 1) .and. DoSomeX) then
          iE = iE+lOff1*NumV
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

#   ifdef _DEBUGPRINT_
    write(u6,*) 'Batch ',iBatch,' of ',nBatch,': NumV = ',NumV
    write(u6,*) 'Total allocated:     ',kRdMem
    write(u6,*) 'iE:              ',iE
    write(u6,*) 'JSYM:                ',jSym
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

        VJ(1:NumV) => Wab%A0(iE+1:iE+NumV)
        VJ(:) = Zero

        do iSymr=1,nSym
          if ((nBas(iSymr) /= 0) .and. (pNocc(jDen)%I1(iSymr) /= 0)) then

            Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

            call DGEMV_('T',Naa,NumV,ONE,Wab%SB(iSymr)%A2,Naa,DLT(jDen)%SB(iSYMR)%A1,1,ONE,VJ,1)

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
    !
    ! WE INCLUDE THE DIAGONAL PART OF THE EXCHANGE HERE!!!!
    !
    ! ********** REORDER and SQUARE ONE VECTOR AT THE TIME *****
    !
    !     CHOVEC(nrs,numv) ---> CHOVEC(nr,numv,ns)

    if ((jSym == 1) .and. DoSomeX) then
      do iSymr=1,nSym

        iSymr_Occ = 0
        do jDen=1,nDen
          iSymr_Occ = iSymr_Occ+pNocc(jDen)%I1(iSymr)
        end do
        nr = NBAS(iSymr)
        if ((nr /= 0) .and. (iSymr_Occ > 0)) then

          call CWTIME(TCREO1,TWREO1)

          LrJs(1:nr,1:NumV,1:nr) => Wab%A0(iE+1:iE+nr*NumV*nr)

          do JVEC=1,NUMV

            do jR=1,NBAS(iSymr)
              do jS=jR,NBAS(iSymr)

                jSR = iTri(jS,jR)

                LrJs(jR,JVEC,jS) = Wab%SB(ISYMR)%A2(jSR,JVEC)

                LrJs(jS,JVEC,jR) = Wab%SB(ISYMR)%A2(jSR,JVEC)

              end do
            end do
          end do

          call CWTIME(TCREO2,TWREO2)
          tread(1) = tread(1)+(TCREO2-TCREO1)
          tread(2) = tread(2)+(TWREO2-TWREO1)

          do jDen=1,nDen

            if (DoExchange(jDen)) then

              kOcc(iSymr) = pNocc(jDen)%I1(iSymr)

              if (kOcc(iSymr) /= 0) then

                call CWTIME(TC1X1,TW1X1)

                ! Calculate intermediate:
                ! X(k,Js) = Sum(r) C(r,k) * L(r,Js).
                ! -----------------------------------
                iSyms = iSymr
                NK = kOcc(iSymr)

                XkJs(1:NK*NumV*nr) => Wab%A0(1:NK*NumV*nr)

                call DGEMM_('T','N',NK,NUMV*NBAS(ISYMR),NBAS(iSYMR),ONE,MSQ(jDen)%SB(ISYMR)%A2,NBAS(ISYMR),LrJs,NBAS(ISYMR),ZERO, &
                            XkJs,NK)

                ! Calculate exchange contribution:
                ! F(r,s) = F(r,s) - FactX * Sum(kJ) X(kJ,r) * X(kJ,s).
                ! ----------------------------------------------------
                !CALL DGEMM_('T','N',NBAS(ISYMR),NBAS(ISYMS),NK*NUMV,-FactX(jDen),XkJs,NK*NUMV,XkJs,NK*NUMV,One, &
                !            FSQ(jDen)%SB(ISYMR)%A2,NBAS(ISYMR))

                ! *** Compute only the LT part of the exchange matrix ***************
                LKV = NK*NUMV
                do jS=1,NBAS(iSymS)
                  NBL = NBAS(iSymR)-(jS-1)
                  jjS = 1+LKV*(jS-1)

                  call DGEMV_('T',LKV,NBL,-FactX(jDen),XkJs(jjS:),LKV,XkJs(jjS:),1,ONE,FSQ(jDen)%SB(ISYMR)%A2(jS:,jS),1)

                end do

                call CWTIME(TC1X2,TW1X2)
                texch(1) = texch(1)+(TC1X2-TC1X1)
                texch(2) = texch(2)+(TW1X2-TW1X1)

                XkJs => null()

              end if ! if kocc /= 0

            end if ! Do Exchange(jDen)

          end do ! loop over the densities
          LrJs => null()

        end if ! nbas /= 0

      end do ! loop over Fock mat symmetries

    end if ! jSym=1 & DoSomeX

    ! ****************** END DIAGONAL EXCHANGE **********************

    ! ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************

    if (jSym /= 1) then

      ! Reordering of the Cholesky vectors (ALL AT ONCE)
      !
      !          Reorder: L(qs,J) -> L(qJ,s).
      !
      ! CHOVEC(nq,ns,numv) ---> CHOVEC(nq,numv,ns)
      ! ------------------------------------------
      call CWTIME(TCREO3,TWREO3)

      jE = iE
      do ISYMS=1,NSYM

        ISYMQ = Mul(ISYMS,JSYM)

        ns = nBas(iSyms)
        nq = nBas(iSymq)
        if (ns*nq <= 0) cycle

        if ((ISYMQ > ISYMS) .and. (iSkip(iSymq) /= 0)) then

          jS = JE+1
          jE = JE+nq*NumV*ns
          LqJs%SB(ISYMQ)%A3(1:nq,1:NumV,1:ns) => Wab%A0(jS:jE)

          do JVEC=1,NUMV

            do jS=1,NBAS(iSyms)

              LqJs%SB(ISYMQ)%A3(:,jVec,jS) = Wab%SB(ISYMQ)%A3(:,jS,jVEC)

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
          ! **** Occupation numbers needed for the Exchange term *****
          do iSym=1,nsym
            kOcc(iSym) = pNocc(jDen)%I1(iSym)
          end do
          ! **********************************************************
          ! --- COMPUTE EXCHANGE FOR OFF-DIAGONAL VECTORS

          call CWTIME(TC2X1,TW2X1)

          do ISYMG=1,NSYM

            ISYMB = Mul(ISYMG,JSYM)

            nb = nBas(iSymb)
            ng = nBas(iSymg)
            if (nb*ng <= 0) cycle

            ISYMA = ISYMB ! block diagonal Fock Matrix
            ISYMD = ISYMG ! total symmetric densities

            if ((ISYMG > ISYMB) .and. (iSkip(iSymg) /= 0)) then

              ! ---------------------------
              ! F(a,b) = - D(g,d) * (ad|gb)
              ! ---------------------------
              NK = kOcc(iSymG)

              XdJb(1:NK*NumV*nb) => Wab%A0(1:NK*NumV*nb)

              if (NK /= 0) then

                ! Calculate intermediate:
                ! X(k,Jb) = Sum(g) C(g,k) * L(g,Jb).
                ! ----------------------------------

                call DGEMM_('T','N',NK,NUMV*NBAS(ISYMB),NBAS(ISYMG),ONE,MSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMG),LqJs%SB(ISYMG)%A3, &
                            NBAS(ISYMG),ZERO,XdJb,NK)

                ! F(a,b) = F(a,b) - Sum(kJ) X(kJ,a) * X(kJ,b).
                ! -------------------------------------------

                !CALL DGEMM_('T','N',NBAS(ISYMA),NBAS(ISYMB),NK*NUMV,-FactX(jDen),XdJb,NK*NUMV,XdJb,NK*NUMV,ONE, &
                !            FSQ(jDen)%SB(iSYMB)%A2,NBAS(ISYMA))

                ! *** Compute only the LT part of the exchange matrix ***************
                LKV = NK*NUMV
                do jB=1,NBAS(iSymB)
                  NBL = NBAS(iSymA)-(jB-1)
                  jjB = 1+LKV*(jB-1)

                  call DGEMV_('T',LKV,NBL,-FactX(jDen),XdJb(jjB:),LKV,XdJb(jjB:),1,ONE,FSQ(jDEN)%SB(ISYMB)%A2(jB:,jB),1)

                end do
                XdJb => null()
                ! ******************************************************************

              end if

              ! ---------------------------
              ! F(g,d) = - D(a,b) * (ad|gb)
              ! ---------------------------
              NK = kOcc(iSymB)

              if (NK /= 0) then

                XgJk(1:ng*NumV*NK) => Wab%A0(1:ng*NumV*NK)

                ! Calculate intermediate:
                ! X(gJ,k) = Sum(a) L(gJ,a) * C(a,k).
                ! ----------------------------------

                call DGEMM_('N','N',NBAS(ISYMG)*NUMV,NK,NBAS(ISYMA),ONE,LqJs%SB(ISYMG)%A3,NBAS(ISYMG)*NUMV,MSQ(jDen)%SB(ISYMB)%A2, &
                            NBAS(ISYMA),ZERO,XgJk,NBAS(ISYMG)*NUMV)

                ! F(g,d) = F(g,d) - Sum(Jk) X(g,Jk) * X(d,Jk).
                ! -------------------------------------------

                !CALL DGEMM_('N','T',NBAS(ISYMG),NBAS(ISYMD),NUMV*NK,-FactX(jDen),XgJk,NBAS(ISYMG),XgJk,NBAS(ISYMD),ONE, &
                !            FSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMG))

                ! *** Compute only the LT part of the exchange matrix ***************
                LVK = NUMV*NK
                do jD=1,NBAS(iSymD)
                  NBL = NBAS(iSymG)-(jD-1)

                  call DGEMV_('N',NBL,LVK,-FactX(jDen),XgJk(jD:),NBAS(iSymG),XgJk(jD:),NBAS(iSymD),ONE, &
                              FSQ(jDen)%SB(ISYMG)%A2(jD:,jD),1)
                end do

                XgJk => null()
              end if

            end if
            LqJs%SB(ISYMG)%A3 => null()

          end do ! loop over orbital symmetries

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

  write(u6,'(6X,A)') 'TEST PRINT FROM CHO_FTWO_MO.'
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

end subroutine CHO_FTWO_MO
