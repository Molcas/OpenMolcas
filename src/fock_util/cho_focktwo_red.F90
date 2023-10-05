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

subroutine CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
!***********************************************************************
!  Author : F. Aquilante
!
!  Purpose:
!          For a given set of total symmetric AO-densities
!          (hermitian matrix representation) this routine
!          computes the Coulomb and Exchange contributions
!          to the Fock matrix from Cholesky vectors in full storage
!
!          The vectors are indeed read in reduced set storage
!          and then reordered.
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

use Cholesky, only: nSym, NumCho, timings
use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Data_Structures, only: Deallocate_DT, DSBA_type, Integer_Pointer, SBA_type
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: nBas(*), nDen, MinMem(*)
logical(kind=iwp), intent(in) :: DoCoulomb(nDen), DoExchange(nDen)
real(kind=wp), intent(in) :: FactC(nDen), FactX(nDen)
type(DSBA_Type), intent(in) :: DLT(nDen), DSQ(nDen)
type(DSBA_Type), intent(inout) :: FLT(nDen), FSQ(nDen)
type(Integer_Pointer), intent(in) :: pNocc(nDen)
integer(kind=iwp) :: iBatch, iE, irc, IREDC, iS, iSkip(8), iSwap, iSym, ISYMA, ISYMB, ISYMD, ISYMG, iSymp, iSymq, iSymr, &
                     iSymr_Occ, iSyms, iVec, jD, jDen, jR, jS, jSR, jSym, JVEC, k, kTOT, l, lScr, LWORK, MaxSym, Naa, nB, nBatch, &
                     nD, nG, Nmax, np, nq, nr, NumB, NumV, nVec
real(kind=wp) :: TC1X1, TC1X2, TC2X1, TC2X2, TCC1, TCC2, tcoul(2), TCR1, TCR2, TCREO1, TCREO2, texch(2), TOTCPU, TOTCPU1, TOTCPU2, &
                 TOTWALL, TOTWALL1, TOTWALL2, tread(2), TW1X1, TW1X2, TW2X1, TW2X2, TWC1, TWC2, TWR1, TWR2, TWREO1, TWREO2
logical(kind=iwp) :: DoSomeC, DoSomeX
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
type(SBA_Type), target :: Wab
real(kind=wp), pointer :: LrJs(:,:,:) => null(), Scr(:) => null(), VJ(:) => null(), XdJb(:,:,:) => null(), XgJb(:,:,:) => null(), &
                          XpJs(:,:,:) => null()
logical(kind=iwp), parameter :: DoRead = .true.
character(len=*), parameter :: SECNAM = 'CHO_FOCKTWO_RED'

!*************************************************
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in SCF-debug
#endif
IREDC = -1 ! unknown reduced set in core

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

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

! *************** BIG LOOP OVER VECTORS SYMMETRY *****************
do jSym=1,MaxSym

  if (NumCho(jSym) < 1) cycle

  ! ------------------------------------------------------

  ! SET UP THE READING
  ! ------------------
  call mma_maxDBLE(LWORK)

  if (MinMem(jSym) > 0) then
    nVec = min(LWORK/MinMem(jSym),NumCho(jSym))
  else
    ! ***QUIT*** bad initialization
    write(u6,*) 'Cho_FockTwo_RED: bad initialization'
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
    write(u6,*) 'Cho_FockTwo_RED: Insufficient memory for batch'
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
    kTOT = MinMem(jSym)*NumV
    call mma_allocate(Wab%A0,kTOT,Label='Wab%A0')
    Wab%iSym = JSYM
    Wab%nSym = nSym
    Wab%iCase = 7

    ! setup the skipping flags according to # Occupied
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

    ! Max dimension of a read symmetry block
    Nmax = 0
    do iSym=1,nSym
      if ((NBAS(iSym) > Nmax) .and. (iSkip(iSym) /= 0)) then
        Nmax = NBAS(iSym)
      end if
    end do

    ! Compute the total amount of memory needed to have the
    ! vectors in core (full storage)
    iE = 0
    do iSymq=1,nSym
      iSymp = Mul(jSym,iSymq)
      nq = nBas(iSymq)
      np = nBas(iSymp)
      if (nq*np <= 0) cycle

      iS = iE+1
      if ((iSymp > iSymq) .and. (iSkip(iSymp) /= 0)) then
        NumB = np*nq
        iE = iE+NumB*NumV
        Wab%SB(iSymp)%A3(1:np,1:nq,1:NumV) => Wab%A0(iS:iE)
        Wab%SB(iSymp)%A2(1:NumB,1:NumV) => Wab%A0(iS:iE)
      else
        if ((iSymp == iSymq) .and. (iSkip(iSymp) /= 0)) then
          NumB = np*(np+1)/2
          ! Special trick for the vector L11 ; used to store X(a,Jb)
          if ((iSymp == 1) .and. (jSym == 1) .and. DoSomeX) then
            iE = iE+(Nmax**2)*NumV
            Wab%SB(iSymp)%A1(1:Nmax**2*NumV) => Wab%A0(iS:iE)
          else
            iE = iE+NumB*NumV
          end if
          Wab%SB(iSymp)%A2(1:NumB,1:NumV) => Wab%A0(iS:iS-1+NumB*NumV)
        end if
      end if
    end do

    lScr = kTOT-iE

    ! Reading of the vectors is done in Reduced sets

    iSwap = min(1,(jSym-1)) ! L(ab,J) --> L(a,J,b) iff jSym /= 1

    call CWTIME(TCR1,TWR1)

    Scr(1:lScr) => Wab%A0(iE+1:iE+lScr)
    call CHO_X_getVfull(irc,Scr,lscr,iVEC,NumV,jSym,iSwap,IREDC,Wab,iSkip,DoRead)
    Scr => null()

    if (irc /= 0) then
      rc = irc
      return
    end if

    call CWTIME(TCR2,TWR2)
    tread(1) = tread(1)+(TCR2-TCR1)
    tread(2) = tread(2)+(TWR2-TWR1)

#   ifdef _DEBUGPRINT_
    write(u6,*) 'Batch ',iBatch,' of   ',nBatch,': NumV = ',NumV
    write(u6,*) 'Total allocated :     ',kTOT
    write(u6,*) 'lScr:                 ',lScr
    write(u6,*) 'JSYM:                 ',jSym
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
        !------------------------------------
        ! V{#J} <- V{#J} + sum_rr L(rr,{#J})*D(rr)
        !==========================================================

        call CWTIME(TCC1,TWC1)

        VJ(1:NumV) => Wab%A0(iE+1:iE+NumV)
        VJ(:) = Zero

        do iSymr=1,nSym
          if ((nBas(iSymr) /= 0) .and. (pNocc(jDen)%I1(iSymr) /= 0)) then

            Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

            call DGEMV_('T',Naa,NumV,ONE,Wab%SB(iSymr)%A2,Naa,DLT(jDen)%SB(ISYMR)%A1,1,ONE,VJ,1)

          end if
        end do

        ! Contraction with the 2nd vector
        ! -----------------------------------
        ! Fpp{#J} <- Fpp{#J} + FactC * sum_J L(pp,{#J})*V{#J}
        !==========================================================

        do iSyms=1,nSym
          if (nBas(iSyms) /= 0) then

            Naa = nBas(iSyms)*(nBas(iSyms)+1)/2

            call DGEMV_('N',Naa,NumV,FactC(jDen),Wab%SB(iSyms)%A2,Naa,VJ,1,ONE,FLT(jDen)%SB(ISYMS)%A1,1)

          end if
        end do

        call CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1)+(TCC2-TCC1)
        tcoul(2) = tcoul(2)+(TWC2-TWC1)

        VJ => null()

      end if !jSym=1 & DoCoulomb

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

          nr = NBAS(iSymr)
          LrJs(1:nr,1:NUMV,1:nr) => Wab%A0(iE+1:iE+nr*NUMV*nr)

          do JVEC=1,NUMV

            do jR=1,NBAS(iSymr)
              do jS=jR,NBAS(iSymr)

                jSR = iTri(jS,jR)

                LrJs(jR,JVEC,jS) = Wab%SB(ISYMR)%A2(jSR,jVEC)

                LrJs(jS,JVEC,jR) = Wab%SB(ISYMR)%A2(jSR,jVEC)

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

                ! Calculate intermediate:
                ! X(p,Js) = Sum(q) D(p,q) * L(q,Js).
                ! -----------------------------------

                NR = NBAS(ISYMR)
                XpJs(1:NR,1:NUMV,1:NR) => Wab%A0(1:NR*NUMV*NR)

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

              end if

            end if ! DoExchange(jDen)

          end do ! loop over the densities
          LrJs => null()

        end if ! nbas /= 0 & nOcc /= 0

      end do ! loop over Fock mat symmetries

    end if ! jSym=1 and DoSomeX

    ! ****************** END DIAGONAL EXCHANGE **********************
    !
    ! ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************

    if (jSym /= 1) then

      do jDen=1,nDen

        if (DoExchange(jDen)) then
          ! COMPUTE EXCHANGE FOR OFF-DIAGONAL VECTORS

          call CWTIME(TC2X1,TW2X1)

          do ISYMG=1,NSYM

            ISYMB = Mul(ISYMG,JSYM)

            if (nBas(iSymb)*nBas(iSymg) <= 0) cycle

            ISYMD = ISYMG ! only total symmetric Density
            ISYMA = ISYMB ! block diagonal Fock Matrix

            if ((ISYMG > ISYMB) .and. (iSkip(iSymg) /= 0)) then

              ! ---------------------------
              ! F(a,b) = - D(g,d) * (ad|gb)
              ! ---------------------------
              if (pNocc(jDen)%I1(iSymg) /= 0) then

                nD = NBAS(ISYMD)
                nB = NBAS(ISYMB)
                XdJb(1:nD,1:NumV,1:nB) => Wab%A0(iE+1:iE+nD*NumV*nB)

                ! Calculate intermediate:
                ! X(d,Jb) = Sum(g) D(d,g) * L(g,Jb).
                ! ----------------------------------

                call DGEMM_('N','N',NBAS(ISYMD),NUMV*NBAS(ISYMB),NBAS(ISYMG),ONE,DSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMD), &
                            Wab%SB(ISYMG)%A2,NBAS(ISYMG),ZERO,XdJb,NBAS(ISYMD))

                ! F(a,b) = F(a,b) - Sum(dJ) L(dJ,a) * X(dJ,b).
                ! -------------------------------------------

                call DGEMM_('T','N',NBAS(ISYMA),NBAS(ISYMB),NBAS(ISYMD)*NUMV,-FactX(jDen),Wab%SB(ISYMG)%A2,NBAS(ISYMD)*NUMV,XdJb, &
                            NBAS(ISYMD)*NUMV,ONE,FSQ(jDen)%SB(ISYMB)%A2,NBAS(ISYMA))

                XdJb => null()

              end if
              ! ---------------------------
              ! F(g,d) = - D(a,b) * (ad|gb)
              ! ---------------------------
              if (pNocc(jDen)%I1(iSyma) /= 0) then

                nG = NBAS(ISYMG)
                nB = NBAS(ISYMB)
                XgJb(1:nG,1:NumV,1:nB) => Wab%A0(iE+1:iE+nG*NumV*nB)

                ! Calculate intermediate:
                ! X(gJ,b) = Sum(a) L(gJ,a)* D(a,b).
                ! ----------------------------------

                call DGEMM_('N','N',NBAS(ISYMG)*NUMV,NBAS(ISYMB),NBAS(ISYMA),ONE,Wab%SB(ISYMG)%A2,NBAS(ISYMG)*NUMV, &
                            DSQ(jDen)%SB(ISYMB)%A2,NBAS(ISYMA),ZERO,XgJb,NBAS(ISYMG)*NUMV)

                ! F(g,d) = F(g,d) - Sum(Jb) X(g,Jb) * L(d,Jb).
                ! -------------------------------------------

                call DGEMM_('N','T',NBAS(ISYMG),NBAS(ISYMD),NUMV*NBAS(ISYMB),-FactX(jDen),XgJb,NBAS(ISYMG),Wab%SB(ISYMG)%A2, &
                            NBAS(ISYMD),ONE,FSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMG))

                XgJb => null()

              end if
            end if

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

end do ! loop over JSYM

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

  write(u6,'(6X,A)') 'TEST PRINT FROM CHO_FOCKTWO_RED.'
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

end subroutine CHO_FOCKTWO_RED
