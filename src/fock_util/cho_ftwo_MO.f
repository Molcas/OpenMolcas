************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,
     &                       lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,
     &                       MinMem,MSQ,pNocc)

************************************************************************
*  Author : F. Aquilante
*
*  Purpose:
*          For a given set of total symmetric AO-densities
*          (hermitian matrix representation) and MO coefficients
*          this routine
*          computes the Coulomb and Exchange contributions
*          to the Fock matrix from Cholesky vectors in full storage
*
*          Coulomb term (Drs is the lower triangular matrix
*                        of r=s symmetry block. The diagonals
*                        are scaled by a factor 1/2 already in
*                        the packed density DLT, given as input)
*
*    F(pq) = F(pq) + FactC * Sum_J Lpq,J * (Sum_rs Drs*Lrs,J)
*
*
*          Exchange term (Crk is the squared matrix of MO coeff.
*                             of r=s symmetry block)
*
*    F(pq) = F(pq) - FactX * Sum_Jk [(Sum_r Crk * Lpr,J) * (Sum_s Csk * Lqs,J)]
*
*
*    Integral-fashion formula:
*
*    F(pq) =  Sum_rs Drs * [(pq|rs)  - 0.5 (ps|rq)]
*
*----------------------------------------------------------------------
*  DoCoulomb(nDen) :  specifies the densities for which coulomb term
*                     has to be computed
*
*  DoExchange(nDen) : specifies the densities for which exchange term
*                     has to be computed
*
*  FactC(nDen) : factor for the coulomb of the corresponding density
*  FactX(nDen) : factor for the exchange of the corresponding density
*
*  ip{X}LT(nDen) : pointer to the array containing {X} in LT storage
*  ip{X}SQ(nDen) : pointer to the array containing {X} in SQ storage
*    {X=D,F,M --- Density, Fock matrix, MOs coeff.}
*
*  MinMem(nSym) : minimum amount of memory required to read
*                 a single Cholesky vector in full storage
*
************************************************************************
      use Data_Structures, only: SBA_Type, Deallocate_SBA
      use Data_Structures, only: DSBA_Type
      Implicit Real*8 (a-h,o-z)

      Integer   rc,nDen,nSym,nBas(nSym),NumCho(nSym),kOcc(nSym)
      Real*8    FactC(nDen),FactX(nDen)
      Integer   Lunit,ISTSQ(nSym),ISTLT(nSym),lOff1
      Integer   MinMem(nSym),iSkip(nSym)
      Integer   ipDLT(2),ipDSQ(3),ipFLT(2),ipFSQ(3), ipMSQ(3)

      Type (DSBA_Type) DLT(nDen), FLT(nDen), FSQ(nDen), DSQ(nDen),
     &                 MSQ(nDen)
      Type Integer_Pointer
          Integer, Pointer :: I1(:)=>Null()
      End Type Integer_Pointer
      Type (Integer_Pointer) :: pNocc(nDen)


      Real*8    tread(2),tcoul(2),texch(2)
#include "chounit.fh"
      Logical DoExchange(nDen),DoCoulomb(nDen),DoSomeX,DoSomeC
#ifdef _DEBUGPRINT_
      Logical Debug
#endif
      Logical Square
#include "chotime.fh"
#include "chodensity.fh"
#include "choscf.fh"
#include "real.fh"
      Character*50 CFmt
      Character(LEN=11), Parameter :: SECNAM = 'CHO_FTWO_MO'

      Character(LEN=4), Parameter :: BaseNm = 'CHFV'
      Character*6 Fname

      Real*8, Parameter :: xone = -One

#include "WrkSpc.fh"
#include "stdalloc.fh"

      Type (SBA_Type) , Target:: Wab, LqJs

      Real*8, Allocatable :: Dchk(:)

      Real*8, Pointer :: VJ(:)=>Null(), LrJs(:,:,:)=>Null(),
     &                    XkJs(:)=>Null(), XdJb(:)=>Null(),
     &                    XgJk(:)=>Null()
*                                                                      *
************************************************************************
*                                                                      *
      Interface
        subroutine dgemv_(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
          Character(LEN=1) TRANS
          Integer M, N
          Real*8 ALPHA, BETA
          Integer LDA, INCX, INCY
          Real*8  A(lda,*), X(*), Y(*)
        End subroutine dgemv_
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      nOcc(jSym,jDen) = pNocc(jDen)%I1(jSym)
**************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in SCF-debug
#endif
      Do iDen = 1, nDen
         ipDLT(iDen) = ip_of_Work(DLT(iDen)%A0(1))
         ipFLT(iDen) = ip_of_Work(FLT(iDen)%A0(1))
         ipFSQ(iDen) = ip_of_Work(FSQ(iDen)%A0(1))
         ipDSQ(iDen) = ip_of_Work(DSQ(iDen)%A0(1))
         ipMSQ(iDen) = ip_of_Work(MSQ(iDen)%A0(1))
      End Do

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      ! 1 --> CPU   2 --> Wall
      tread(:) = zero  !time read/rreorder vectors
      tcoul(:) = zero  !time for computing Coulomb
      texch(:) = zero  !time for computing Exchange

c ISTSQ: Offsets to full symmetry block in DSQ,FSQ
c ISTLT: Offsets to packed LT symmetry blocks in DLT,FLT
      ISTSQ(1)=0
      ISTLT(1)=0
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NB2=NB*NB
        NB3=(NB2+NB)/2
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NB2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NB3
      END DO

      if(DensityCheck)then
        Thr=1.0d-12
        Square=.true.
        xf=2.0d0
        if (DECO) xf=1.0d0
        Do jDen=1,nDen
          Do jSym=1,nSym
           if (nBas(jSym).ne.0.and.nOcc(jSym,jDen).ne.0) then
             Call mma_allocate(Dchk,nBas(jSym)**2,Label='Dchk')
             Call Cho_X_Test(Work(ipDSQ(jDen)+ISTSQ(jSym)),nBas(jSym),
     &                       Square,Work(ipMSQ(jDen)+ISTSQ(jSym)),
     &                       nOcc(jSym,jDen),xf,Dchk,
     &                       nBas(jSym)**2,Thr,irc)
             Call mma_deallocate(Dchk)
             if(irc.eq.0)then
                write(6,*)'*** DENSITY CHECK : OK! *** SYMM= ',jSym
             else
                write(6,*)'*** DENSITY CHECK FAILED IN SYMMETRY: ',
     &                    jSym,'FOR THE DENSITY NR. ',jDen
             endif
           endif
          End Do
        End Do
      endif


C --- Tests on the type of calculation
        DoSomeC=.false.
        DoSomeX=.false.
        MaxSym=0
c
        jD=0
        do while (jD.lt.nDen .and. .not.DoSomeX)
         jD=jD+1
         DoSomeX=DoExchange(jD)
        end do
        jD=0
        do while (jD.lt.nDen .and. .not.DoSomeC)
         jD=jD+1
         DoSomeC=DoCoulomb(jD)
        end do

        if (DoSomeX) Then
           MaxSym=nSym !we want to do some exchange
        Else
           If (DoSomeC) Then
                   MaxSym = 1 !we want to do only Coulomb
           Else
                   Return  !we are lazy and won''t do anything
           End If
        End If

C ********* Data from the Runfile
      Call get_iarray('NumCho',NumCho,nSym)

C *************** BIG LOOP OVER VECTORS SYMMETRY *****************
      DO jSym=1,MaxSym

         If (NumCho(jSym).lt.1) Cycle

C --- Reading of the vectors is done in the full dimension
         Lunit(:) = -1

C Open Files
       Do ksym=1,nSym
        if (nBas(ksym).ne.0) then
          iSymp=MulD2h(ksym,jSym)
          If (iSymp.ge.ksym) then
             Lunit(iSymp) = 7
             Write(Fname,'(A4,I1,I1)') BaseNm,iSymp,ksym
             Call DANAME_MF_WA(Lunit(iSymp),Fname)
          End If
        End If
       End Do

C ------------------------------------------------------


C SET UP THE READING
C ------------------
      Call mma_maxDBLE(LWORK)

      If (MinMem(jSym) .gt. 0) Then
         nVec = Min(LWORK/MinMem(jSym),NumCho(jSym))
      Else
C         ***QUIT*** bad initialization
         WRITE(6,*) 'Cho_FTwo_MO: bad initialization'
         rc=99
         CALL Abend()
         nVec = -9999  ! dummy assignment - avoid compiler warnings
      End If

C...this is ONLY for debug:
c      nVec = min(nVec,1)

      If (nVec .gt. 0) Then
         nBatch = (NumCho(jSym) - 1)/nVec + 1
      Else
C         ***QUIT*** insufficient memory
         WRITE(6,*) 'Cho_FTwo_MO: Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'min. mem. need= ',MinMem(jSym)
         WRITE(6,*) 'NumCho= ',NumCho(jsym)
         WRITE(6,*) 'jsym= ',jsym
         rc = 205
         CALL Abend()
         nBatch = -9999  ! dummy assignment
      End If


C *************** BATCHING  *****************

      DO iBatch = 1,nBatch

         If (iBatch .eq. nBatch) Then
            NumV = NumCho(jSym) - nVec*(nBatch - 1)
         Else
            NumV = nVec
         End If

         iVec = nVec*(iBatch - 1) + 1

C --- Allocate memory for reading the vectors
         kRdMem = MinMem(jSym)*NumV
         Call mma_allocate(Wab%A0,kRdMem,Label='Wab%A0')
         Wab%iSym=JSYM
         Wab%nSym=nSym
         Wab%iCase=7

c--- setup the skipping flags according to # of Occupied
      do k=1,nSym
         iSkip(k)=0
         l=Muld2h(k,jsym) !L(kl) returned if nOcc(k or l) .ne.0
         if (k.eq.l) then
            iSkip(k) = 666 ! always contribute to Coulomb
         else
          do jDen=1,nDen
               iSkip(k) = iSkip(k) + nOcc(k,jDen) + nOcc(l,jDen)
          end do
         endif
      end do

C Max dimension of a symmetry block
        Nmax=0
        Do iSym=1,nSym
           if(NBAS(iSym).gt.Nmax .and. iSkip(iSym).ne.0)then
           Nmax = NBAS(iSym)
           endif
        End Do

       CALL CWTIME(TCR1,TWR1)

       iE=0
       Do ksym=1,nSym

         iSymp=MulD2h(ksym,jSym)

         nk=nBas(kSym)
         np=nBas(iSymp)

         IF (nk*np<=0) Cycle
         iS = iE + 1

         If (iSymp.eq.ksym .and. iSkip(iSymp).ne.0) then
           NumB = nk*(nk+1)/2
C --- Special trick for the vector L11 ; used to store X(a,Jb)
           if(ksym.eq.1.and.jSym.eq.1.and.DoSomeX)then
              iE = iE + lOff1 * NumV
           else
              iE = iE + NumB * NumV
           endif
           Wab%SB(iSymp)%A2(1:NumB,1:NumV) => Wab%A0(iS:iS-1+NumB*NumV)
         Else
           If (iSymp.gt.ksym .and. iSkip(iSymp).ne.0) then
              NumB = nk*np
              iE = iE + NumB * NumV
              Wab%SB(iSymp)%A3(1:np,1:nk,1:NumV) => Wab%A0(iS:iE)
              Wab%SB(iSymp)%A2(1:NumB,1:NumV) => Wab%A0(iS:iE)
           End If
         End If

       End Do ! ends the loop over symmetries

       Do kSym = 1, nSym
         iSymp=MulD2h(ksym,jSym)
         If (.NOT.Associated(Wab%SB(iSymp)%A2)) Cycle
         NumB = SIZE(Wab%SB(iSymp)%A2,1)
         Call RdChoVec(Wab%SB(iSymp)%A2,NumB,NumV,iVec,Lunit(iSymp))
       End Do


       CALL CWTIME(TCR2,TWR2)
       tread(1) = tread(1) + (TCR2 - TCR1)
       tread(2) = tread(2) + (TWR2 - TWR1)

#ifdef _DEBUGPRINT_
       write(6,*) 'Batch ',iBatch,' of ',nBatch,': NumV = ',NumV
       write(6,*) 'Total allocated:     ',kRdMem
       write(6,*) 'iE:              ',iE
       write(6,*) 'JSYM:                ',jSym
#endif

C
C --- Reading Done!
C
C ********************************************************************
      DO jDen=1,nDen

       IF (jSym.eq.1.and.DoCoulomb(jDen)) THEN
C
C *************** BEGIN COULOMB CONTRIBUTION  ****************
C
C---- Computing the intermediate vector V(J)
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J} + sum_rr L(rr,{#J})*D(rr)
C==========================================================
C
         CALL CWTIME(TCC1,TWC1)

         VJ(1:NumV) => Wab%A0(iE+1:iE+NumV)
         VJ(:) = Zero

         DO iSymr=1,nSym
         IF(nBas(iSymr).ne.0.and.nOcc(iSymr,jDen).ne.0)THEN

         ISDLT = ipDLT(jDen) + ISTLT(ISYMR)
         Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

         CALL DGEMV_('T',Naa,NumV,
     &              ONE,Wab%SB(iSymr)%A2,Naa,
     &              Work(ISDLT),1,ONE,VJ,1)

         ENDIF
         End DO

C --- Contraction with the 2nd vector
C ---------------------------------------
C --- Fss{#J} <- Fss{#J} + FactC * sum_J L(ss,{#J})*V{#J}
C==========================================================
C
          DO iSyms=1,nSym
          IF(nBas(iSyms).ne.0)THEN

           ISFLT = ipFLT(jDen) + ISTLT(ISYMS)
           Naa = nBas(iSyms)*(nBas(iSyms)+1)/2

           CALL DGEMV_('N',Naa,NumV,
     &              FactC(jDen),Wab%SB(iSyms)%A2,Naa,
     &              VJ,1,ONE,Work(ISFLT),1)

c           WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMS
c           CALL TRIPRT('Coulomb FLT',' ',Work(ISFLT),nBas(iSyms))
          ENDIF
          End DO

          VJ => Null()

        CALL CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1) + (TCC2 - TCC1)
        tcoul(2) = tcoul(2) + (TWC2 - TWC1)

       ENDIF  ! jSym=1 & DoCoulomb

      END DO ! loop over densities for the Coulomb term

C *************** END COULOMB CONTRIBUTION  ****************
C
C WE INCLUDE THE DIAGONAL PART OF THE EXCHANGE HERE!!!!
C
C ********** REORDER and SQUARE ONE VECTOR AT THE TIME *****
C
C     CHOVEC(nrs,numv) ---> CHOVEC(nr,numv,ns)

       IF (jSym.eq.1.and.DoSomeX) THEN
         ISFSQ=0
         ISMSQ=0
         DO iSymr=1,nSym

          iSymr_Occ=0
          do jDen=1,nDen
             iSymr_Occ = iSymr_Occ + nOcc(iSymr,jDen)
          end do
          nr=NBAS(iSymr)
          IF(nr.ne.0.and.iSymr_Occ.gt.0)THEN

             CALL CWTIME(TCREO1,TWREO1)

             LrJs(1:nr,1:NumV,1:nr) => Wab%A0(iE+1:iE+nr*NumV*nr)

             Do JVEC=1,NUMV

                Do jR=1,NBAS(iSymr)
                   Do jS=jR,NBAS(iSymr)

                      jSR = iTri(jS,jR)

                      LrJs(jR,JVEC,jS) = Wab%SB(ISYMR)%A2(jSR,JVEC)

                      LrJs(jS,JVEC,jR) = Wab%SB(ISYMR)%A2(jSR,JVEC)

                   End Do
                End Do
             End Do

             CALL CWTIME(TCREO2,TWREO2)
             tread(1) = tread(1) + (TCREO2 - TCREO1)
             tread(2) = tread(2) + (TWREO2 - TWREO1)

        DO jDen=1,nDen

           IF (DoExchange(jDen)) THEN

              kOcc(iSymr) = nOcc(iSymr,jDen)

           IF (kOcc(iSymr).ne.0) THEN

               CALL CWTIME(TC1X1,TW1X1)

C              Calculate intermediate:
C              X(k,Js) = Sum(r) C(r,k) * L(r,Js).
C              -----------------------------------
               iSyms=iSymr
               ISFSQ = ISTSQ(ISYMR) + ipFSQ(jDen)
               ISMSQ = ISTSQ(ISYMR) + ipMSQ(jDen)
               NK = kOcc(iSymr)

               XkJs(1:NK*NumV*nr) => Wab%A0(1:NK*NumV*nr)

               CALL DGEMM_('T','N',
     &               NK,NUMV*NBAS(ISYMR),NBAS(iSYMR),
     &               ONE,WORK(ISMSQ),NBAS(ISYMR),
     &                   LrJs, NBAS(ISYMR),
     &               ZERO,XkJs,NK)

C              Calculate exchange contribution:
C              F(r,s) = F(r,s) - FactX * Sum(kJ) X(kJ,r) * X(kJ,s).
C              ----------------------------------------------------
c               CALL DGEMM_('T','N',
c     &                NBAS(ISYMR),NBAS(ISYMS),NK*NUMV,
c     &             -FactX(jDen),XkJs,NK*NUMV,
c     &                          XkJs,NK*NUMV,
c     &                      One,Work(ISFSQ),NBAS(ISYMR))

c *** Compute only the LT part of the exchange matrix ***************
               ipF=0
               LKV=NK*NUMV
               DO jS=1,NBAS(iSymS)
                     NBL = NBAS(iSymR) - (jS-1)
                     ipF = ISFSQ + NBAS(iSymR)*(jS-1) + (jS-1)
                     jjS = 1 + LKV*(jS-1)

                     CALL DGEMV_('T',LKV,NBL,
     &                    -FactX(jDen),XkJs(jjS:),LKV,
     &                                 XkJs(jjS:),1,
     &                     ONE,Work(ipF),1)

               END DO

             CALL CWTIME(TC1X2,TW1X2)
             texch(1) = texch(1) + (TC1X2 - TC1X1)
             texch(2) = texch(2) + (TW1X2 - TW1X1)

             XkJs=>Null()

         ENDIF  ! if kocc.ne.0

       ENDIF  ! Do Exchange(jDen)

      END DO  ! loop over the densities
      LrJs => Null()

      ENDIF   ! nbas.ne.0

c         write(6,*)'Symmetry block of FSQ= ',isymr
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMR),NBAS(ISYMR))

      END DO  !loop over Fock mat symmetries

      END IF  ! jSym=1 & DoSomeX

C ****************** END DIAGONAL EXCHANGE **********************
C
C ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************
C
       IF (jSym.ne.1) THEN

C
C  Reordering of the Cholesky vectors (ALL AT ONCE)
C
C              Reorder: L(qs,J) -> L(qJ,s).
C
C     CHOVEC(nq,ns,numv) ---> CHOVEC(nq,numv,ns)
C     ------------------------------------------
            CALL CWTIME(TCREO3,TWREO3)

            jE = iE
            DO ISYMS = 1,NSYM

               ISYMQ = MULD2H(ISYMS,JSYM)

            ns=nBas(iSyms)
            nq=nBas(iSymq)
            IF(ns*nq<=0) Cycle

             IF (ISYMQ.gt.ISYMS .and. iSkip(iSymq).ne.0) THEN

                  jS = JE + 1
                  jE = JE + nq*NumV*ns
                  LqJs%SB(ISYMQ)%A3(1:nq,1:NumV,1:ns) => Wab%A0(jS:jE)

              DO JVEC = 1,NUMV

                DO jS=1,NBAS(iSyms)

                  LqJs%SB(ISYMQ)%A3(:,jVec,jS) =
     &               Wab%SB(ISYMQ)%A3(:,jS,jVEC)

                END DO

              END DO   !loop over the vectors

             END IF  ! sym(Q) > sym(S)

            END DO  !loop over symmetry blocks

            CALL CWTIME(TCREO4,TWREO4)
            tread(1) = tread(1) + (TCREO4 - TCREO3)
            tread(2) = tread(2) + (TWREO4 - TWREO3)

C     -------------- REORDERING DONE !! --------

       DO jDen=1,nDen

       IF (DoExchange(jDen)) THEN
C **** Occupation numbers needed for the Exchange term *****
              do iSym=1,nsym
                 kOcc(iSym) = nOcc(iSym,jDen)
              end do
C **********************************************************
C --- COMPUTE EXCHANGE FOR OFF-DIAGONAL VECTORS

            CALL CWTIME(TC2X1,TW2X1)

            DO ISYMG = 1,NSYM

               ISYMB = MULD2H(ISYMG,JSYM)

            nb=nBas(iSymb)
            ng=nBas(iSymg)
            IF(nb*ng<=0) Cycle

               ISYMA = ISYMB     !block diagonal Fock Matrix
               ISYMD = ISYMG     !total symmetric densities

             IF (ISYMG.gt.ISYMB .and. iSkip(iSymg).ne.0) THEN

C -------------------------------
C --- F(a,b) = - D(g,d) * (ad|gb)
C -------------------------------
               NK = kOcc(iSymG)

               XdJb(1:NK*NumV*nb) => Wab%A0(1:NK*NumV*nb)

               if(NK.ne.0)then

               ISFSQ = ISTSQ(ISYMB) + ipFSQ(jDen)
               ISMSQ = ISTSQ(ISYMG) + ipMSQ(jDen)


C              Calculate intermediate:
C              X(k,Jb) = Sum(g) C(g,k) * L(g,Jb).
C              ----------------------------------

               CALL DGEMM_('T','N',
     &                    NK,NUMV*NBAS(ISYMB),NBAS(ISYMG),
     &                    ONE,Work(ISMSQ),NBAS(ISYMG),
     &                        LqJs%SB(ISYMG)%A3,NBAS(ISYMG),
     &                    ZERO,XdJb,NK)


C              F(a,b) = F(a,b) - Sum(kJ) X(kJ,a) * X(kJ,b).
C              -------------------------------------------

c               CALL DGEMM_('T','N',
c     &              NBAS(ISYMA),NBAS(ISYMB),NK*NUMV,
c     &              -FactX(jDen),XdJb,NK*NUMV,
c     &                           XdJb,NK*NUMV,
c     &                       ONE,Work(ISFSQ),NBAS(ISYMA))

c *** Compute only the LT part of the exchange matrix ***************
               ipF=0
               LKV=NK*NUMV
               DO jB=1,NBAS(iSymB)
                     NBL = NBAS(iSymA) - (jB-1)
                     ipF = ISFSQ + NBAS(iSymA)*(jB-1) + (jB-1)
                     jjB = 1 + LKV*(jB-1)

                     CALL DGEMV_('T',LKV,NBL,
     &                    -FactX(jDen),XdJb(jjB:),LKV,
     &                                 XdJb(jjB:),1,
     &                             ONE,Work(ipF),1)

               END DO
               XdJb=>Null()
c ******************************************************************

c         write(6,*)'Symmetry block of FSQ= ',isyma
c         write(6,*)'Symmetry block of DSQ= ',isymg
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMA),NBAS(ISYMA))
               endif

C -------------------------------
C --- F(g,d) = - D(a,b) * (ad|gb)
C -------------------------------
               NK = kOcc(iSymB)

               if(NK.ne.0)then

               ISFSQ = ISTSQ(ISYMG) + ipFSQ(jDen)
               ISMSQ = ISTSQ(ISYMB) + ipMSQ(jDen)

               XgJk(1:ng*NumV*NK) => Wab%A0(1:ng*NumV*NK)

C              Calculate intermediate:
C              X(gJ,k) = Sum(a) L(gJ,a) * C(a,k).
C              ----------------------------------

               CALL DGEMM_('N','N',
     &                  NBAS(ISYMG)*NUMV,NK,NBAS(ISYMA),
     &                  ONE,LqJs%SB(ISYMG)%A3,NBAS(ISYMG)*NUMV,
     &                      WORK(ISMSQ),NBAS(ISYMA),
     &                  ZERO,XgJk,NBAS(ISYMG)*NUMV)


C              F(g,d) = F(g,d) - Sum(Jk) X(g,Jk) * X(d,Jk).
C              -------------------------------------------

c               CALL DGEMM_('N','T',
c     &              NBAS(ISYMG),NBAS(ISYMD),NUMV*NK,
c     &              -FactX(jDen),XgJk,NBAS(ISYMG),
c     &                           XgJk,NBAS(ISYMD),
c     &                       ONE,Work(ISFSQ),NBAS(ISYMG))

c *** Compute only the LT part of the exchange matrix ***************
               ipF=0
               LVK=NUMV*NK
               DO jD=1,NBAS(iSymD)
                     NBL = NBAS(iSymG) - (jD-1)
                     ipF = ISFSQ + NBAS(iSymG)*(jD-1) + (jD-1)

                     CALL DGEMV_('N',NBL,LVK,
     &                    -FactX(jDen),XgJk(jD:),NBAS(iSymG),
     &                                 XgJk(jD:),NBAS(iSymD),
     &                             ONE,Work(ipF),1)
               END DO

c         write(6,*)'Symmetry block of FSQ= ',isymg
c         write(6,*)'Symmetry block of DSQ= ',isyma
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMG),NBAS(ISYMG))
               XgJk=>Null()
               endif

             ENDIF
             LqJs%SB(ISYMG)%A3=>Null()

            END DO  ! loop over orbital symmetries

            CALL CWTIME(TC2X2,TW2X2)
            texch(1) = texch(1) + (TC2X2 - TC2X1)
            texch(2) = texch(2) + (TW2X2 - TW2X1)


C ************ END "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ***********


       END IF  ! DoExchange section

      END DO  ! end of the loop over densities

      ENDIF ! jSym.ne.1


C --- Free the memory
         Call Deallocate_SBA(Wab)

      END DO  ! end of the batch procedure

C -- Close Files
      Do ksym=1,nSym
       if (Lunit(ksym).ne.-1) then
         Call DACLOS(Lunit(ksym))
         Lunit(ksym) = -1
       end if
      End do

      END DO  ! loop over JSYM

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1


*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky SCF timing from '//SECNAM
      Write(6,CFmt)'------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ/REORDER VECTORS             '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'EXCHANGE                         '
     &                           //'         ',texch(1),texch(2)

         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

c Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in SCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM CHO_FTWO_MO.'
      WRITE(6,'(6X,A)')'***** EXCHANGE MATRIX ***** '
      Do jDen=1,nDen
          WRITE(6,'(6X,A,I2)')'DENSITY TYPE: ',jDen
          WRITE(6,'(6X,A,I2)')'DoExchange: ',DoExchange(jDen)
          WRITE(6,*)
      if(DoExchange(jDen))then
      icount=0
      DO ISYM=1,NSYM
      ISQ=ipFSQ(jDen)+icount
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call cho_output(Work(ISQ),1,NB,1,NB,NB,NB,1,6)
          icount=icount+NB**2
        END IF
      END DO
      endif
      End Do

      endif

#endif

      rc  = 0

      Return
      END

**************************************************************
