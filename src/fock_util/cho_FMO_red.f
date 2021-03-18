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
      SUBROUTINE CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                       lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                       MinMem,ipMSQ,ipNocc)

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
*          The vectors are indeed read in reduced set storage
*          and then reordered.
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
*  ipNocc(nDen) : pointer to the array of the Occupation numbers
*                 for the corresponding density
*
*  MinMem(nSym) : minimum amount of memory required to read
*                 a single Cholesky vector in full storage
*
************************************************************************
      use Data_structures, only: SBA_Type
      Implicit Real*8 (a-h,o-z)

      Integer   rc,nDen,kOcc(8),KSQ1(8)
      Real*8    FactC(nDen),FactX(nDen)
      Integer   ISTSQ(8),ISTLT(8),iSkip(8),lOff1
      Real*8    tread(2),tcoul(2),texch(2)
      Integer   ipDLT(nDen),ipDSQ(nDen),ipFLT(nDen),ipFSQ(nDen)
      Integer   ipMSQ(nDen),ipNocc(nDen),MinMem(*)
      Logical DoExchange(nDen),DoCoulomb(nDen),DoSomeX,DoSomeC
#ifdef _DEBUGPRINT_
      Logical Debug
#endif
      Logical Square
      Character*50 CFmt
      Character(LEN=11), Parameter :: SECNAM = 'CHO_FMO_RED'
#include "chodensity.fh"
#include "chotime.fh"
#include "choscf.fh"
#include "real.fh"

      Real*8, parameter :: xone = -One
      Logical, Parameter :: DoRead = .true.

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Real*8, Allocatable:: DChk(:)

      Type (SBA_Type), Target :: Wab
      Real*8, Pointer :: LrJs(:,:,:)=>Null(), VJ(:)=>Null(),
     &                   Scr(:)=>Null(), XkJb(:,:,:)=>Null(),
     &                   XgJk(:,:,:)=>Null(), XkJs(:,:,:)=>Null()

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      nOcc(jSym,jDen) = iWork(ipNocc(jDen)-1+jSym)
**************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in SCF-debug
      DensityCheck=.true.
#endif
      IREDC = -1  ! unknown reduced set in core

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
        xf = 2.0d0
        if (DECO) xf=1.0d0
        Do jDen=1,nDen
          Do jSym=1,nSym
           if (nBas(jSym).ne.0.and.nOcc(jSym,jDen).ne.0) then
             Call mma_allocate(Dchk,nBas(jSym)**2,Label='Dchk')
             Call Cho_X_Test(Work(ipDSQ(jDen)+ISTSQ(jSym)),nBas(jSym),
     &                       Square,Work(ipMSQ(jDen)+ISTSQ(jSym)),
     &                       nOcc(jSym,jDen),xf,Dchk,
     &                       nBas(jSym)**2,Thr,irc)
             if(irc.eq.0)then
                write(6,*)'*** DENSITY CHECK : OK! *** SYMM= ',jSym
             else
                write(6,*)'*** DENSITY CHECK FAILED IN SYMMETRY: ',
     &                    jSym,'FOR THE DENSITY NR. ',jDen
                xnormY = sqrt(ddot_(nBas(jSym)**2,Dchk,1,Dchk,1))
                write(6,*)'    Norm of the residual vector = ',xnormY
             endif
             Call mma_deallocate(Dchk)
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

C *************** BIG LOOP OVER VECTORS SYMMETRY *****************
      DO jSym=1,MaxSym

      If (NumCho(jSym).lt.1) GOTO 1000

C Total length of the vectors
      do i=1,nSym
         kSQ1(i)   = -6666
      end do
C ------------------------------------------------------


C SET UP THE READING
C ------------------
      Call mma_maxDBLE(LWORK)

      If (MinMem(jSym) .gt. 0) Then
         nVec = Min(LWORK/MinMem(jSym),NumCho(jSym))
      Else
C         ***QUIT*** bad initialization
         WRITE(6,*) 'Cho_FMO_red: bad initialization'
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
         WRITE(6,*) 'Cho_FMO_red: Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'min. mem. need= ',MinMem(jSym)
         WRITE(6,*) 'NumCho= ',NumCho(jsym)
         WRITE(6,*) 'jsym= ',jsym
         rc = 33
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

         kTOT = MinMem(jSym)*NumV
         Call mma_allocate(Wab%A0,kTOT,Label='Wab')
         kWab = ip_of_Work(Wab%A0(1))

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

C --- Compute the total amount of memory needed to have the
C --- vectors in core (full storage)
         iE=0
         iadd = 0
         Do iSymq=1,nSym
            iSymp = muld2h(jSym,iSymq)
            nq=nBas(iSymq)
            np=nBas(iSymp)
            IF (nq*np<=0) Cycle
            iS = iE + 1

            If(iSymp.gt.iSymq .and. iSkip(iSymp).ne.0) then
               KSQ1(iSymp) = kWab + iE
               iE = iE + np*nq*NumV

               Wab%SB(iSymq)%A2(1:np*nq,1:NumV) => Wab%A0(iS:iE)
            Else
               If(iSymp.eq.iSymq .and. iSkip(iSymp).ne.0) then
                 KSQ1(iSymp) = kWab + iE
C --- Special trick for the vector L11 ; used to store X(a,Jb)
                 npp=np*(np+1)/2
                 if(iSymp.eq.1.and.jSym.eq.1.and.DoSomeX)then
                    iadd = Max(0,lOff1 - npp)
                 endif
                 iE = iE + np*(np+1)/2*NumV
                 Wab%SB(iSymq)%A2(1:npp,1:NumV) => Wab%A0(iS:iE)
               Endif
            Endif
         End Do
         iE = iE + iadd*NumV

         lChoV = iE
         lScr = kTOT - lChoV
         iS = iE + 1
         Scr(1:lScr) => Wab%A0(iS:kTOT)


C --- Reading of the vectors is done in Reduced sets
      iSwap = min(1,(jSym-1)) ! L(ab,J) --> L(a,J,b) iff jSym.ne.1

      CALL CWTIME(TCR1,TWR1)

      Call CHO_X_getVfull(irc,Scr,lScr,iVEC,NumV,jSym,iSwap,
     &                    IREDC,KSQ1,iSkip,DoRead)

      CALL CWTIME(TCR2,TWR2)
      tread(1) = tread(1) + (TCR2 - TCR1)
      tread(2) = tread(2) + (TWR2 - TWR1)

         Scr => Null()

#ifdef _DEBUGPRINT_
       write(6,*) 'Batch ',iBatch,' of   ',nBatch,': NumV = ',NumV
       write(6,*) 'Total allocated :     ',kTOT
       write(6,*) 'Memory pointers KSQ1: ',(KSQ1(i),i=1,nSym)
       write(6,*) 'lScr:                 ',lScr
       write(6,*) 'lOff1:                ',lOff1
       write(6,*) 'JSYM:                 ',jSym
#endif
       if (irc.ne.0) then
          rc = irc
          RETURN
       endif

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

         VJ(1:NumV) => Wab%A0(iS:iS-1+NumV)
         VJ(:)=Zero

         DO iSymr=1,nSym
            IF(nBas(iSymr).ne.0.and.nOcc(isymr,jDen).ne.0)THEN

            ISDLT = ipDLT(jDen) + ISTLT(ISYMR)
            Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

            CALL DGEMV_('T',Naa,NumV,
     &                 ONE,Wab%SB(iSymr)%A2,Naa,
     &                 Work(ISDLT),1,ONE,VJ,1)

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
     &                 FactC(jDen),Wab%SB(iSyms)%A2,Naa,
     &                 VJ,1,ONE,Work(ISFLT),1)

c              WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMS
c              CALL TRIPRT('Coulomb FLT',' ',Work(ISFLT),nBas(iSyms))

             ENDIF
          End DO

         CALL CWTIME(TCC2,TWC2)
         tcoul(1) = tcoul(1) + (TCC2 - TCC1)
         tcoul(2) = tcoul(2) + (TWC2 - TWC1)

         VJ=>Null()

       ENDIF  ! jSym=1 & DoCoulomb

      END DO ! loop over densities for the Coulomb term

C *************** END COULOMB CONTRIBUTION  ****************

C **********************************************************
C WE INCLUDE THE DIAGONAL PART OF THE EXCHANGE HERE!!!!
C
C ********** REORDER and SQUARE ONE VECTOR AT THE TIME *****
C
C     CHOVEC(nrs,numv) ---> CHOVEC(nr,numv,ns)

      IF (jSym.eq.1 .and. DoSomeX) THEN

        ISFSQ=0
        ISMSQ=0
        DO iSymr=1,nSym

          iSymr_Occ=0
          do jDen=1,nDen
             iSymr_Occ = iSymr_Occ + nOcc(iSymr,jDen)
          end do
          IF(nBas(iSymr).ne.0.and.iSymr_Occ.ne.0)THEN

             CALL CWTIME(TCREO1,TWREO1)

             nr=NBAS(iSymr)
             LrJs(1:nr,1:NUMV,1:nr) => Wab%A0(iS:iS-1+nr*NUMV*nr)

             Do JVEC=1,NUMV

                Do jR=1,NBAS(iSymr)
                   Do jS=jR,NBAS(iSymr)

                      jSR = iTri(jS,jR)

                      LrJs(jR,JVEC,jS) = Wab%SB(iSymr)%A2(jSR,JVEC)


                      LrJs(jS,JVEC,jR) = Wab%SB(iSymr)%A2(jSR,JVEC)

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

               XkJs(1:NK,1:NUMV,1:NBAS(ISYMR)) =>
     &            Wab%A0(1:NK*NUMV*NBAS(ISYMR))

               CALL DGEMM_('T','N',
     &               NK,NUMV*NBAS(ISYMR),NBAS(iSYMR),
     &               ONE,WORK(ISMSQ),NBAS(ISYMR),
     &                   LrJs, NBAS(ISYMR),
     &               ZERO,XkJs,NK)

C              F(s,t) = F(s,t) - Sum(kJ) X(kJ,s)*X(kJ,t).
C              ------------------------------------------
c               CALL DGEMM_('T','N',
c     &               NBAS(ISYMR),NBAS(iSYMR),NK*NUMV,
c     &               -FactX(jDen),XkJs,NK*NUMV,
c     &                   XkJs,NK*NUMV,
c     &                   One,Work(ISFSQ),NBAS(ISYMR))

c *** Compute only the LT part of the exchange matrix ***************
C
               ipF=0
               LKV=NK*NUMV
               DO jS=1,NBAS(iSymS)
                     NBL = NBAS(iSymR) - (jS-1)
                     ipF = ISFSQ + NBAS(iSymR)*(jS-1) + (jS-1)

                     CALL DGEMV_('T',LKV,NBL,
     &                    -FactX(jDen),XkJs(1,1,jS),LKV,
     &                                 XkJs(1,1,jS),1,
     &                        ONE,Work(ipF),1)

               END DO


             CALL CWTIME(TC1X2,TW1X2)
             texch(1) = texch(1) + (TC1X2 - TC1X1)
             texch(2) = texch(2) + (TW1X2 - TW1X1)

c         write(6,*)'Symmetry block of FSQ= ',isymr
c         CALL RECPRT('FSQ',' ',Work(ISFSQ),NBAS(iSYMR),NBAS(iSYMR))
             XkJs=>Null()

           ENDIF   ! if kocc.ne.0

           END IF  ! DoExchange(jDen)

          END DO  ! end of the loop over densities
          LrJs=> Null()

         ENDIF  ! if nbas.ne.0 & iSymr_nOcc.ne.0

        END DO  !loop over Fock mat symmetries

      ENDIF  ! jSym=1 & DoSomeX

C ****************** END DIAGONAL EXCHANGE **********************
C
C ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************
C
       IF (jSym.ne.1) THEN

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

            IF(nBas(iSymb)*nBas(iSymg).ne.0)THEN

               ISYMA = ISYMB     !block diagonal Fock Matrix
               ISYMD = ISYMG     !total symmetric densities

             IF (ISYMG.gt.ISYMB .and. iSkip(iSymg).ne.0) THEN

C -------------------------------
C --- F(a,b) = - D(g,d) * (ad|gb)
C -------------------------------
               NK = kOcc(iSymG)

               if(NK.ne.0)then

               ISFSQ = ISTSQ(ISYMB) + ipFSQ(jDen)
               ISMSQ = ISTSQ(ISYMG) + ipMSQ(jDen)

               XkJb(1:NK,1:NUMV,1:NBAS(ISYMB)) =>
     &             Wab%A0(iS:iS-1+NK*NUMV*NBAS(ISYMB))

C              Calculate intermediate:
C              X(k,Jb) = Sum(g) C(g,k) * L(g,Jb).
C              ----------------------------------

               CALL DGEMM_('T','N',
     &                    NK,NUMV*NBAS(ISYMB),NBAS(ISYMG),
     &                    ONE,Work(ISMSQ),NBAS(ISYMG),
     &                        Wab%SB(ISYMB)%A2,NBAS(ISYMG),
     &                    ZERO,XkJb,NK)


C              F(a,b) = F(a,b) - Sum(kJ) X(kJ,a) * X(kJ,b).
C              -------------------------------------------

c               CALL DGEMM_('T','N',
c     &               NBAS(ISYMA),NBAS(iSYMB),NK*NUMV,
c     &               -FactX(jDen),XkJb,NK*NUMV,
c     &                     XkJb,NK*NUMV,
c     &                One,Work(ISFSQ),NBAS(ISYMA))

c *** Compute only the LT part of the exchange matrix ***************
               ipF=0
               LKV=NK*NUMV
               DO jB=1,NBAS(iSymB)
                     NBL = NBAS(iSymA) - (jB-1)
                     ipF = ISFSQ + NBAS(iSymA)*(jB-1) + (jB-1)

                     CALL DGEMV_('T',LKV,NBL,
     &                    -FactX(jDen),XkJb(1,1,jB),LKV,
     &                                 XkJb(1,1,jB),1,
     &                      ONE,Work(ipF),1)

               END DO
c ******************************************************************

c         write(6,*)'Symmetry block of FSQ= ',isyma
c         write(6,*)'Symmetry block of DSQ= ',isymg
c         CALL RECPRT('FSQ',' ',Work(ISFSQ),NBAS(iSYMA),NBAS(iSYMA))

               XkJb=>Null()

               endif

C -------------------------------
C --- F(g,d) = - D(a,b) * (ad|gb)
C -------------------------------
               NK = kOcc(iSymB)

               if(NK.ne.0)then

               ISFSQ = ISTSQ(ISYMG) + ipFSQ(jDen)
               ISMSQ = ISTSQ(ISYMB) + ipMSQ(jDen)

               XgJk(1:NBAS(ISYMG),1:NUMV,1:NK) =>
     &            Wab%A0(iS:iS-1+NBAS(ISYMG)*NUMV*NK)

C              Calculate intermediate:
C              X(gJ,k) = Sum(a) L(gJ,a) * C(a,k).
C              ----------------------------------

               CALL DGEMM_('N','N',
     &                  NBAS(ISYMG)*NUMV,NK,NBAS(ISYMA),
     &                  ONE,Wab%SB(ISYMB)%A2,NBAS(ISYMG)*NUMV,
     &                      WORK(ISMSQ),NBAS(ISYMA),
     &                  ZERO,XgJk,NBAS(ISYMG)*NUMV)


C              F(g,d) = F(g,d) - Sum(Jk) X(g,Jk) * X(d,Jk).
C              -------------------------------------------

c               CALL DGEMM_('N','T',
c     &               NBAS(ISYMG),NBAS(ISYMD),NK*NUMV,
c     &               -FactX(jDen),XgJk,NBAS(ISYMG),
c     &                   XgJk,NBAS(ISYMD),
c     &                   One,Work(ISFSQ),NBAS(ISYMG))

c *** Compute only the LT part of the exchange matrix ***************
               ipF=0
               LVK=NUMV*NK
               DO jD=1,NBAS(iSymD)
                     NBL = NBAS(iSymG) - (jD-1)
                     ipF = ISFSQ + NBAS(iSymG)*(jD-1) + (jD-1)

                     CALL DGEMV_('N',NBL,LVK,
     &                    -FactX(jDen),XgJk(jD,1,1),NBAS(iSymG),
     &                           XgJk(jD,1,1),NBAS(iSymD),
     &                           ONE,Work(ipF),1)

               END DO
c ******************************************************************

c         write(6,*)'Symmetry block of FSQ= ',isymg
c         write(6,*)'Symmetry block of DSQ= ',isyma
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMG),NBAS(ISYMG))

               XgJk=>Null()

               endif

             ENDIF

            ENDIF
            END DO  ! loop over orbital symmetries

            CALL CWTIME(TC2X2,TW2X2)
            texch(1) = texch(1) + (TC2X2 - TC2X1)
            texch(2) = texch(2) + (TW2X2 - TW2X1)

C ************ END "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ***********


         END IF  ! DoExchange section

         END DO  ! end of the loop over densities

       ENDIF ! jSym.ne.1

C --- Free the memory
         Call mma_deallocate(Wab%A0)

      END DO  ! end of the batch procedure

1000  CONTINUE

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

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
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
      Endif
      End Do

      endif

#endif

      rc  = 0

      Return
      END

**************************************************************
