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
      SUBROUTINE CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,
     &           FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)

************************************************************************
*  Author : F. Aquilante
*
*  Purpose:
*          For a given set of total symmetric AO-densities
*          (hermitian matrix representation) this routine
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
*          Exchange term (Drs is the squared matrix
*                        of r=s symmetry block.)
*
*    F(pq) = F(pq) - FactX * Sum_rJ  Lpr,J * (Sum_s Drs*Lqs,J)
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
*  ip{X}FS(nDen) : pointer to the array containing {X} in SQ storage
*    {X=D,F --- Density, Fock matrix}
*
*  ipNocc(nDen) : pointer to the array of the Occupation numbers
*                 for the corresponding density
*
*  MinMem(nSym) : minimum amount of memory required to read
*                 a single Cholesky vector in full storage
*
************************************************************************

      Implicit Real*8 (a-h,o-z)

      Integer  rc,nDen,nBas(*)
      Real*8   FactC(nDen),FactX(nDen)
      Integer  KSQ1(8),ISTSQ(8),ISTLT(8),iSkip(8),MinMem(*)
      Integer  ipDLT(nDen),ipDSQ(nDen),ipNocc(nDen)
      Integer  ipFLT(nDen),ipFSQ(nDen)
      Logical  DoExchange(nDen),DoCoulomb(nDen),DoSomeX,DoSomeC,Debug
      Real*8   tread(2),tcoul(2),texch(2)
      Logical  timings

      COMMON   /CHOTIME /timings
      Character*50 CFmt
      Character*15 SECNAM
      Parameter (SECNAM = 'CHO_FOCKTWO_RED')

      parameter (zero = 0.0D0, one = 1.0D0)
      Logical DoRead
      parameter (DoRead = .true.)


#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
**************************************************


#ifdef _DEBUG_
c      Debug=.true.
      Debug=.false.! to avoid double printing in SCF-debug
#else
      Debug=.false.
#endif
      IREDC = -1  ! unknown reduced set in core

        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time read/rreorder vectors
           tcoul(i) = zero  !time for computing Coulomb
           texch(i) = zero  !time for computing Exchange
        end do

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

C *************** BIG LOOP OVER VECTORS SYMMETRY *****************
      DO jSym=1,MaxSym

      If (NumCho(jSym).lt.1) GOTO 1000


C ------------------------------------------------------

C --- Pointers to be used for the vectors
         do i=1,nSym
            KSQ1(i)=-6666
         end do

C SET UP THE READING
C ------------------
      Call GetMem('Maxmem','MAX ','REAL',KDUM,LWORK)

      If (MinMem(jSym) .gt. 0) Then
         nVec = Min(LWORK/MinMem(jSym),NumCho(jSym))
      Else
C         ***QUIT*** bad initialization
         WRITE(6,*) 'Cho_FockTwo_RED: bad initialization'
         rc=99
         CALL QTrace()
         CALL Abend()
         nVec = -9999  ! dummy assignment - avoid compiler warnings
      End If

C...this is ONLY for debug:
c      nVec = min(nVec,1)

      If (nVec .gt. 0) Then
         nBatch = (NumCho(jSym) - 1)/nVec + 1
      Else
C         ***QUIT*** insufficient memory
         WRITE(6,*) 'Cho_FockTwo_RED: Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'min. mem. need= ',MinMem(jSym)
         WRITE(6,*) 'NumCho= ',NumCho(jsym)
         WRITE(6,*) 'jsym= ',jsym
         rc = 205
         CALL QTrace()
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
         kTOT = MinMem(jSym)*NumV
         Call GetMem('MemR','ALLO','REAL',kWab,kTOT)

c--- setup the skipping flags according to # Occupied
      do k=1,nSym
         iSkip(k)=0
         l=Muld2h(k,jsym) !L(kl) returned if nOcc(k or l) .ne.0
         if (k.eq.l) then
            iSkip(k) = 666 ! always contribute to Coulomb
         else
          do jDen=1,nDen
               iSkip(k) = iSkip(k) + iWork(ipNocc(jDen)+k-1)
     &                             + iWork(ipNocc(jDen)+l-1)
          end do
         endif
      end do


C Max dimension of a read symmetry block
        Nmax=0
        Do iSym=1,nSym
           if(NBAS(iSym).gt.Nmax .and. iSkip(iSym).ne.0)then
           Nmax = NBAS(iSym)
           iSymMax= iSym
           endif
        End Do

C --- Compute the total amount of memory needed to have the
C --- vectors in core (full storage)
         lChoV=0
         Do iSymq=1,nSym
            iSymp = muld2h(jSym,iSymq)
            IF (nBas(iSymq)*nBas(iSymp).ne.0) THEN
             If(iSymp.gt.iSymq .and. iSkip(iSymp).ne.0) then
               KSQ1(iSymp) = kWab + lChoV
               lChoV = lChoV + Nbas(iSymp)*Nbas(iSymq)*NumV
             Else
               if(iSymp.eq.iSymq .and. iSkip(iSymp).ne.0) then
                 KSQ1(iSymp) = kWab + lChoV
C --- Special trick for the vector L11 ; used to store X(a,Jb)
                 if(iSymp.eq.1.and.jSym.eq.1.and.DoSomeX)then
                    lChoV = lChoV + (Nmax**2) * NumV
                 else
                    lChoV = lChoV + Nbas(iSymp)*(Nbas(iSymp)+1)/2*NumV
                 endif
               endif
             Endif
            ENDIF
         End Do

         kLab = kWab + lChoV
         lScr = kTOT - lChoV

C --- Reading of the vectors is done in Reduced sets

      iSwap = min(1,(jSym-1)) ! L(ab,J) --> L(a,J,b) iff jSym.ne.1

      CALL CWTIME(TCR1,TWR1)

      Call CHO_X_getVfull(irc,Work(kLab),lscr,iVEC,NumV,jSym,iSwap,
     &                    IREDC,KSQ1,iSkip,DoRead)

       if (irc.ne.0) then
          rc = irc
          RETURN
       endif

       CALL CWTIME(TCR2,TWR2)
       tread(1) = tread(1) + (TCR2 - TCR1)
       tread(2) = tread(2) + (TWR2 - TWR1)

#ifdef _DEBUG_
       write(6,*) 'Batch ',iBatch,' of   ',nBatch,': NumV = ',NumV
       write(6,*) 'Total allocated :     ',kTOT,' at ',kWab
       write(6,*) 'Memory pointers KSQ1: ',(KSQ1(i),i=1,nSym)
       write(6,*) 'kLab:                 ',kLab
       write(6,*) 'lScr:                 ',lScr
       write(6,*) 'JSYM:                 ',jSym
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

         CALL FZERO(Work(kLab),NumV)
         DO iSymr=1,nSym
         IF(nBas(iSymr).ne.0.and.iWork(ipNocc(jDen)+iSymr-1).ne.0)THEN

         ISDLT = ipDLT(jDen) + ISTLT(ISYMR)
         Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

         CALL DGEMV_('T',Naa,NumV,
     &              ONE,Work(KSQ1(iSymr)),Naa,
     &              Work(ISDLT),1,ONE,Work(kLab),1)

         ENDIF
         End DO

C --- Contraction with the 2nd vector
C ---------------------------------------
C --- Fpp{#J} <- Fpp{#J} + FactC * sum_J L(pp,{#J})*V{#J}
C==========================================================
C
          DO iSyms=1,nSym
          IF(nBas(iSyms).ne.0)THEN

           ISFLT = ipFLT(jDen) + ISTLT(ISYMS)
           Naa = nBas(iSyms)*(nBas(iSyms)+1)/2

           CALL DGEMV_('N',Naa,NumV,
     &              FactC(jDen),Work(KSQ1(iSyms)),Naa,
     &              Work(kLab),1,ONE,Work(ISFLT),1)

          ENDIF
          End DO

        CALL CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1) + (TCC2 - TCC1)
        tcoul(2) = tcoul(2) + (TWC2 - TWC1)

       ENDIF  ! jSym=1 & DoCoulomb

      END DO ! loop over densities for the Coulomb term

C *************** END COULOMB CONTRIBUTION  ****************

       IF (jSym.eq.1.and.DoSomeX) THEN
C WE INCLUDE THE DIAGONAL PART OF THE EXCHANGE HERE!!!!
C
C ********** REORDER AND SQUARE ONE VECTOR AT THE TIME *****
C     Reorder: L(rs,J) -> L(r,J,s).
C     CHOVEC(nrs,numv) ---> CHOVEC(nr,numv,ns)

         KQS1=0
         KQS2=0
         ISFSQ=0
         ISDSQ=0
         DO iSymr=1,nSym

          iSymr_Occ=0
          do jDen=1,nDen
             iSymr_Occ = iSymr_Occ + iWork(ipNocc(jDen)+iSymr-1)
          end do
          IF(nBas(iSymr).ne.0.and.iSymr_Occ.gt.0)THEN

             CALL CWTIME(TCREO1,TWREO1)

             KQS1 = KSQ1(ISYMR)
             KQS2 = kLab
             NRS  = NBAS(iSymr)*(NBAS(iSymr)+1)/2

             KOFF1=0
             KOFF2=0
             KOFF3=0
             iSQ = NBAS(iSymr)*NUMV
             Do JVEC=1,NUMV

                iLT = NRS*(JVEC-1)
                iVR = NBAS(iSymr)*(JVEC-1)

                Do jR=1,NBAS(iSymr)
                   Do jS=jR,NBAS(iSymr)

                      jSR = iTri(jS,jR)
                      KOFF1 = (KQS1-1) + iLT + jSR

                      KOFF2 = (KQS2-1) + iSQ*(jS-1) + iVR + jR

                      WORK(KOFF2) = WORK(KOFF1)

                      KOFF3 = (KQS2-1) + iSQ*(jR-1) + iVR + jS

                      WORK(KOFF3) = WORK(KOFF1)

                   End Do
                End Do
             End Do

             CALL CWTIME(TCREO2,TWREO2)
             tread(1) = tread(1) + (TCREO2 - TCREO1)
             tread(2) = tread(2) + (TWREO2 - TWREO1)

          DO jDen=1,nDen

           IF (DoExchange(jDen)) THEN

               IF (iWork(ipNocc(jDen)+iSymr-1).ne.0) THEN

               CALL CWTIME(TC1X1,TW1X1)

C              Calculate intermediate:
C              X(p,Js) = Sum(q) D(p,q) * L(q,Js).
C              -----------------------------------
               ISDSQ = ISTSQ(ISYMR) + ipDSQ(jDen)
               ISFSQ = ISTSQ(ISYMR) + ipFSQ(jDen)


               CALL DGEMM_('N','N',
     &               NBAS(ISYMR),NUMV*NBAS(ISYMR),NBAS(ISYMR),
     &               ONE,WORK(ISDSQ),NBAS(ISYMR),
     &                   Work(KQS2), NBAS(ISYMR),
     &               ZERO,WORK(kWab),NBAS(ISYMR))

C              Calculate exchange contribution:
C              F(p,q) = F(p,q) - FactX Sum(sJ) L(sJ,p) * X(sJ,q).
C              --------------------------------------------------

               CALL DGEMM_('T','N',
     &                NBAS(ISYMR),NBAS(ISYMR),NBAS(ISYMR)*NUMV,
     &             -FactX(jDen),WORK(KQS2),NBAS(ISYMR)*NUMV,
     &                          WORK(kWab),NBAS(ISYMR)*NUMV,
     &                      One,Work(ISFSQ),NBAS(ISYMR))


               CALL CWTIME(TC1X2,TW1X2)
               texch(1) = texch(1) + (TC1X2 - TC1X1)
               texch(2) = texch(2) + (TW1X2 - TW1X1)

               ENDIF

          ENDIF  ! DoExchange(jDen)

         END DO  ! loop over the densities

c         write(6,*)'Symmetry block of FSQ= ',isymr
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMR),NBAS(ISYMR))

        ENDIF  ! nbas.ne.0 & nOcc.ne.0

        END DO  !loop over Fock mat symmetries

       END IF  ! jSym=1 and DoSomeX

C ****************** END DIAGONAL EXCHANGE **********************
C
C ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************
C
       IF (jSym.ne.1) THEN

       DO jDen=1,nDen

       IF (DoExchange(jDen)) THEN
C --- COMPUTE EXCHANGE FOR OFF-DIAGONAL VECTORS

            CALL CWTIME(TC2X1,TW2X1)

            DO ISYMG = 1,NSYM

               ISYMB = MULD2H(ISYMG,JSYM)

            IF(nBas(iSymb)*nBas(iSymg).ne.0)THEN

               ISYMD = ISYMG     !only total symmetric Density
               ISYMA = ISYMB     !block diagonal Fock Matrix

             IF (ISYMG.gt.ISYMB .and. iSkip(iSymg).ne.0) THEN

               KQS1 = KSQ1(ISYMG)

C -------------------------------
C --- F(a,b) = - D(g,d) * (ad|gb)
C -------------------------------
               if (iWork(ipNocc(jDen)+iSymg-1).ne.0) then

               ISFSQ = ISTSQ(ISYMB) + ipFSQ(jDen)
               ISDSQ = ISTSQ(ISYMG) + ipDSQ(jDen)

C              Calculate intermediate:
C              X(d,Jb) = Sum(g) D(d,g) * L(g,Jb).
C              ----------------------------------

               CALL DGEMM_('N','N',
     &                    NBAS(ISYMD),NUMV*NBAS(ISYMB),NBAS(ISYMG),
     &                    ONE,Work(ISDSQ),NBAS(ISYMD),Work(KQS1),
     &                    NBAS(ISYMG),ZERO,WORK(kLab),NBAS(ISYMD))


C              F(a,b) = F(a,b) - Sum(dJ) L(dJ,a) * X(dJ,b).
C              -------------------------------------------

              CALL DGEMM_('T','N',
     &              NBAS(ISYMA),NBAS(ISYMB),NBAS(ISYMD)*NUMV,
     &              -FactX(jDen),WORK(KQS1),NBAS(ISYMD)*NUMV,
     & WORK(kLab),NBAS(ISYMD)*NUMV,ONE,Work(ISFSQ),NBAS(ISYMA))

c         write(6,*)'Symmetry block of FSQ= ',isyma
c         write(6,*)'Symmetry of the density= ',isymg
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMA),NBAS(ISYMA))
               endif
C -------------------------------
C --- F(g,d) = - D(a,b) * (ad|gb)
C -------------------------------
               if (iWork(ipNocc(jDen)+iSyma-1).ne.0) then

               ISFSQ = ISTSQ(ISYMG) + ipFSQ(jDen)
               ISDSQ = ISTSQ(ISYMB) + ipDSQ(jDen)

C              Calculate intermediate:
C              X(gJ,b) = Sum(a) L(gJ,a)* D(a,b).
C              ----------------------------------

               CALL DGEMM_('N','N',
     &                  NBAS(ISYMG)*NUMV,NBAS(ISYMB),NBAS(ISYMA),
     &                  ONE,Work(KQS1),NBAS(ISYMG)*NUMV,WORK(ISDSQ),
     &                  NBAS(ISYMA),ZERO,WORK(kLab),NBAS(ISYMG)*NUMV)


C              F(g,d) = F(g,d) - Sum(Jb) X(g,Jb) * L(d,Jb).
C              -------------------------------------------

              CALL DGEMM_('N','T',
     &              NBAS(ISYMG),NBAS(ISYMD),NUMV*NBAS(ISYMB),
     &              -FactX(jDen),WORK(kLab),NBAS(ISYMG),
     & WORK(KQS1),NBAS(ISYMD),ONE,Work(ISFSQ),NBAS(ISYMG))

c         write(6,*)'Symmetry block of FSQ= ',isymg
c         write(6,*)'Symmetry of the density= ',isyma
c         call recprt('FSQ','',Work(ISFSQ),NBAS(ISYMG),NBAS(ISYMG))
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
         Call GetMem('MemR','FREE','REAL',kWab,kTOT)

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
      Write(6,CFmt)'-----------------------------------------'
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
#ifdef _DEBUG_

      if(Debug) then !to avoid double printing in SCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM CHO_FOCKTWO_RED.'
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
