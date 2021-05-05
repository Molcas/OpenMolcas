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
      SUBROUTINE CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,
     &           FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)

************************************************************************
*  Author : F. Aquilante
*
*  Purpose:
*          For a given set of total symmetric AO-densities
*          (hermitian matrix representation) this routine
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
*  pNocc(nDen) : pointer to the array of the Occupation numbers
*                 for the corresponding density
*
*  MinMem(nSym) : minimum amount of memory required to read
*                 a single Cholesky vector in full storage
*
************************************************************************
      use Data_Structures, only: SBA_Type, Deallocate_SBA
      use Data_Structures, only: DSBA_Type
      Implicit Real*8 (a-h,o-z)

      Integer   rc,nDen,nSym,nBas(nSym),NumCho(nSym),iSkip(nSym)
      Real*8    FactC(nDen),FactX(nDen)
      Integer   MinMem(nSym)

      Type (DSBA_Type) DLT(nDen), FLT(nDen), FSQ(nDen), DSQ(nDen)

      Type Integer_Pointer
          Integer, Pointer :: I1(:)=>Null()
      End Type Integer_Pointer
      Type (Integer_Pointer) :: pNocc(nDen)


#include "chounit.fh"
#include "real.fh"
#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   DoExchange(nDen),DoCoulomb(nDen),DoSomeX,DoSomeC
      Real*8    tread(2),tcoul(2),texch(2)

#include "chotime.fh"
      Character*50 CFmt
      Character(LEN=11), Parameter :: SECNAM = 'CHO_FOCKTWO'

      Character*6 Fname
      Character(LEN=4), Parameter :: BaseNm = 'CHFV'

#include "stdalloc.fh"

      Type (SBA_Type), Target:: Wab, LqJs
      Real*8, Pointer :: VJ(:)=>Null(), LrJs(:,:,:)=>Null(),
     &                   XpJs(:)=>Null(), XdJb(:)=>Null()

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      nOcc(k,jDen)= pNocc(jDen)%I1(k)
**************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in SCF-debug
#endif

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      ! 1 --> CPU   2 --> Wall
      tread(:) = zero  !time read/rreorder vectors
      tcoul(:) = zero  !time for computing Coulomb
      texch(:) = zero  !time for computing Exchange

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

      Call get_iarray('NumCho',NumCho,nSym)


C *************** BIG LOOP OVER VECTORS SYMMETRY *****************
      DO jSym=1,MaxSym

      If (NumCho(jSym).lt.1) Cycle

C Open Files
       LUnit(:) = -1
       Do ksym=1,nSym
        if (nBas(ksym).ne.0) then
          iSymp=MulD2h(ksym,jSym)
          If (iSymp.ge.ksym) then
             lu = 7
             lu = isfreeunit(lu)
             Lunit(iSymp) = lu
             Write(Fname,'(A4,I1,I1)') BaseNm,iSymp,ksym
             Call DANAME_MF_WA(Lunit(iSymp),Fname)
          End If
        End If
       End Do

C SET UP THE READING
C ------------------
      Call mma_maxDBLE(LWORK)

      If (MinMem(jSym) .gt. 0) Then
         nVec = Min(LWORK/MinMem(jSym),NumCho(jSym))
      Else
C         ***QUIT*** bad initialization
         WRITE(6,*) 'Cho_FockTwo: bad initialization'
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
         WRITE(6,*) 'Cho_FockTwo: Insufficient memory for batch'
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
        Write (6,*) 'iBatch=',iBatch

         If (iBatch .eq. nBatch) Then
            NumV = NumCho(jSym) - nVec*(nBatch - 1)
         Else
            NumV = nVec
         End If

         iVec = nVec*(iBatch - 1) + 1

C --- Allocate memory for reading the vectors
C     We allocate the array as an element of a SBA_type.
         kRdMem = MinMem(jSym)*NumV
         Call mma_allocate(Wab%A0,kRdMem,Label='Wab%A0')
         Wab%nsym=nSym
         Wab%isym=jSym
         Wab%icase=7

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

         iS = iE + 1

         nk=nBas(kSym)
         np=nBas(iSymp)
         IF (nk*np<=0) Cycle

         If (iSymp.eq.ksym .and. iSkip(iSymp).ne.0) then
           NumB = nk*(nk+1)/2
C --- Special trick for the vector L11 ; used to store X(a,Jb)
           if(ksym.eq.1.and.jSym.eq.1.and.DoSomeX)then
              iE = iE + (Nmax**2) * NumV
              Wab%SB(iSymp)%A1(1:Nmax**2 * NumV) => Wab%A0(iS:iE)
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

       ! first "free" position in Wab%A0
       iS = iE + 1

#ifdef _DEBUGPRINT_
       write(6,*) 'Batch ',iBatch,' of  ',nBatch,': NumV = ',NumV
       write(6,*) 'Total allocated:     ',kRdMem
       write(6,*) 'iE:                  ',iE
       write(6,*) 'iS:                ',iS
       write(6,*) 'JSYM:                ',JSYM
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

         VJ(1:NumV) => Wab%A0(iS:iE+NumV)
         VJ(:)=Zero

         DO iSymr=1,nSym
         IF(nBas(iSymr).ne.0.and.nOcc(iSymr,jDen).ne.0)THEN

         Naa = nBas(iSymr)*(nBas(iSymr)+1)/2

         CALL DGEMV_('T',Naa,NumV,
     &              ONE,Wab%SB(iSymr)%A2,Naa,
     &                  DLT(jDen)%SB(ISYMR)%A1,1,
     &              ONE,VJ,1)

         ENDIF
         End DO

C --- Contraction with the 2nd vector
C ---------------------------------------
C --- Fss{#J} <- Fss{#J} + FactC * sum_J L(ss,{#J})*V{#J}
C==========================================================
C
          DO iSyms=1,nSym
          IF(nBas(iSyms).ne.0)THEN

           Naa = nBas(iSyms)*(nBas(iSyms)+1)/2

           CALL DGEMV_('N',Naa,NumV,
     &              FactC(jDen),Wab%SB(iSyms)%A2,Naa,
     &                          VJ,1,
     &                      ONE,FLT(jDen)%SB(ISYMS)%A1,1)


          ENDIF
          End DO

          VJ=>Null()

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

         ISFSQ=0
         DO iSymr=1,nSym

          iSymr_Occ=0
          do jDen=1,nDen
             iSymr_Occ = iSymr_Occ + nOcc(iSymr,jDen)
          end do
          IF(nBas(iSymr).ne.0.and.iSymr_Occ.gt.0)THEN

             CALL CWTIME(TCREO1,TWREO1)

             LrJs(1:nBas(iSymr),1:NUMV,1:nBas(iSymr)) =>
     &           Wab%A0(iS:iE+nBas(iSymr)**2*NUMV)

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

               IF (nOcc(iSymr,jDen).ne.0) THEN

               CALL CWTIME(TC1X1,TW1X1)

               ndim3=NBAS(ISYMR)*NUMV*NBAS(ISYMR)
               XpJs(1:nDim3) => Wab%A0(1:nDim3)

C              Calculate intermediate:
C              X(p,Js) = Sum(q) D(p,q) * L(q,Js).
C              -----------------------------------


               CALL DGEMM_('N','N',
     &               NBAS(ISYMR),NUMV*NBAS(ISYMR),NBAS(ISYMR),
     &               ONE,DSQ(jDen)%SB(ISYMR)%A2,NBAS(ISYMR),
     &                   LrJs, NBAS(ISYMR),
     &               ZERO,XpJs,NBAS(ISYMR))

C              Calculate exchange contribution:
C              F(p,q) = F(p,q) - FactX Sum(sJ) L(sJ,p) * X(sJ,q).
C              --------------------------------------------------

               CALL DGEMM_('T','N',
     &                NBAS(ISYMR),NBAS(ISYMR),NBAS(ISYMR)*NUMV,
     &             -FactX(jDen),LrJs,NBAS(ISYMR)*NUMV,
     &                          XpJs,NBAS(ISYMR)*NUMV,
     &                      One,FSQ(jDen)%SB(ISYMR)%A2,NBAS(ISYMR))


             CALL CWTIME(TC1X2,TW1X2)
             texch(1) = texch(1) + (TC1X2 - TC1X1)
             texch(2) = texch(2) + (TW1X2 - TW1X1)

             XpJs=>Null()

             ENDIF   ! nOcc .ne. 0

          ENDIF  ! DoExchange(jDen)

         END DO  ! loop over the densities

        ENDIF  ! nbas.ne.0 & nOcc.ne.0

        LrJs=>Null()

        END DO  !loop over Fock mat symmetries

       END IF  ! jSym=1 and DoSomeX

C ****************** END DIAGONAL EXCHANGE **********************
C
C ************ "OFF-DIAGONAL" EXCHANGE CONTRIBUTION  ************
C
       IF (jSym.ne.1) THEN

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

               IF (nBas(iSyms)*nBas(iSymq)<=0) Cycle

               IF (ISYMQ.gt.ISYMS .and. iSkip(iSymq).ne.0 ) THEN

                  jS = jE + 1
                  jE = jE + nBas(iSymq)*NumV*nBas(iSyms)
                  LqJs%SB(ISYMQ)%A3(1:nBas(iSymq),1:NumV,1:nBas(iSyms))
     &             => Wab%A0(jS:jE)

              DO JVEC = 1,NUMV

                DO jS=1,NBAS(iSyms)

                   LqJs%SB(ISYMQ)%A3(:,JVEC,jS) =
     &               Wab%SB(ISYMQ)%A3(:,jS,JVEC)

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
C --- COMPUTE EXCHANGE FOR OFF-DIAGONAL VECTORS

            CALL CWTIME(TC2X1,TW2X1)

            DO ISYMG = 1,NSYM

               ISYMB = MULD2H(ISYMG,JSYM)

            IF(nBas(iSymb)*nBas(iSymg).ne.0)THEN

               ISYMD = ISYMG     !only total symmetric Density
               ISYMA = ISYMB     !block diagonal Fock Matrix

               !mem used to store the intermediate
               nd = NBAS(ISYMD)
               nb = NBAS(ISYMB)
               XdJb(1:nd*NUMV*nb) => Wab%A0(1:nd*NUMV*nb)

             IF (ISYMG.gt.ISYMB .and. iSkip(iSymg).ne.0) THEN


C -------------------------------
C --- F(a,b) = - D(g,d) * (ad|gb)
C -------------------------------
               if (nOcc(iSymg,jDen).ne.0) then

C              Calculate intermediate:
C              X(d,Jb) = Sum(g) D(d,g) * L(g,Jb).
C              ----------------------------------

               CALL DGEMM_('N','N',
     &                    NBAS(ISYMD),NUMV*NBAS(ISYMB),NBAS(ISYMG),
     &                    ONE,DSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMD),
     &                        LqJs%SB(ISYMG)%A3,NBAS(ISYMG),
     &                   ZERO,XdJb,NBAS(ISYMD))


C              F(a,b) = F(a,b) - Sum(dJ) L(dJ,a) * X(dJ,b).
C              -------------------------------------------

              CALL DGEMM_('T','N',
     &              NBAS(ISYMA),NBAS(ISYMB),NBAS(ISYMD)*NUMV,
     &              -FactX(jDen),LqJs%SB(ISYMG)%A3,NBAS(ISYMD)*NUMV,
     &                          XdJb,NBAS(ISYMD)*NUMV,
     &                      ONE,FSQ(jDen)%SB(ISYMB)%A2,NBAS(ISYMA))

              endif
C -------------------------------
C --- F(g,d) = - D(a,b) * (ad|gb)
C -------------------------------
               if (nOcc(iSyma,jDen).ne.0) then

C              Calculate intermediate:
C              X(gJ,b) = Sum(a) L(gJ,a)* D(a,b).
C              ----------------------------------

               CALL DGEMM_('N','N',
     &                  NBAS(ISYMG)*NUMV,NBAS(ISYMB),NBAS(ISYMA),
     &                  ONE,LqJs%SB(ISYMG)%A3,NBAS(ISYMG)*NUMV,
     &                      DSQ(jDen)%SB(ISYMB)%A2,NBAS(ISYMA),
     &                 ZERO,XdJb,NBAS(ISYMG)*NUMV)


C              F(g,d) = F(g,d) - Sum(Jb) X(g,Jb) * L(d,Jb).
C              -------------------------------------------

              CALL DGEMM_('N','T',
     &              NBAS(ISYMG),NBAS(ISYMD),NUMV*NBAS(ISYMB),
     &             -FactX(jDen),XdJb,NBAS(ISYMG),
     &                      LqJs%SB(ISYMG)%A3,NBAS(ISYMD),
     &               ONE,FSQ(jDen)%SB(ISYMG)%A2,NBAS(ISYMG))

              endif

             ENDIF

             XdJb=>Null()

            ENDIF
             LqJs%SB(ISYMG)%A3=>Null()

            END DO  ! loop over iorbital symmetries

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

      WRITE(6,'(6X,A)')'TEST PRINT FROM CHO_FOCKTWO.'
      WRITE(6,'(6X,A)')'***** EXCHANGE MATRIX ***** '
      Do jDen=1,nDen
          WRITE(6,'(6X,A,I2)')'DENSITY TYPE: ',jDen
          WRITE(6,'(6X,A,I2)')'DoExchange: ',DoExchange(jDen)
          WRITE(6,*)
      if(DoExchange(jDen))then
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call cho_output(FSQ(jDen)%SB(ISYM)%A2,1,NB,1,NB,NB,NB,1,6)
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
