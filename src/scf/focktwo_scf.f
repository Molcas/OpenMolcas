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
* Copyright (C) Markus P. Fuelscher                                    *
*               1992, Per Ake Malmqvist                                *
*               2002, Roland Lindh                                     *
************************************************************************
      SUBROUTINE FOCKTWO_scf(NSYM,NBAS,NFRO,KEEP,
     &                   DLT,DSQ,FLT,nFlt,FSQ,LBUF,X1,X2,ExFac,iUHF,
     &                   DLT_ab,DSQ_ab,FLT_ab,FSQ_ab)
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 FSQ(*),FLT(nFlt),DSQ(*),DLT(*),X1(*),X2(*)
      Real*8 DLT_ab(*),DSQ_ab(*),FLT_ab(*),FSQ_ab(*)
      Integer ISTLT(8),ISTSQ(8),KEEP(8),NBAS(8),NFRO(8)
      Logical myDebug
c
c This routine has been nicked from the MOTRA package. It was
c originally written by Marcus Fuelscher, and has been slightly
c modified by P-A Malmqvist 1992-12-04.
c Further modifications by RL 2002-8-30.
c Purpose: Return FLT, which is the effective one-electron
c Hamiltonian when frozen orbitals are removed. It is the
c same as a closed-shell Fock matrix computed using the
c two-electron contribution only from frozen orbitals.
c FLT is returned as symmetry-blocked lower triangles. FSQ
c contains the same matrix, as symmetry-blocked square matrices.
c DSQ and DLT are computed in the calling routine.
c It is assumed that the SEWARD integral file was opened before
c call to this routine.
c ISTSQ, ISTLT: Offsets to symmetry blocks of FSQ, FLT etc.
c
************************************************************************
*                                                                      *
      MUL(I,J)=IEOR(I-1,J-1)+1
*                                                                      *
************************************************************************
       myDebug=.false.
       Factor=0.5D0
        if(iUHF.eq.1) Factor=1.0D0
#ifdef _BUGPRINT_
c      myDebug=.true. ! very extensive print out
#endif
      ISTSQ(1)=0
      ISTLT(1)=0
      DO 20 ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NB2=NB*NB
        NB3=(NB2+NB)/2
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NB2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NB3
20    CONTINUE

      DO 110 IS=1,NSYM
        IB=NBAS(IS)
        IK=KEEP(IS)
        NFI=NFRO(IS)
        DO 120 JS=1,IS
          JB=NBAS(JS)
          JK=KEEP(JS)
          NFJ=NFRO(JS)
          IJS=MUL(IS,JS)
          IJB=IB*JB
          IF( IS.EQ.JS ) IJB=(IB*(IB+1))/2
          DO 130 KS=1,IS
            KB=NBAS(KS)
            KK=KEEP(KS)
            NFK=NFRO(KS)
            LSMAX=KS
            IF ( KS.EQ.IS ) LSMAX=JS
            LS=MUL(IJS,KS)
            IF(LS.GT.LSMAX) GOTO 130
            LB=NBAS(LS)
            LK=KEEP(LS)
            NFL=NFRO(LS)
            KLB=KB*LB
            IF( KS.EQ.LS ) KLB=(KB*(KB+1))/2
C INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
            IF((IK+JK+KK+LK).NE.0) GOTO 130
C NO FROZEN ORBITALS?
            IF((NFI+NFJ+NFK+NFL).EQ.0) GOTO 130
C NO BASIS FUNCTIONS?
            IF((IJB*KLB).EQ.0 ) GOTO 130

            IF ( IS.EQ.JS .AND. IS.EQ.KS ) THEN
c CASE 1: Integrals are of symmetry type (II/II)
c Coulomb and exchange terms need to be accumulated
c Option code 1: Begin reading at first integral.
c NPQ: Nr of submatrices in buffer X1.
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO 200 IP=1,IB
                  DO 201 JQ=1,IP
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
                      IF(IRC.GT.1) GOTO 999
c Option code 2: Continue reading at next integral.
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    ISF=ISTLT(IS)+LPQ
                    ISD=ISTLT(IS)+1
                    TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD),1)
                    FLT(ISF)=FLT(ISF)+TEMP
                    if(iUHF.eq.1) then
                     TEMP_ab=DDOT_(KLB,X1(ISX),1,DLT_ab(ISD),1)
                     FLT(ISF)=FLT(ISF)+TEMP_ab
                     FLT_ab(ISF)=FLT(ISF)
                    endif
        if(myDebug) then
          write (6,'(a,i5,a,f12.6)') '00 Flt(',isf,')=',FLT(ISF)
          if(iUHF.eq.1) then
          write (6,'(a,i5,a,f12.6)') '00 Flt_ab(',isf,')=',FLT_ab(ISF)
          endif
        endif
                    CALL SQUARE (X1(ISX),X2(1),1,KB,LB)
                    ISF=ISTSQ(IS)+(JQ-1)*JB+1
                    ISD=ISTSQ(IS)+(IP-1)*IB+1
c
                 if(iUHF.eq.0) then
*                    CALL DGEMX  (KB,LB,-Factor*ExFac,X2(1),KB,
*     &                           DSQ(ISD),1,FSQ(ISF),1)
                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                           DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                 else
*                    CALL DGEMX  (KB,LB,-Factor*ExFac,X2(1),KB,
*     &                           DSQ(ISD),1,FSQ(ISF),1)
*
*                    CALL DGEMX  (KB,LB,-Factor*ExFac,X2(1),KB,
*     &                           DSQ_ab(ISD),1,FSQ_ab(ISF),1)
                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                           DSQ(ISD),1,1.0D0,FSQ(ISF),1)

                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                           DSQ_ab(ISD),1,1.0D0,FSQ_ab(ISF),1)
                 endif
                    IF ( IP.NE.JQ ) THEN
                      ISF=ISTSQ(IS)+(IP-1)*IB+1
                      ISD=ISTSQ(IS)+(JQ-1)*JB+1
c
                 if(iUHF.eq.0) then
*                      CALL DGEMX  (KB,LB,-Factor*ExFac,X2(1),KB,
*     &                             DSQ(ISD),1,FSQ(ISF),1)
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                             DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                 else
*                      CALL DGEMX  (KB,LB,-Factor*ExFac,X2(1),KB,
*     &                             DSQ(ISD),1,FSQ(ISF),1)
*                      CALL DGEMX  (KB,LB,-Factor*ExFac,X2(1),KB,
*     &                             DSQ_ab(ISD),1,FSQ_ab(ISF),1)
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                             DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                             DSQ_ab(ISD),1,1.0D0,FSQ_ab(ISF),1)
                 endif
                    ENDIF
        if(myDebug) then
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,')=',FSQ(ISF+ivv-1),ivv=1,kb)
          if(iUHF.eq.1) then
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq_ab(',isf+ivv-1,')=',FSQ_ab(ISF+ivv-1),ivv=1,kb)
          endif
        endif

201               CONTINUE
200             CONTINUE

              ELSE IF ( IS.EQ.JS .AND. IS.NE.KS ) THEN
c CASE 2: Integrals are of symmetry type (II/JJ)
c Coulomb terms need to be accumulated only
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO 210 IP=1,IB
                  DO 211 JQ=1,IP
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
                      IF(IRC.GT.1) GOTO 999
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    IF ( NFI.NE.0 ) THEN
                      ISF=ISTLT(KS)+1
                      ISD=ISTLT(IS)+LPQ
                      TEMP=DLT(ISD)
                    if(iUHF.eq.1) then
                      TEMP=DLT(ISD)+DLT_ab(ISD)
                    endif
                      CALL DAXPY_(KLB,TEMP,X1(ISX),1,FLT(ISF),1)
                    if(iUHF.eq.1) then
                      CALL DAXPY_(KLB,TEMP,X1(ISX),1,FLT_ab(ISF),1)
                    endif
                    ENDIF
                    IF ( NFK.NE.0 ) THEN
                      ISF=ISTLT(IS)+LPQ
                      ISD=ISTLT(KS)+1
                      TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD),1)
                      FLT(ISF)=FLT(ISF)+TEMP
                if(iUHF.eq.1) then
                      TEMP_ab=DDOT_(KLB,X1(ISX),1,DLT_ab(ISD),1)
                      FLT(ISF)=FLT(ISF)+TEMP_ab
                      FLT_ab(ISF)=FLT(ISF)
                endif
        if(myDebug) then
          write (6,'(a,i5,a,f12.6)') '02 Flt(',isf,')=',FLT(ISF)
        endif

                    ENDIF
211               CONTINUE
210             CONTINUE
              ELSE IF ( IS.EQ.KS .AND. JS.EQ.LS ) THEN
c CASE 3: Integrals are of symmetry type (IJ/IJ)
c Exchange terms need to be accumulated only
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO 220 IP=1,IB
                  DO 221 JQ=1,JB
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
                      IF(IRC.GT.1) GOTO 999
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    IF ( NFI.NE.0 ) THEN
                      ISD=ISTSQ(IS)+(IP-1)*IB+1
                      ISF=ISTSQ(JS)+(JQ-1)*JB+1
                      if(iUHF.eq.0) then
*                      CALL DGEMX  (LB,KB,-Factor*ExFac,X1(ISX),LB,
*     &                             DSQ(ISD),1,FSQ(ISF),1)
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                      else
*                      CALL DGEMX  (LB,KB,-Factor*ExFac,X1(ISX),LB,
*     &                             DSQ(ISD),1,FSQ(ISF),1)
*                      CALL DGEMX  (LB,KB,-Factor*ExFac,X1(ISX),LB,
*     &                             DSQ_ab(ISD),1,FSQ_ab(ISF),1)
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ_ab(ISD),1,1.0D0,FSQ_ab(ISF),1)
                      endif
                    ENDIF
                    IF ( NFJ.NE.0 ) THEN
                      ISD=ISTSQ(JS)+(JQ-1)*JB+1
                      ISF=ISTSQ(IS)+(IP-1)*IB+1
                      if(iUHF.eq.0) then
*                      CALL DGEMTX (LB,KB,-Factor*ExFac,X1(ISX),LB,
*     &                             DSQ(ISD),1,FSQ(ISF),1)
                      CALL DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                      else
*                      CALL DGEMTX (LB,KB,-Factor*ExFac,X1(ISX),LB,
*     &                             DSQ(ISD),1,FSQ(ISF),1)
*                      CALL DGEMTX (LB,KB,-factor*ExFac,X1(ISX),LB,
*     &                             DSQ_ab(ISD),1,FSQ_ab(ISF),1)
                      CALL DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD),1,1.0D0,FSQ(ISF),1)
                      CALL DGEMV_('T',LB,KB,-factor*ExFac,X1(ISX),LB,
     &                             DSQ_ab(ISD),1,1.0D0,FSQ_ab(ISF),1)

                      endif
                    ENDIF
        if(myDebug) then
          write (6,'(a,i5,a,f20.6)')
     &          ('03 Fsq(',isf+ivv-1,')=',FSQ(ISF+ivv-1),ivv=1,kb)
          if(iUHF.eq.1) then
          write (6,'(a,i5,a,f20.6)')
     &          ('03 Fsq_ab(',isf+ivv-1,')=',FSQ_ab(ISF+ivv-1),ivv=1,kb)
          endif
        endif

221               CONTINUE
220             CONTINUE
            ENDIF
130       CONTINUE
120     CONTINUE
110   CONTINUE

c Accumulate the contributions
      DO 300 ISYM=1,NSYM
        NB=NBAS(ISYM)
        K1=ISTLT(ISYM)
        K2=ISTSQ(ISYM)
        DO 310 IB=1,NB
          DO 315 JB=1,IB
c         write (6,'(a,i5,a,f12.6)') ' >> Flt(',K1+JB,')=',FLT(K1+JB)
c         write (6,'(a,i5,a,f12.6)') ' >> Fsq(',K2+JB,')=',FSQ(K2+JB)

            FLT(K1+JB)=FLT(K1+JB)+FSQ(K2+JB)
            if(iUHF.eq.1) then
             FLT_ab(K1+JB)=FLT_ab(K1+JB)+FSQ_ab(K2+JB)
            endif
        if(myDebug) then
         if(iUHF.eq.0)then
          write (6,'(a,i5,a,f12.6)') 'Flt(',K1+JB,')=',FLT(K1+JB)
          else
          write (6,'(a,i5,a,2f12.6)') 'Flt_ab(',K1+JB,')=',
     &               FLT(K1+JB),FLT_ab(K1+JB)
          endif
        endif

315       CONTINUE
          K1=K1+IB
          K2=K2+NB
310     CONTINUE
300   CONTINUE
*
      Call GADSum(Flt,nFlt)
      If (iUHF.eq.1) Call GADSum(Flt_ab,nFlt)
*
c Print the Fock-matrix
#ifdef _DEBUGPRINT_
      WRITE(6,'(6X,A)')'TEST PRINT FROM FTWOI.'
      WRITE(6,'(6X,A)')'FROZEN FOCK MATRIX IN AO BASIS:'
      ISTLTT=1
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT(' ',' ',FLT(ISTLTT),NB)
          if(iUHF.eq.1) then
          CALL TRIPRT(' ',' ',FLT_ab(ISTLTT),NB)
          endif
          ISTLTT=ISTLTT+NB*(NB+1)/2
        END IF
      END DO
      WRITE(6,'(6X,A)')'----------------------------'
#endif

      RETURN
 999  CONTINUE
      WRITE(6,*)' Error return code IRC=',IRC
      WRITE(6,*)' from RDORD call, in FTWOI.'
      CALL Abend
      END
