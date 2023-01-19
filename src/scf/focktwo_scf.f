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
*               2002,2023, Roland Lindh                                *
************************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE FOCKTWO_scf(NSYM,NBAS,NFRO,KEEP,DLT,DSQ,FLT,nFlt,FSQ,
     &                       X1,nX1,X2,nX2,ExFac,nD,nBSQT)
      use RICD_Info, only: Do_DCCD
      IMPLICIT None
      Integer nSym, nFlt, nD, nBSQT
      Integer NBAS(8),NFRO(8), KEEP(8)
      Real*8 DLT(nFlt,nD)
      Real*8 DSQ(nBSQT,nD)
      Real*8 FSQ(nBSQT,nD)
      Real*8 FLT(nFlt,nD)
      Integer nX1, nX2
      Real*8 X1(nX1),X2(nX2)
      Real*8 ExFac

      Integer ISTLT(8),ISTSQ(8)
      Integer MUL, I, J
      Real*8  Factor, temp, temp_ab
      Integer IB, JB, IJB, IJS, LB, KB, KLB,  iOpt, IP, JQ, IPQ, IRC
      Integer LPQ, NPQ
      Integer IK, JK, KK, LK
      Integer IS, JS, KS, LS, ISYM
      Integer ISF, ISX, ISD
      Integer K1, K2, LSMAX, NB, NB2, NB3, NFI, NFJ, NFK, NFL
      Real*8, External:: DDot_
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
      If (Do_DCCD.and.NSYM/=1) Then
         Write (6,*) 'DCCD not implemented for nSym/=1'
         Call Abend()
      End If

      Factor=DBLE(nD)*0.5D0
      ISTSQ(:)=0
      ISTLT(:)=0

      IF (Do_DCCD) THEN
         If (NSYM/=1) Then
            Write (6,*) 'DCCD not implemented for nSym/=1'
            Call Abend()
         End If
         Call FOCKTWO_scf_DCCD()
      ELSE
         IF (NSYM==1) THEN
            Call FOCKTWO_scf_NoSym()
         ELSE
            Call FOCKTWO_scf_Sym()
         END IF
      END IF

      IF (IRC/=0) THEN
         WRITE(6,*)' Error return code IRC=',IRC
         WRITE(6,*)' from RDORD call, in FTWOI.'
         CALL Abend()
      END IF

c Accumulate the contributions
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        K1=ISTLT(ISYM)
        K2=ISTSQ(ISYM)
        DO IB=1,NB
          DO JB=1,IB

            FLT(K1+JB,1)=FLT(K1+JB,1)+FSQ(K2+JB,1)
            if(nD==2) then
             FLT(K1+JB,2)=FLT(K1+JB,2)+FSQ(K2+JB,2)
            endif
#ifdef _DEBUGPRINT_
         if(nD==1)then
          write (6,'(a,i5,a,f12.6)') 'Flt(',K1+JB,',1)=',FLT(K1+JB,1)
          else
          write (6,'(a,i5,a,2f12.6)') 'Flt(',K1+JB,',:)=',
     &               FLT(K1+JB,1),FLT(K1+JB,2)
          endif
#endif

          END DO
          K1=K1+IB
          K2=K2+NB
        END DO
      END DO
*
      Call GADSum(Flt,nFlt*nD)
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
          CALL TRIPRT(' ',' ',FLT(ISTLTT,1),NB)
          if(nD==2) then
          CALL TRIPRT(' ',' ',FLT(ISTLTT,2),NB)
          endif
          ISTLTT=ISTLTT+NB*(NB+1)/2
        END IF
      END DO
      WRITE(6,'(6X,A)')'----------------------------'
#endif


      CONTAINS

      Subroutine FOCKTWO_scf_Sym()

      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NB2=NB*NB
        NB3=(NB2+NB)/2
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NB2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NB3
      END DO

!     Loop over the symmetry blocks (IS,JS|KS,LS)

      DO IS=1,NSYM
        IB=NBAS(IS)
        IK=KEEP(IS)
        NFI=NFRO(IS)
        DO JS=1,IS
          JB=NBAS(JS)
          JK=KEEP(JS)
          NFJ=NFRO(JS)
          IJS=MUL(IS,JS)
          IJB=IB*JB
          IF( IS.EQ.JS ) IJB=(IB*(IB+1))/2
          DO KS=1,IS
            KB=NBAS(KS)
            KK=KEEP(KS)
            NFK=NFRO(KS)
            LSMAX=KS
            IF ( KS.EQ.IS ) LSMAX=JS
            LS=MUL(IJS,KS)
            IF(LS.GT.LSMAX) CYCLE
            LB=NBAS(LS)
            LK=KEEP(LS)
            NFL=NFRO(LS)
            KLB=KB*LB
            IF( KS.EQ.LS ) KLB=(KB*(KB+1))/2
C INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
            IF((IK+JK+KK+LK).NE.0) CYCLE
C NO FROZEN ORBITALS?
            IF((NFI+NFJ+NFK+NFL).EQ.0) CYCLE
C NO BASIS FUNCTIONS?
            IF((IJB*KLB).EQ.0 ) CYCLE

! Process the different symmetry cases

            IF ( IS.EQ.JS .AND. IS.EQ.KS ) THEN
c CASE 1: Integrals are of symmetry type (II/II)
c Coulomb and exchange terms need to be accumulated
c Option code 1: Begin reading at first integral.
c NPQ: Nr of submatrices in buffer X1.
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO IP=1,IB
                  DO JQ=1,IP
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
                      IF(IRC.GT.1) Return
c Option code 2: Continue reading at next integral.
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    ISF=ISTLT(IS)+LPQ
                    ISD=ISTLT(IS)+1
                    TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD,1),1)
                    FLT(ISF,1)=FLT(ISF,1)+TEMP
                    if(nD==2) then
                     TEMP_ab=DDOT_(KLB,X1(ISX),1,DLT(ISD,2),1)
                     FLT(ISF,1)=FLT(ISF,1)+TEMP_ab
                     FLT(ISF,2)=FLT(ISF,1)
                    endif
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)') '00 Flt(',isf,',1)=',FLT(ISF,1)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)') '00 Flt(',isf,',2)=',FLT(ISF,2)
          endif
#endif
                    CALL SQUARE (X1(ISX),X2(1),1,KB,LB)
                    ISF=ISTSQ(IS)+(JQ-1)*JB+1
                    ISD=ISTSQ(IS)+(IP-1)*IB+1
c
                 if(nD==1) then
                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                           DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                 else
                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                           DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)

                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                           DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
                 endif
                    IF ( IP.NE.JQ ) THEN
                      ISF=ISTSQ(IS)+(IP-1)*IB+1
                      ISD=ISTSQ(IS)+(JQ-1)*JB+1
c
                 if(nD==1) then
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                             DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                 else
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                             DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                             DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
                 endif
                    ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
          endif
#endif

                  END DO  ! JQ
                END DO    ! IP

              ELSE IF ( IS.EQ.JS .AND. IS.NE.KS ) THEN
c CASE 2: Integrals are of symmetry type (II/JJ)
c Coulomb terms need to be accumulated only
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO IP=1,IB
                  DO JQ=1,IP
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
                      IF(IRC.GT.1) Return
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    IF ( NFI.NE.0 ) THEN
                      ISF=ISTLT(KS)+1
                      ISD=ISTLT(IS)+LPQ
                      TEMP=DLT(ISD,1)
                    if(nD==2) then
                      TEMP=DLT(ISD,1)+DLT(ISD,2)
                    endif
                      CALL DAXPY_(KLB,TEMP,X1(ISX),1,FLT(ISF,1),1)
                    if(nD==2) then
                      CALL DAXPY_(KLB,TEMP,X1(ISX),1,FLT(ISF,2),1)
                    endif
                    ENDIF
                    IF ( NFK.NE.0 ) THEN
                      ISF=ISTLT(IS)+LPQ
                      ISD=ISTLT(KS)+1
                      TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD,1),1)
                      FLT(ISF,1)=FLT(ISF,1)+TEMP
                      if (nD==2) then
                         TEMP_ab=DDOT_(KLB,X1(ISX),1,DLT(ISD,2),1)
                         FLT(ISF,1)=FLT(ISF,1)+TEMP_ab
                         FLT(ISF,2)=FLT(ISF,1)
                      endif
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)') '02 Flt(',isf,',1)=',FLT(ISF,1)
#endif

                    ENDIF
                  END DO! JQ
                END DO  ! IP
              ELSE IF ( IS.EQ.KS .AND. JS.EQ.LS ) THEN
c CASE 3: Integrals are of symmetry type (IJ/IJ)
c Exchange terms need to be accumulated only
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO IP=1,IB
                  DO JQ=1,JB
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
                      IF(IRC.GT.1) Return
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    IF ( NFI.NE.0 ) THEN
                      ISD=ISTSQ(IS)+(IP-1)*IB+1
                      ISF=ISTSQ(JS)+(JQ-1)*JB+1
                      if(nD==1) then
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                      else
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
                      endif
                    ENDIF
                    IF ( NFJ.NE.0 ) THEN
                      ISD=ISTSQ(JS)+(JQ-1)*JB+1
                      ISF=ISTSQ(IS)+(IP-1)*IB+1
                      if(nD==1) then
                      CALL DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                      else
                      CALL DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                      CALL DGEMV_('T',LB,KB,-factor*ExFac,X1(ISX),LB,
     &                             DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)

                      endif
                    ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f20.6)')
     &          ('03 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
          if(nD==2) then
          write (6,'(a,i5,a,f20.6)')
     &          ('03 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
          endif
#endif

                  END DO ! JQ
                END DO   ! IP

            ENDIF

            END DO! KS
         END DO   ! JS
      END DO      ! IS

      END SUBROUTINE FOCKTWO_scf_Sym

      Subroutine FOCKTWO_scf_NoSym()

      IS=1
      IB=NBAS(IS)
      IK=KEEP(IS)
      NFI=NFRO(IS)

      JS=1
      JB=NBAS(JS)
      JK=KEEP(JS)
      NFJ=NFRO(JS)
      IJS=MUL(IS,JS)
      IJB=(IB*(IB+1))/2

      KS=1
      KB=NBAS(KS)
      KK=KEEP(KS)
      NFK=NFRO(KS)
      LSMAX=JS

      LS=1
      LB=NBAS(LS)
      LK=KEEP(LS)
      NFL=NFRO(LS)
      KLB=(KB*(KB+1))/2

C INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?

      IF((IK+JK+KK+LK)/=0) Return
C NO FROZEN ORBITALS?
      IF((NFI+NFJ+NFK+NFL)==0) Return
C NO BASIS FUNCTIONS?
      IF((IJB*KLB)==0) Return

! Process the different symmetry cases

c CASE 1: Integrals are of symmetry type (II/II)
c Coulomb and exchange terms need to be accumulated
c Option code 1: Begin reading at first integral.
c NPQ: Nr of submatrices in buffer X1.
      IOPT=1
      LPQ=0
      IPQ=0
      NPQ=0
      DO IP=1,IB
         DO JQ=1,IP
            IPQ=IPQ+1
            LPQ=LPQ+1
            IF ( IPQ.GT.NPQ ) THEN
               CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
               IF(IRC.GT.1) Return
c Option code 2: Continue reading at next integral.
               IOPT=2
               IPQ=1
            ENDIF
            ISX=(IPQ-1)*KLB+1
            ISF=LPQ
            ISD=1
            TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD,1),1)
            FLT(ISF,1)=FLT(ISF,1)+TEMP
            if(nD==2) then
              TEMP_ab=DDOT_(KLB,X1(ISX),1,DLT(ISD,2),1)
              FLT(ISF,1)=FLT(ISF,1)+TEMP_ab
              FLT(ISF,2)=FLT(ISF,1)
            endif
#ifdef _DEBUGPRINT_
            write (6,'(a,i5,a,f12.6)') '00 Flt(',isf,',1)=',FLT(ISF,1)
            if(nD==2) then
            write (6,'(a,i5,a,f12.6)') '00 Flt(',isf,',2)=',FLT(ISF,2)
            endif
#endif
            CALL SQUARE (X1(ISX),X2(1),1,KB,LB)
            ISF=(JQ-1)*JB+1
            ISD=(IP-1)*IB+1
c
            if(nD==1) then
              CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                    DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
            else
              CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                    DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)

              CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                    DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
            endif
            IF ( IP.NE.JQ ) THEN
               ISF=(IP-1)*IB+1
               ISD=(JQ-1)*JB+1
c
               if(nD==1) then
                 CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                       DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
               else
                 CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                       DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                 CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,
     &                       DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
               endif
            ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
          endif
#endif

         END DO  ! JQ
      END DO    ! IP

      END SUBROUTINE FOCKTWO_scf_NoSym

      Subroutine FOCKTWO_scf_DCCD()
      use stdalloc, only: mma_allocate, mma_deallocate
      Integer nData
      Integer, Allocatable:: Basis_IDs(:,:)
      Logical Found
      Integer IP, JQ, IPQ, KR, LS, IRS


      Call Init_GetInt(IRC)
      Call Qpg_iArray('Basis IDs',Found,nData)
      Call mma_allocate(Basis_IDs,4,nData/4,Label='Basis_IDs')
      call Get_iArray('Basis IDs',Basis_IDs,nData)

      IS=1
      IB=NBAS(IS)
      IK=KEEP(IS)
      NFI=NFRO(IS)

      IJB=(IB*(IB+1))/2

      Call Get_Int_Open(IS,IS,IS,IS)

C INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?

      IF(IK/=0) Return
C NO FROZEN ORBITALS?
      IF(NFI==0) Return
C NO BASIS FUNCTIONS?
      IF(IJB==0) Return

! Process the different symmetry cases

c CASE 1: Integrals are of symmetry type (II/II)
c Coulomb and exchange terms need to be accumulated
c Option code 1: Begin reading at first integral.
c NPQ: Nr of submatrices in buffer X1.
      IOPT=1
      LPQ=0
      IPQ=0
      NPQ=0
      DO IP=1,IB
         DO JQ=1,IP
            IPQ=IPQ+1
            LPQ=LPQ+1
            IF ( IPQ.GT.NPQ ) THEN
               CALL Get_Int_DCCD(IRC,IOPT,IS,IS,IS,IS,X1,IJB+1,NPQ)
               IF(IRC.GT.1) Return
! Option code 2: Continue reading at next integral.
               IOPT=2
               IPQ=1
            ENDIF
! Skip processing (P,Q|... if they do not share the same center
            IF (Basis_IDs(1,IP)/=Basis_IDs(1,JQ)) CYCLE
! Do the Coulomb contribution
            IF (nD==1) Then
               TEMP=0.0D0
               DO KR=1,IB
                  DO LS=1,KR
                     IRS=KR*(KR-1)/2 + LS
                     TEMP=TEMP+X1(IRS)*DLT(IRS,1)
                  END DO
               END DO
               FLT(LPQ,1)=FLT(LPQ,1)+TEMP
            ELSE
               TEMP=0.0D0
               TEMP_ab=0.0D0
               DO KR=1,IB
                  DO LS=1,KR
                     IRS=KR*(KR-1)/2 + LS
                     TEMP=TEMP+X1(IRS)*DLT(IRS,1)
                     TEMP_ab=TEMP_ab+X1(IRS)*DLT(IRS,2)
                  END DO
               END DO
               FLT(LPQ,1)=FLT(LPQ,1)+TEMP+TEMP_ab
               FLT(LPQ,2)=FLT(LPQ,1)
            END IF
#ifdef _DEBUGPRINT_
            write (6,'(a,i5,a,f12.6)') '00 Flt(',LPQ,',1)=',FLT(LPQ,1)
            if(nD==2) then
            write (6,'(a,i5,a,f12.6)') '00 Flt(',LPQ,',2)=',FLT(LPQ,2)
            endif
#endif
! Do the exchange contribution
            CALL SQUARE (X1(:),X2(:),1,IB,IB)
            ISF=(JQ-1)*IB+1
            ISD=(IP-1)*IB+1
c
            if(nD==1) then
              CALL DGEMV_('N',IB,IB,-Factor*ExFac,X2(1),IB,
     &                    DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
            else
              CALL DGEMV_('N',IB,IB,-Factor*ExFac,X2(1),IB,
     &                    DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)

              CALL DGEMV_('N',IB,IB,-Factor*ExFac,X2(1),IB,
     &                    DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
            endif
            IF ( IP.NE.JQ ) THEN
               ISF=(IP-1)*IB+1
               ISD=(JQ-1)*IB+1
c
               if(nD==1) then
                 CALL DGEMV_('N',IB,IB,-Factor*ExFac,X2(1),IB,
     &                       DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
               else
                 CALL DGEMV_('N',IB,IB,-Factor*ExFac,X2(1),IB,
     &                       DSQ(ISD,1),1,1.0D0,FSQ(ISF,1),1)
                 CALL DGEMV_('N',IB,IB,-Factor*ExFac,X2(1),IB,
     &                       DSQ(ISD,2),1,1.0D0,FSQ(ISF,2),1)
               endif
            ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,Ib)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,Ib)
          endif
#endif

         END DO  ! JQ
      END DO    ! IP

      Call Get_Int_Close()
      Call mma_deallocate(Basis_IDs)

      END SUBROUTINE FOCKTWO_scf_DCCD

      END SUBROUTINE FOCKTWO_scf
