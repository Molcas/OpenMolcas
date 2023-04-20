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
! Copyright (C) Markus P. Fuelscher                                    *
!               1992, Per Ake Malmqvist                                *
!               2002,2023, Roland Lindh                                *
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE FOCKTWO_scf(NSYM,NBAS,NFRO,KEEP,DLT,DSQ,FLT,nFlt,FSQ,X1,nX1,X2,nX2,ExFac,nD,nBSQT)
      use RICD_Info, only: Do_DCCD
      use Constants, only: Zero, Half, One
!     IMPLICIT None
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
      Integer IB, JB, IJB, IJS, LB, KB, KLB,  iOpt, IPQ, IRC
      Integer LPQ, NPQ
      Integer IK, JK, KK, LK
      Integer IS, JS, KS, LS, ISYM
      Integer ISF, ISX, ISD
      Integer K1, K2, LSMAX, NB, NB2, NB3, NFI, NFJ, NFK, NFL
      Real*8, External:: DDot_
!
! This routine has been nicked from the MOTRA package. It was
! originally written by Marcus Fuelscher, and has been slightly
! modified by P-A Malmqvist 1992-12-04.
! Further modifications by RL 2002-8-30.
! Purpose: Return FLT, which is the effective one-electron
! Hamiltonian when frozen orbitals are removed. It is the
! same as a closed-shell Fock matrix computed using the
! two-electron contribution only from frozen orbitals.
! FLT is returned as symmetry-blocked lower triangles. FSQ
! contains the same matrix, as symmetry-blocked square matrices.
! DSQ and DLT are computed in the calling routine.
! It is assumed that the SEWARD integral file was opened before
! call to this routine.
! ISTSQ, ISTLT: Offsets to symmetry blocks of FSQ, FLT etc.
!
!***********************************************************************
!                                                                      *
      MUL(I,J)=IEOR(I-1,J-1)+1
!                                                                      *
!***********************************************************************
      If (Do_DCCD.and.NSYM/=1) Then
         Write (6,*) 'DCCD not implemented for nSym/=1'
         Call Abend()
      End If

      Factor=DBLE(nD)*Half
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

! Accumulate the contributions
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
          write (6,'(a,i5,a,2f12.6)') 'Flt(',K1+JB,',:)=',FLT(K1+JB,1),FLT(K1+JB,2)
          endif
#endif

          END DO
          K1=K1+IB
          K2=K2+NB
        END DO
      END DO
!
      Call GADSum(Flt,nFlt*nD)
!
! Print the Fock-matrix
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
      Integer ISYM, IS, JS, KS, IP, JQ

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
! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
            IF((IK+JK+KK+LK).NE.0) CYCLE
! NO FROZEN ORBITALS?
            IF((NFI+NFJ+NFK+NFL).EQ.0) CYCLE
! NO BASIS FUNCTIONS?
            IF((IJB*KLB).EQ.0 ) CYCLE

! Process the different symmetry cases

            IF ( IS.EQ.JS .AND. IS.EQ.KS ) THEN
! CASE 1: Integrals are of symmetry type (II/II)
! Coulomb and exchange terms need to be accumulated
! Option code 1: Begin reading at first integral.
! NPQ: Nr of submatrices in buffer X1.
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
! Option code 2: Continue reading at next integral.
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
!
                 if(nD==1) then
                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                 else
                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)

                    CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
                 endif
                    IF ( IP.NE.JQ ) THEN
                      ISF=ISTSQ(IS)+(IP-1)*IB+1
                      ISD=ISTSQ(IS)+(JQ-1)*JB+1
!
                 if(nD==1) then
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                 else
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                      CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
                 endif
                    ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
          endif
#endif

                  END DO  ! JQ
                END DO    ! IP

              ELSE IF ( IS.EQ.JS .AND. IS.NE.KS ) THEN
! CASE 2: Integrals are of symmetry type (II/JJ)
! Coulomb terms need to be accumulated only
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
! CASE 3: Integrals are of symmetry type (IJ/IJ)
! Exchange terms need to be accumulated only
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
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                      else
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                      CALL DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
                      endif
                    ENDIF
                    IF ( NFJ.NE.0 ) THEN
                      ISD=ISTSQ(JS)+(JQ-1)*JB+1
                      ISF=ISTSQ(IS)+(IP-1)*IB+1
                      if(nD==1) then
                      CALL DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                      else
                      CALL DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                      CALL DGEMV_('T',LB,KB,-factor*ExFac,X1(ISX),LB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)

                      endif
                    ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f20.6)') ('03 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
          if(nD==2) then
          write (6,'(a,i5,a,f20.6)') ('03 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
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
      Integer IP, JQ

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

! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?

      IF((IK+JK+KK+LK)/=0) Return
! NO FROZEN ORBITALS?
      IF((NFI+NFJ+NFK+NFL)==0) Return
! NO BASIS FUNCTIONS?
      IF((IJB*KLB)==0) Return

! Process the different symmetry cases

! CASE 1: Integrals are of symmetry type (II/II)
! Coulomb and exchange terms need to be accumulated
! Option code 1: Begin reading at first integral.
! NPQ: Nr of submatrices in buffer X1.
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
! Option code 2: Continue reading at next integral.
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
            CALL SQUARE (X1(ISX),X2(:),1,KB,LB)
            ISF=(JQ-1)*JB+1
            ISD=(IP-1)*IB+1
!
            if(nD==1) then
              CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
            else
              CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)

              CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
            endif
            IF ( IP.NE.JQ ) THEN
               ISF=(IP-1)*IB+1
               ISD=(JQ-1)*JB+1
!
               if(nD==1) then
                 CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
               else
                 CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                 CALL DGEMV_('N',KB,LB,-Factor*ExFac,X2(1),KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
               endif
            ENDIF
#ifdef _DEBUGPRINT_
          write (6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
          endif
#endif

         END DO  ! JQ
      END DO    ! IP

      END SUBROUTINE FOCKTWO_scf_NoSym

      Subroutine FOCKTWO_scf_DCCD()
      use stdalloc, only: mma_allocate, mma_deallocate
      use GetInt_mod, only: Basis_IDs, ID_IP, hash_table, lists, I, NumV, nVec, Vec2, NumCho, LuCVec, nPQ
      use Index_Functions, only: iTri
      Integer nData
      Logical Found
      Integer IS, IB, IP, JQ, IPQ, KR, LS, IRS
      Integer ISR, ISP, IRQ, IRP, ISQ
      Integer IP_, JQ_, KR_, LS_
      Integer iVec1, J

      IS=1
      IB=NBAS(IS)
      IK=KEEP(IS)
      NFI=NFRO(IS)

      Call Init_GetInt(IRC)
      Call Qpg_iArray('Basis IDs',Found,nData)
      Call mma_allocate(Basis_IDs,4,nData/4,Label='Basis_IDs')
      call Get_iArray('Basis IDs',Basis_IDs,nData)

      call mma_allocate(hash_table,Size(Basis_IDs,2),Label='hash_table')
      Do IP=1,IB
         hash_table(IP)=IP
      End Do
      Do IP=1,IB-1
         KR=hash_table(IP)
         Do JQ=IP+1,IB
            LS=hash_table(JQ)
            If (Basis_IDs(1,KR)>Basis_IDs(1,LS)) Then
               IPQ=KR
               KR=LS
               LS=IPQ
               hash_table(IP)=KR
               hash_table(JQ)=LS
            End If
         End Do
      End Do
      IP=1
      ID_IP=Basis_IDs(1,hash_table(IP))
      LS=1
      Do IP=2,IB
         If (ID_IP/=Basis_IDs(1,hash_table(IP))) Then
            ID_IP=Basis_IDs(1,hash_table(IP))
            LS=LS+1
         End If
      End Do
      Call mma_allocate(lists,4,LS)
      IP=1
      ID_IP=Basis_IDs(1,hash_table(IP))
      LS=1
      lists(1,LS)=1
      lists(2,LS)=ID_IP
      lists(3,LS)=IP
      lists(4,LS)=IP
      Do IP=2,IB
         If (ID_IP/=Basis_IDs(1,hash_table(IP))) Then
            ID_IP=Basis_IDs(1,hash_table(IP))
            LS=LS+1
            lists(1,LS)=1
            lists(2,LS)=ID_IP
            lists(3,LS)=IP
            lists(4,LS)=IP
         Else
            lists(1,LS)=lists(1,LS)+1
            lists(4,LS)=IP
         End If
      End Do

      IJB=(IB*(IB+1))/2
      If (nX1<IJB) Then
         Write (6,*) 'FOCKTWO_SCF_DCCD: nX1<IJB'
         Call Abend()
      End If
      Call Get_Int_Open(IS,IS,IS,IS)

! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?

      IF(IK/=0) Return
! NO FROZEN ORBITALS?
      IF(NFI==0) Return
! NO BASIS FUNCTIONS?
!     IF(IJB==0) Return

! Process the different symmetry cases
      do iVec1=1,NumCho(1),nVec
         NumV=Min(nVec,NumCho(1)-iVec1+1)
         call RdChoVec(Vec2,nPQ,NumV,iVec1,LuCVec(1))

! CASE 1: Integrals are of symmetry type (II/II)
! Coulomb and exchange terms need to be accumulated
! Option code 1: Begin reading at first integral.
      Do J = 1, Size(lists,2)
         I = J
         ID_IP=lists(2,I)

      DO IP_=lists(3,I),lists(4,I)
         IP=hash_table(IP_)
         DO JQ_=lists(3,I),IP_
            JQ=hash_table(JQ_)
! Skip processing (P,Q|... if they do not share the same center
            IPQ=iTri(IP,JQ)
            ! do batches of integrals for a single fixed pair of pq
            CALL Get_Int_DCCD(IRC,X1,IPQ,IJB+1)
            IF(IRC.GT.1) Return
! Do the Coulomb contribution
            IF (nD==1) Then
               TEMP=Zero
               DO KR_=lists(3,I),lists(4,I)
                  KR=hash_table(KR_)
                  DO LS_=lists(3,I),KR_
                     LS=hash_table(LS_)
                     IRS=iTri(KR,LS)
                     TEMP=TEMP+X1(IRS)*DLT(IRS,1)
                  END DO
               END DO
               FLT(IPQ,1)=FLT(IPQ,1)+TEMP
            ELSE
               TEMP=Zero
               TEMP_ab=Zero
               DO KR_=lists(3,I),lists(4,I)
                  KR=hash_table(KR_)
                  DO LS_=lists(3,I),KR_
                     LS=hash_table(LS_)
                     IRS=iTri(KR,LS)
                     TEMP=TEMP+X1(IRS)*DLT(IRS,1)
                     TEMP_ab=TEMP_ab+X1(IRS)*DLT(IRS,2)
                  END DO
               END DO
               FLT(IPQ,1)=FLT(IPQ,1)+TEMP+TEMP_ab
               FLT(IPQ,2)=FLT(IPQ,1)
            END IF
#ifdef _DEBUGPRINT_
            write (6,'(a,i5,a,f12.6)') '00 Flt(',IPQ,',1)=',FLT(IPQ,1)
            if(nD==2) then
            write (6,'(a,i5,a,f12.6)') '00 Flt(',IPQ,',2)=',FLT(IPQ,2)
            endif
#endif
! Do the exchange contribution
            CALL SQUARE (X1(:),X2(:),1,IB,IB)
!
            if(nD==1) then
              DO KR_=lists(3,I),lists(4,I)
                 KR=hash_table(KR_)
                 IRQ=(JQ-1)*IB+KR
                 TEMP=Zero
                 DO LS_=lists(3,I),lists(4,I)
                    LS=hash_table(LS_)
                    ISR=(KR-1)*IB+LS
                    ISP=(IP-1)*IB+LS
                    TEMP=TEMP-Factor*ExFac*X2(ISR)*DSQ(ISP,1)
                 END DO
                 FSQ(IRQ,1)=FSQ(IRQ,1)+TEMP
              END DO
            else
              DO KR_=lists(3,I),lists(4,I)
                 KR=hash_table(KR_)
                 IRQ=(JQ-1)*IB+KR
                 TEMP=Zero
                 TEMP_ab=Zero
                 DO LS_=lists(3,I),lists(4,I)
                    LS=hash_table(LS_)
                    ISR=(KR-1)*IB+LS
                    ISP=(IP-1)*IB+LS
                    TEMP=TEMP-Factor*ExFac*X2(ISR)*DSQ(ISP,1)
                    TEMP_ab=TEMP_ab-Factor*ExFac*X2(ISR)*DSQ(ISP,2)
                 END DO
                 FSQ(IRQ,1)=FSQ(IRQ,1)+TEMP
                 FSQ(IRQ,2)=FSQ(IRQ,2)+TEMP_ab
              END DO
            endif
            IF ( IP.NE.JQ ) THEN
               ISF=(IP-1)*IB+1
               ISD=(JQ-1)*IB+1
!
            if(nD==1) then
              DO KR_=lists(3,I),lists(4,I)
                 KR=hash_table(KR_)
                 IRP=(IP-1)*IB+KR
                 TEMP=Zero
                 DO LS_=lists(3,I),lists(4,I)
                    LS=hash_table(LS_)
                    ISR=(KR-1)*IB+LS
                    ISQ=(JQ-1)*IB+LS
                    TEMP=TEMP-Factor*ExFac*X2(ISR)*DSQ(ISQ,1)
                 END DO
                 FSQ(IRP,1)=FSQ(IRP,1)+TEMP
              END DO
            else
              DO KR_=lists(3,I),lists(4,I)
                 KR=hash_table(KR_)
                 IRP=(IP-1)*IB+KR
                 TEMP=Zero
                 TEMP_ab=Zero
                 DO LS_=lists(3,I),lists(4,I)
                    LS=hash_table(LS_)
                    ISR=(KR-1)*IB+LS
                    ISQ=(JQ-1)*IB+LS
                    TEMP=TEMP-Factor*ExFac*X2(ISR)*DSQ(ISQ,1)
                    TEMP_ab=TEMP_ab-Factor*ExFac*X2(ISR)*DSQ(ISQ,2)
                 END DO
                 FSQ(IRP,1)=FSQ(IRP,1)+TEMP
                 FSQ(IRP,2)=FSQ(IRP,2)+TEMP_ab
              END DO
            endif

            ENDIF
#ifdef _DEBUGPRINT_
          ISF=(JQ-1)*IB
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv,',1)=',FSQ(ISF+ivv,1),ivv=1,Ib)
          if(nD==2) then
          write (6,'(a,i5,a,f12.6)')
     &          ('01 Fsq(',isf+ivv,',2)=',FSQ(ISF+ivv,2),ivv=1,Ib)
          endif
#endif

         END DO  ! JQ
      END DO    ! IP
      END DO    ! I
      end do  ! end of the batch procedure


      call mma_deallocate(lists)
      Call Get_Int_Close()
      call mma_deallocate(hash_table)
      Call mma_deallocate(Basis_IDs)

      END SUBROUTINE FOCKTWO_scf_DCCD

      END SUBROUTINE FOCKTWO_scf
