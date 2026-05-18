!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Process_RHS_Block(ITI,ITP,ITK,ITQ,                     &
     &                             Case,                                &
     &                             Cho_Bra,nBra,Cho_Ket,nKet,           &
     &                             nSh,JSYM,IVEC,NV)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use caspt2_global, only: iPrGlb, PIQK, BUFF, idxb
      use PrintLevel, only: DEBUG
      use caspt2_module, only: NSYM
      use AddRHS, only: ADDRHSA, ADDRHSB, ADDRHSC, ADDRHSD1, ADDRHSD2,  &
     &                  ADDRHSE, ADDRHSF, ADDRHSG, ADDRHSH

      IMPLICIT None
      integer(kind=iwp), Intent(in):: ITI,ITP,ITK,ITQ
      Character(LEN=2), intent(in)::  Case
      integer(kind=iwp), intent(in):: nBra, nKet
      real(kind=wp), intent(in):: Cho_Bra(nBra), Cho_Ket(nKet)
      integer(kind=iwp), intent(in):: nSh(8,3), JSYM, iVec, nV

      integer(kind=iwp) ISYI, ISYK, ISYP, ISYQ, KPI, KQK, LBRASM,       &
     &                  LKETSM, NBRASM, NI, NK, NKETSM, NP, NPI, NPIQK, &
     &                  NQ, NQK
      integer(kind=iwp) mxPIQK, nBuff
      mxPIQK=size(PIQK)
      nBuff=size(BUFF)
!
!
      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'Processing RHS block '//Case
      END IF

      LBRASM=1
      DO ISYI=1,NSYM
         NI=NSH(ISYI,ITI)
         IF(NI.EQ.0) Cycle
         ISYP=Mul(ISYI,JSYM)
         NP=NSH(ISYP,ITP)
         IF(NP.EQ.0) Cycle
         NPI=NP*NI
         NBRASM=NPI*NV
!
         LKETSM=1
         DO ISYK=1,NSYM
            NK=NSH(ISYK,ITK)
            IF(NK.EQ.0) Cycle
            ISYQ=Mul(ISYK,JSYM)
            NQ=NSH(ISYQ,ITQ)
            IF(NQ.EQ.0) Cycle
            NQK=NQ*NK
            NKETSM=NQK*NV
!
! SVC: we need an NPI*NQK to store the 2-electron integrals, and 2
! buffers (values+indices) for sorting them.  Later, we can try to get
! rid of the buffer that stores the values and only use an index buffer
! and the two-electron integrals for the scatter operation.  For the
! buffer, any size can be taken, but assuming there is enough memory
! available, it's set to the size of the two-electron integrals unless
! larger than some predefined maximum buffer size.
            NPIQK=NPI*NQK
            IF (NPIQK.GT.MXPIQK) THEN
              IF (Case.eq.'H') THEN
                KPI=MXPIQK/NQK
                NPIQK=KPI*NQK
              ELSE IF (Case.eq.'G') THEN
                KQK=MXPIQK/NPI
                NPIQK=NPI*KQK
              ELSE
                WRITE(6,*) ' NPIQK > MXPIQK and case != G or H'
                WRITE(6,'(A,A2)')  ' CASE =   ', Case
                WRITE(6,'(A,I12)') ' NPIQK =  ', NPIQK
                WRITE(6,'(A,I12)') ' MXPIQK = ', MXPIQK
                WRITE(6,*) ' This should not happen, please report.'
                CALL AbEnd()
              END IF
            END IF

!-SVC: sanity check
            IF (NPIQK.LE.0) THEN
              WRITE(6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
              CALL AbEnd()
            END IF
!
            If (Case.eq.'A ') Then
               CALL ADDRHSA(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,                           &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'B ') Then
               CALL ADDRHSB(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,                           &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'D1') Then
               CALL ADDRHSD1(IVEC,JSYM,ISYI,ISYK,                       &
     &                       NP,NI,NQ,NK,PIQK,                          &
     &                       nBuff,Buff,idxb,                           &
     &                       Cho_Bra(LBRASM),                           &
     &                       Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'H ') Then
               CALL ADDRHSH(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,NPIQK,                     &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'C ') Then
               CALL ADDRHSC(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,                           &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'F ') Then
               CALL ADDRHSF(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,                           &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'D2') Then
               CALL ADDRHSD2(IVEC,JSYM,ISYI,ISYK,                       &
     &                       NP,NI,NQ,NK,PIQK,                          &
     &                       nBuff,Buff,idxb,                           &
     &                       Cho_Bra(LBRASM),                           &
     &                       Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'G ') Then
               CALL ADDRHSG(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,NPIQK,                     &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'E ') Then
               CALL ADDRHSE(IVEC,JSYM,ISYI,ISYK,                        &
     &                      NP,NI,NQ,NK,PIQK,                           &
     &                      nBuff,Buff,idxb,                            &
     &                      Cho_Bra(LBRASM),                            &
     &                      Cho_Ket(LKETSM),NV)
            Else
               Call Abend()
            End If
!
         LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
      END DO
!
      End Subroutine Process_RHS_Block
