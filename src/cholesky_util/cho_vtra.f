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
      SUBROUTINE CHO_VTRA(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,IREDC,
     &                   iSwap,nDen,kDen,MOs,ChoT,iSkip)

*********************************************************
*   Author: F. Aquilante
*
*   Purpose:  SCR(lscr) contains JNUM cholesky vectors
*             starting from JVEC1 and stored in reduced
*             sets. The routine performs an MOs half transformation
*             of these elements in a set of target
*             arrays, ChoT, identified the data type Laq_Type.
*             In the target arrays, the vectors are
*             stored in full dimension and as a
*             subset of a a given NUMV number of vectors.
*
*   Input:
*       jVref =  index of the first vector to be transformed
*                computed wrt to the first vector in the
*                target arrays (see calling routine)
*
*       Jvec1 =  first vector to be transformed
*       JNum  =  # of vectors to be transformed
*
*       NumV  =  total # of vectors in the target arrays
*
*       nDen  =  total # of densities to which MOs refer
*       kDen  =  first density to be treated
*
*       iSwap :   = 0   L(k,b,J) is returned
*                 = 1   L(a,k,J) is returned
*                 = 2   L(k,J,b) is returned
*                 = 3   L(a,J,k) is returned
*
*       iSkip(syma)=0 : skip the symmetry block a.
*                    Any vector L(ak) with syma x symk=JSYM
*                    won't be returned in the target array
*
*       IREDC :  reduced set in core at the moment of
*                the first call to the routine.
*                Can be set to -1 by the calling routine
*
*********************************************************
      use ChoArr, only: nDimRS, iRS2F
      use ChoSwp, only: InfVec, IndRed
      use Data_Structures, only: CMO_Type, Laq_Type
      use Data_Structures, only: Map_to_Laq
      Implicit Real*8 (a-h,o-z)

      Integer irc, nDen,kDen,lScr
      Real*8  Scr(lscr)
      Type (CMO_Type) MOs(nDen)
      Type (Laq_Type) ChoT(nDen)

      Integer iSkip(*)

      Integer, External:: cho_isao

      Character(Len=8), Parameter:: SECNAM = 'CHO_VTRA'

      Real*8::  Fac(0:1)=[0.5D0,1.0D0]

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Integer, Allocatable:: nPorb(:,:)
      Integer, Allocatable:: ipChoT(:,:)

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
************************************************************************

      Call mma_allocate(nPorb,8,nDen,Label='nPorb')
      Do iDen = 1, nDen
        Do iSym = 1, nSym
          nPorb(iSym,iDen)=SIZE(MOs(iDen)%pA(iSym)%A,1)
        End Do
      End Do
!     Temporry code while cleaning up.
      Call mma_allocate(ipChoT,8,nDen,Label='ipChoT')
      Do iDen = 1, nDen
         Call Map_to_Laq(ChoT(iDen),ipChoT(:,iDen))
      End Do

**********************************************************
C
C    From Reduced sets to half-MOs full storage
C    ------------------------------------------
C
C    iSwap = 0
C
C     L{a,b,J} ---> L(p,b,J)   ! stride-1 transformation
C
C    iSwap = 1
C
C     L{a,b,J} ---> L(a,q,J)
C
C    iSwap = 2
C
C     L{a,b,J} ---> L(p,J,b)   ! stride-1 transformation
C
C    iSwap = 3
C
C     L{a,b,J} ---> L(a,J,q)
**********************************************************

      iLoc = 3 ! use scratch location in reduced index arrays

      IF (iSwap.eq.0) THEN     ! L(pb,J) storage

         NREAD = 0
         kchot = 0
         DO JVEC=1,JNUM   ! Relative index in the JNUM batch

            LVEC = JVEC - 1 + JVREF  ! Relative index in the NUMV batch
            kVEC = JVEC - 1 + JVEC1  ! Absolute index
            JRED = InfVec(KVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b presnt

               Do jRab=1,nnBstR(jSym,iLoc)

                  kRab = iiBstr(jSym,iLoc) + jRab
                  iRab = IndRed(kRab,iLoc)

                  iag   = iRS2F(1,iRab)  !global address
                  ibg   = iRS2F(2,iRab)

                  iSyma = cho_isao(iag) !symmetry block; S(b)=S(a)=S(p)

                  kscr  = kscr + 1

                  IF (iSkip(iSyma).ne.0) THEN

                     ias   = iag - ibas(iSyma) !address within that sym
                     ibs   = ibg - ibas(iSyma)
                     xfd   = Fac(min(abs(ias-ibs),1)) !fac for diag
C
C     L(p,b,J) = sum_a  xfd* L(a,b,J) * C(p,a)
C     ----------------------------------------
                     DO jDen=kDen,nDen

                        ichot = nPorb(iSyma,jDen)*nBas(iSyma)*(LVEC-1)
     &                        + ipChoT(iSyma,jDen)

*                       kchot = ichot + nPorb(iSyma,jDen)*(ias-1)

                        ! C(1,b)
                        CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                              MOs(JDen)%pA(iSyma)%A(:,ibs),1,
*    &                              Work(kchot),1)
     &                             ChoT(jDen)%pA(iSyma)%A(:,ias,LVEC),1)

*                       kchot = ichot + nPorb(iSyma,jDen)*(ibs-1)

                        ! C(1,a)
                        CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                              MOs(JDen)%pA(iSyma)%A(:,ias),1,
*    &                              Work(kchot),1)
     &                             ChoT(jDen)%pA(iSyma)%A(:,ibs,LVEC),1)

                     END DO  ! loop over densities

                  ENDIF  ! skipping blocks check

               End Do  ! jRab loop

            ELSE  ! jSym.ne.1

               Do jRab=1,nnBstR(jSym,iLoc)

                  kRab = iiBstr(jSym,iLoc) + jRab
                  iRab = IndRed(kRab,iLoc)

                  iag   = iRS2F(1,iRab)  !global address
                  ibg   = iRS2F(2,iRab)

                  iSyma = cho_isao(iag)  !symmetry block
                  iSymb = muld2h(jSym,iSyma)

                  kscr  = kscr + 1

                  ias   = iag - ibas(iSyma)  !address within that sym
                  ibs   = ibg - ibas(iSymb)

                  IF (iSkip(iSyma).ne.0) THEN
C
C     L(p,b,J) = sum_a  L(a,b,J) * C(p,a)
C     -----------------------------------
                     DO jDen=kDen,nDen

*                       kchot = nPorb(iSyma,jDen)*nBas(iSymb)*(LVEC-1)
*    &                        + nPorb(iSyma,jDen)*(ibs-1)
*    &                        + ipChoT(iSyma,jDen)

                        CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),
     &                              MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                              Work(kchot),1)
     &                             ChoT(jDen)%pA(iSyma)%A(:,ibs,LVEC),1)

                     END DO

                  ENDIF  ! skipping block

                  IF (iSkip(iSymb).ne.0) THEN
C
C     L(p,a,J) = sum_b  L(a,b,J) * C(p,b)
C     -----------------------------------
                     DO jDen=kDen,nDen

*                       kchot = nPorb(iSymb,jDen)*nBas(isyma)*(LVEC-1)
*    &                        + nPorb(iSymb,jDen)*(ias-1)
*    &                        + ipChoT(iSymb,jDen)

                        CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),
     &                              MOs(jDen)%pA(iSymb)%A(:,ibs),1,
*    &                              Work(kchot),1)
     &                             ChoT(jDen)%pA(iSymb)%A(:,ias,LVEC),1)

                     END DO

                  ENDIF  ! skipping blocks check

               End Do  ! jRab loop

            ENDIF ! total symmetric vectors check

         END DO

      ELSEIF (iSwap.eq.1) THEN     ! L(ap,J) storage

      NREAD = 0
      kchot = 0
      DO JVEC=1,JNUM

         LVEC = JVEC - 1 + JVREF
         kVEC = JVEC - 1 + JVEC1  ! Absolute index
         JRED = InfVec(KVEC,2,jSym)

         IF (JRED .NE. IREDC) THEN ! JRED is not the reduced set in core
            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            IREDC = JRED
         END IF

         kscr = NREAD
         NREAD= NREAD + nDimRS(jSym,JRED)

        IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b present

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)=Sym(p)

            kscr  = kscr + 1

            IF (iSkip(iSyma).ne.0) THEN

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            xfd   = Fac(min(abs(ias-ibs),1)) !fac for diagonal elements
C
C     L(a,p,J) = sum_b  xfd* L(a,b,J) * C(p,b)
C     ----------------------------------------
            DO jDen=kDen,nDen

*              kchot = nBas(iSyma)*nPorb(iSyma,jDen)*(LVEC-1)
*    &               + ipChoT(iSyma,jDen) - 1

               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                     MOs(jDen)%pA(iSyma)%A(:,ibs),1,
*    &                     Work(kchot+ias),nBas(iSyma))
     &                ChoT(jDen)%pA(iSyma)%A(ias,1,LVEC),nBas(iSyma))

               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                     MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                     Work(kchot+ibs),nBas(iSyma))
     &                ChoT(jDen)%pA(iSyma)%A(ibs,1,LVEC),nBas(iSyma))

            END DO  ! loop over densities

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ELSE  ! jSym.ne.1

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = muld2h(jSym,iSyma)

            kscr  = kscr + 1

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSymb)

            IF (iSkip(iSyma).ne.0) THEN
C
C     L(a,q,J) = sum_b  L(a,b,J) * C(q,b)
C     -----------------------------------
             DO jDen=kDen,nDen

*              kchot = nBas(iSyma)*nPorb(iSymb,jDen)*(LVEC-1)
*    &               + ias
*    &               + ipChoT(iSyma,jDen) - 1

               CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),
     &                     MOs(jDen)%pA(iSymb)%A(:,ibs),1,
*    &                     Work(kchot),nBas(iSyma))
     &                ChoT(jDen)%pA(iSyma)%A(ias,1,LVEC),nBas(iSyma))

             END DO

            ENDIF  ! skipping block

            IF (iSkip(iSymb).ne.0) THEN
C
C     L(b,q,J) = sum_a  L(a,b,J) * C(q,a)
C     -----------------------------------
             DO jDen=kDen,nDen

*              kchot = nBas(iSymb)*nPorb(iSyma,jDen)*(LVEC-1)
*    &               + ibs
*    &               + ipChoT(iSymb,jDen) - 1

               CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),
     &                     MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                     Work(kchot),nBas(iSyma))
     &                ChoT(jDen)%pA(iSymb)%A(ibs,1,LVEC),nBas(iSyma))

             END DO

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ENDIF ! total symmetric vectors check

      END DO


      ELSEIF (iSwap.eq.2) THEN     ! L(pJ,b) storage

             NREAD = 0
             kchot = 0
             DO JVEC=1,JNUM

                LVEC = JVEC - 1 + JVREF
                kVEC = JVEC - 1 + JVEC1  ! Absolute index
                JRED = InfVec(KVEC,2,jSym)

                IF (JRED .NE. IREDC) THEN
                   Call Cho_X_SetRed(irc,iLoc,JRED)
                   IREDC = JRED
                END IF

             kscr = NREAD
             NREAD= NREAD + nDimRS(jSym,JRED)

             IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b

                Do jRab=1,nnBstR(jSym,iLoc)

                   kRab = iiBstr(jSym,iLoc) + jRab
                   iRab = IndRed(kRab,iLoc)

                   iag   = iRS2F(1,iRab)  !global address
                   ibg   = iRS2F(2,iRab)

                   iSyma = cho_isao(iag)

                   kscr  = kscr + 1

                   IF (iSkip(iSyma).ne.0) THEN

                      ias   = iag - ibas(iSyma)
                      ibs   = ibg - ibas(iSyma)
                      xfd   = Fac(min(abs(ias-ibs),1))
C
C     L(p,J,b) = sum_a  xfd* L(a,b,J) * C(p,a)
C     ----------------------------------------
                      DO jDen=kDen,nDen

*                        ichot = nPorb(iSyma,jDen)*(LVEC-1)
*    &                         + ipChoT(iSyma,jDen)

*                        kchot = ichot + nPorb(iSyma,jDen)*NUMV*(ias-1)

                         CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                               MOs(jDen)%pA(iSyma)%A(:,ibs),1,
*    &                               Work(kchot),1)
     &                             ChoT(jDen)%pA(iSyma)%A(:,LVEC,ias),1)

*                        kchot = ichot + nPorb(iSyma,jDen)*NUMV*(ibs-1)

                         CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                               MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                               Work(kchot),1)
     &                             ChoT(jDen)%pA(iSyma)%A(:,LVEC,ibs),1)

                      END DO  ! loop over densities

                   ENDIF  ! skipping blocks check

                End Do  ! jRab loop

             ELSE  ! jSym.ne.1

                Do jRab=1,nnBstR(jSym,iLoc)

                   kRab = iiBstr(jSym,iLoc) + jRab
                   iRab = IndRed(kRab,iLoc)

                   iag   = iRS2F(1,iRab)  !global address
                   ibg   = iRS2F(2,iRab)

                   iSyma = cho_isao(iag)  !symmetry block
                   iSymb = muld2h(jSym,iSyma)  !(syma>symb)

                   kscr  = kscr + 1

                   ias   = iag - ibas(iSyma)
                   ibs   = ibg - ibas(iSymb)

                   IF (iSkip(iSyma).ne.0) THEN
C
C     L(p,J,b) = sum_a  L(a,b,J) * C(p,a)
C     -----------------------------------
                      DO jDen=kDen,nDen

*                        kchot = nPorb(iSyma,jDen)*NUMV*(ibs-1)
*    &                         + nPorb(iSyma,jDen)*(LVEC-1)
*    &                         + ipChoT(iSyma,jDen)

                         CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),
     &                               MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                               Work(kchot),1)
     &                             ChoT(jDen)%pA(iSyma)%A(:,LVEC,ibs),1)

                      END DO

                   ENDIF  ! skipping block

                   IF (iSkip(iSymb).ne.0) THEN
C
C     L(p,J,a) = sum_b  L(a,b,J) * C(p,b)
C     -----------------------------------
                      DO jDen=kDen,nDen


*                        kchot = nPorb(iSymb,jDen)*NUMV*(ias-1)
*    &                         + nPorb(iSymb,jDen)*(LVEC-1)
*    &                         + ipChoT(iSymb,jDen)

                         CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),
     &                               MOs(jDen)%pA(iSymb)%A(:,ibs),1,
*    &                               Work(kchot),1)
     &                             ChoT(jDen)%pA(iSymb)%A(:,LVEC,ias),1)

                      END DO

                   ENDIF  ! skipping blocks check

                End Do  ! jRab loop

              ENDIF ! total symmetric vectors check

           END DO

      ELSEIF (iSwap.eq.3) THEN     ! L(aJ,p) storage

      NREAD = 0
      kchot = 0
      DO JVEC=1,JNUM

         LVEC = JVEC - 1 + JVREF
         KVEC = JVEC - 1 + JVEC1  ! Absolute index
         JRED = InfVec(KVEC,2,jSym)

         IF (JRED .NE. IREDC) THEN ! JRED is not the reduced set in core
            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            IREDC = JRED
         END IF

         kscr = NREAD
         NREAD= NREAD + nDimRS(jSym,JRED)

        IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b present

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)=Sym(p)

            kscr  = kscr + 1

            IF (iSkip(iSyma).ne.0) THEN

            ias   = iag - ibas(iSyma)  !address within that sym block
            ibs   = ibg - ibas(iSyma)
            xfd   = Fac(min(abs(ias-ibs),1)) !scale fac for diag
C
C     L(a,J,p) = sum_b  xfd* L(a,b,J) * C(p,b)
C     ----------------------------------------
            DO jDen=kDen,nDen

*              kchot = nBas(iSyma)*NUMV
*    &               + nBas(iSyma)*(LVEC-1)
*    &               + ipChoT(iSyma,jDen) - 1

               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                     MOs(jDen)%pA(iSyma)%A(:,ibs),1,
*    &                     Work(kchot+ias),nBas(iSyma)*NUMV)
     &                     ChoT(jDen)%pA(iSyma)%A(ias,LVEC,1),
     &                     nBas(iSyma)*NUMV)

               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
     &                     MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                     Work(kchot+ibs),nBas(iSyma)*NUMV)
     &                     ChoT(jDen)%pA(iSyma)%A(ibs,LVEC,1),
     &                     nBas(iSyma)*NUMV)

            END DO  ! loop over densities

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ELSE  ! jSym.ne.1

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = muld2h(jSym,iSyma)

            kscr  = kscr + 1

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSymb)

            IF (iSkip(iSyma).ne.0) THEN
C
C     L(a,J,q) = sum_b  L(a,b,J) * C(q,b)
C     -----------------------------------
             DO jDen=kDen,nDen

*              kchot = nBas(iSyma)*NUMV
*    &               + nBas(iSyma)*(LVEC-1) + ias
*    &               + ipChoT(iSyma,jDen) - 1

               CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),
     &                     MOs(jDen)%pA(iSymb)%A(:,ibs),1,
*    &                     Work(kchot),nBas(iSyma)*NUMV)
     &                     ChoT(jDen)%pA(iSyma)%A(ias,LVEC,1),
     &                     nBas(iSyma)*NUMV)

             END DO

            ENDIF  ! skipping block

            IF (iSkip(iSymb).ne.0) THEN
C
C     L(b,J,q) = sum_a  L(a,b,J) * C(q,a)
C     -----------------------------------
             DO jDen=kDen,nDen

*              kchot = nBas(iSymb)*NUMV
*    &               + nBas(iSymb)*(LVEC-1) + ibs
*    &               + ipChoT(iSymb,jDen) - 1

               CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),
     &                     MOs(jDen)%pA(iSyma)%A(:,ias),1,
*    &                     Work(kchot),nBas(iSyma)*NUMV)
     &                     ChoT(jDen)%pA(iSymb)%A(ibs,LVEC,1),
     &                     nBas(iSyma)*NUMV)

             END DO

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ENDIF ! total symmetric vectors check

      END DO


      ELSE   ! iSwap check

         write(6,*)SECNAM//': invalid argument. Iswap= ',Iswap
         irc = 66
         Return

      ENDIF  ! iSwap check

      Call mma_deallocate(ipChoT)
      Call mma_deallocate(nPorb)
      irc=0

      Return
      END

**************************************************************
