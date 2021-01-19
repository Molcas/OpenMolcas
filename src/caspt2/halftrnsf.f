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
*               2007, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE HALFTRNSF(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,
     &                     JREDC,CMO,ISTART,NUSE,ipChoT)
*********************************************************
*   Author: F. Aquilante as subroutine cho_vtra
*   Modified PAM 2007: Use ordinary CMO array without restructuring
*   This amounts to (1) having symmetry blocks of MO:s accessed
*   untransposed, (2) not through pointers to workspace array Work()
*   but giving the index of start orbital in each symmetry instead,
*   and (3) having CMO array in call parameter list.
*
*   Purpose:  SCR(lscr) contains JNUM cholesky vectors
*             starting from JVEC1 and stored in reduced
*             sets. The routine performs an MOs half transformation
*             of these elements in a set of target
*             arrays identified by the pointers ipChoT.
*             In the target arrays, the vectors are
*             stored in full dimension and as a
*             subset of a a given NUMV number of vectors.
*             Each pointer should thereby point to a
*             location where the corresponding Cholesky
*             vector of a given unique symmetry pair
*             of indices has to be stored
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
*       JREDC :  reduced set in core at the moment of
*                the first call to the routine.
*                Can be set to -1 by the calling routine
*
*********************************************************
      use ChoArr, only: iRS2F, nDimRS
      use ChoSwp, only: InfVec, IndRed
      Implicit Real*8 (a-h,o-z)
      Real*8  Scr(lscr)
      Integer ipChoT(8)
      Real*8 CMO(*)
      Integer IOFFC(8),ISTART(8),NUSE(8)

* Integer function cho_isao
      Integer  cho_isao
      External cho_isao

* Uses MAXVEC, NSYM, IIBSTR(8,3), NNBSTR(8,3), NNBSTRT(3) in
* cholesky.fh commons /CHOLEV/, /CHORST/, /CHOSHL/
#include "cholesky.fh"
#include "choptr.fh"
* Uses ibas(8), nbas(8) in choorb.fh common /CHOORB/
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
************************************************************************

* iLoc = 3 means 'use scratch location in reduced index arrays'
      iLoc = 3

* Offset counter into CMO array:
      IOC=0
      DO ISYM=1,NSYM
       IOFFC(ISYM)=IOC
       IOC=IOC+NBAS(ISYM)**2
      END DO

      Do iSymp=1,nSym
         IF (nUse(iSymp).ne.0) THEN
            iSymb = muld2h(JSYM,iSymp)
            nElem=nUse(iSymp)*nBas(iSymb)*NUMV
            Call dCopy_(nElem,[0.0D0],0,Work(ipChoT(iSymp)),1)
         ENDIF
      End Do

      NREAD = 0
      DO JVEC=1,JNUM
* JVTRNS=Cholesky vector to be transformed.
         JVTRNS=JVEC1-1+JVEC
* JRED: Which reduced set does it belong to:
         JRED = InfVec(JVTRNS,2,JSYM)

* Is it the same still?
         IF (JRED .NE. JREDC) THEN
* It is not. Tables must be regenerated with information that is
* common to this reduced set.
        write(6,*)' Rats! It was assumed that the Cholesky vectors'
        write(6,*)' in HALFTRNSF all belonged to a given reduced'
        write(6,*)' set, but they don''t!'
        write(6,*)' JRED, JREDC:',JRED,JREDC
        write(6,*)' Back to the drawing board?'
        write(6,*)' Let the program continue and see what happens.'
            Call Cho_X_SetRed(irc,iLoc,JRED)
            JREDC = JRED
         END IF

         kscr = NREAD
         NREAD = NREAD + nDimRS(JSYM,JRED)

         IF (JSYM.eq.1) THEN
* L(a,b,J)=L(b,a,J); only a.ge.b stored

            Do jRab=1,nnBstR(jSym,iLoc)

               kRab = iiBstr(jSym,iLoc) + jRab
               iRab = IndRed(kRab,3)

* Global address:
               iag = iRS2F(1,iRab)
               ibg = iRS2F(2,iRab)

               iSyma = cho_isao(iag)

               kscr  = kscr + 1

               ISYMB=ISYMA

               NUSEA=nUse(iSyma)
               NUSEB=nUse(iSymb)

               IF (NUSEA.ne.0) THEN

                 ias   = iag - ibas(iSyma)
                 ibs   = ibg - ibas(iSyma)
*  L(p,J,b) = sum( C(a,p)* L(a,b,J), a=1..NBA), where p=1..NUSEA
*  L(p,J,a) = sum( C(b,p)* L(b,a,J), b=1..NBB), where p=1..NUSEB
*  ----------------------------------------

                   NBA=NBAS(ISYMA)
                   ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA)-1)
                   NBB=NBAS(ISYMB)
                   ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB)-1)

                   kchot = ipChoT(iSymb)+NUSEB*
     &                               (jVref+JVEC-2+NUMV*(ias-1))

                   CALL DAXPY_(NUSEA,Scr(kscr),
     &                           CMO(ISCB),NBB,Work(kchot),1)

                 IF(IAS.NE.IBS) THEN
                   kchot = ipChoT(iSyma)+NUSEA*
     &                               (jVref+JVEC-2+NUMV*(ibs-1))

                   CALL DAXPY_(NUSEA,Scr(kscr),
     &                           CMO(ISCA),NBA,Work(kchot),1)
                 END IF


* End of NUSE  test
              ENDIF

* End of loop over basis function pair index JRAB
           End Do

        ELSE
* jSym.ne.1

           Do jRab=1,nnBstR(jSym,iLoc)

              kRab = iiBstr(jSym,iLoc) + jRab
              iRab = IndRed(kRab,3)

* Global address:
              iag = iRS2F(1,iRab)
              ibg = iRS2F(2,iRab)

* iSyma = cho_isao(iag) = symmetry block of basis function iag
              iSyma = cho_isao(iag)
* iSyma > isymb since jsym.ne.1 and a.ge.b
              iSymb = muld2h(jSym,iSyma)
              NUSEA=nUse(iSyma)
              NUSEB=nUse(iSymb)

              kscr  = kscr + 1

              ias   = iag - ibas(iSyma)
              ibs   = ibg - ibas(iSymb)

              IF (NUSEA.ne.0) THEN
*  L(p,J,b) = sum( C(a,p)*L(a,b,J), a=1..NBA), where p=1..NUSEA

                   NBA=NBAS(ISYMA)
                   ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA)-1)

                   kchot = ipChoT(iSyma)+nUse(iSyma)*
     &                               (jVref+JVEC-2+NUMV*(ibs-1))

                   CALL DAXPY_(NUSEA,Scr(kscr),
     &                              CMO(ISCA),NBA,Work(kchot),1)


* End of NUSE  test
              ENDIF

              IF (NUSEB.ne.0) THEN
* L(p,J,a) = sum( C(b,p)*L(b,a,J), b=1..NBA), where p=1..NUSEB

                   NBB=NBAS(ISYMB)
                   ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB)-1)

                   kchot = ipChoT(iSymb)+NUSEB*
     &                               (jVref+JVEC-2+NUMV*(ias-1))

                   CALL DAXPY_(NUSEB,Scr(kscr),
     &                              CMO(ISCB),NBB,Work(kchot),1)


* End of NUSE  test
              END IF

* End of loop over basis function pair index JRAB
           End Do

* End of JSYM.ne.1 test
         ENDIF

* End of JSYM loop
      END DO



      irc=0


      Return
      END
