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
      Subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                         iSwap,IREDC,nDen,kDen,ipMOs,nPorb,ipChoT,
     &                         iSkip,DoRead)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_getVtra</Name>
*     <Syntax>Call Cho\_X\_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,iSwap,IREDC,nDen,ipMOs,nPorb,ipChoT,iSkip,DoRead)</Syntax>
*     <Arguments>
*       \Argument{irc}{return code}{Integer}{out}
*       \Argument{RedVec}{Vectors stored in reduced set(s) [DoRead
*       option off] or scratch space for reading reduced vectors
*       [DoRead option on]}{Real*8}{in}
*       \Argument{lRedVec}{size of the RedVec}{Integer}{in}
*       \Argument{IVEC1}{first vector to read}{Integer}{in}
*       \Argument{NUMV}{number of vectors to transform starting from IVEC1}{Integer}{in}
*       \Argument{ISYM}{compound symmetry of the Cholesky vectors}{Integer}{in}
*       \Argument{iSwap}{type of the full storage for the half transformed Cholesky vectors}{Integer}{in}
*       \Argument{IREDC}{reduced set in core}{Integer}{in}
*       \Argument{nDen}{total number of densities to which MOs refer}{Integer}{in}
*       \Argument{kDen}{first density for which the MO transformation has to be performed}{Integer}{in}
*       \Argument{ipMOs}{matrix (8 x nDen) of pointers to the MOs coefficients}{Integer}{in}
*       \Argument{nPorb}{number of orbitals in the primary space for a
*       given symmetry and density}{Integer}{in}
*       \Argument{ipChoT}{pointers to the half transformed vectors}{Integer}{in}
*       \Argument{iSkip}{skipping parameters for each symmetry block
*       (ab) of compound symmety ISYM.}{Integer}{in}
*       \Argument{DoRead}{flag for reading the reduced vectors}{Logical}{in}
*     </Arguments>
*     <Purpose>
*             This routine performs a half-MO-transformation of Cholesky vectors stored in reduced
*             storage. For DoRead=.true. the vectors are read from
*             disk using array RedVec as scratch space, whereas for
*             DoRead=.false. the reduced vectors must be supplied in
*             array RedVec.
*             Given a set of pointers (ipChoT) to target
*             arrays, the routine performs a half-MO-transformation of NUMV Cholesky
*             vectors of compound symmetry ISYM starting with
*             vector IVEC1 and returns them in the target arrays.
*             Each pointer should thereby point to a
*             location where the corresponding Cholesky
*             vector of a given unique symmetry pair
*             of indices has to be stored
*         </Purpose>
*     <Dependencies></Dependencies>
*     <Author>F. Aquilante</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*
*       iSwap :   = 0   L(k,b,J) is returned
*                 = 1   L(a,k,J) is returned
*                 = 2   L(k,J,b) is returned
*                 = 3   L(a,J,k) is returned

*
*       iSkip(syma)=0 : skip the symmetry block a.
*                    Any vector L(ab) or L(ba) with syma x symb=ISYM
*                    won't be returned in the target array
*
*       IREDC :  reduced set in core at the moment of
*                the call to the routine.
*                Can be set to -1 (= unknown or undefined)
*                by the calling routine.
*     </Description>
*    </DOC>
*
************************************************************

      Implicit Real*8 (a-h,o-z)
      Dimension RedVec(lRedVec)
      Integer   ipChoT(8,*),nDen,kDen
      Integer   iSkip(*),ipMOs(8,*),nPorb(8,*)
      Logical   DoRead
      Character*13 SECNAM
      Parameter (SECNAM = 'Cho_X_GetVtra')

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      ipVec(i,j) = iWork(ip_Vec-1+8*(j-1)+i)
**************************************************

      MXUSD = 0
      MUSED = 0

C zeroing the target arrays
C--------------------------
      Do iSymp=1,nSym
         IF (iSkip(iSymp).ne.0) THEN
            iSymb = muld2h(ISYM,iSymp)
            Do jDen=kDen,nDen
               If (iSwap.eq.0 .or. iSwap.eq.2) then ! Lpb,J or LpJ,b
                  Call FZero(Work(ipChoT(iSymp,jDen)),
     &                 nPorb(iSymp,jDen)*nBas(iSymb)*NUMV)
               ElseIf (iSwap.eq.1 .or. iSwap.eq.3) then !Laq,J or LaJ,q
                  Call FZero(Work(ipChoT(iSymp,jDen)),
     &                 nPorb(iSymb,jDen)*nBas(iSymp)*NUMV)
               EndIf
            End Do
         ENDIF
      End Do

C --- define local pointers to the target arrays
C ----------------------------------------------
      Call GetMem('ip_Vec','Allo','Inte',ip_Vec,8*nDen)

      Do jDen=kDen,nDen
         Do i=1,nSym
            iWork(ip_Vec-1+8*(jDen-1)+i) = ipChoT(i,jDen)
         End Do
      End do


C ===============================================
      IF (DoRead) THEN

       JVEC1 = IVEC1
       IVEC2 = JVEC1 + NUMV - 1

       Do While (jVec1.le.iVec2)

        Call CHO_VECRD(RedVec,lRedVec,JVEC1,IVEC2,ISYM,
     &                 JNUM,IREDC,MUSED)

        MXUSD = MAX(MXUSD,MUSED)

          If (JNUM.le.0 .or. JNUM.gt.(IVEC2-JVEC1+1)) then
             irc=77
             RETURN
          End If

        jVref = JVEC1 - IVEC1 + 1

        Call cho_vTra(irc,RedVec,lRedVec,jVref,JVEC1,JNUM,NUMV,ISYM,
     &            IREDC,iSwap,nDen,kDen,ipMOs,nPorb,iWork(ip_Vec),iSkip)

        if (irc.ne.0) then
           return
        endif

        jVec1 = jVec1 + JNUM

C --- Updating the local pointers to the target arrays (iff iSwap=0,1)
C --------------------------------------------------------------------
        IF(iSwap.eq.0)THEN    ! Lpb,J

         Do jDen=kDen,nDen
          do iSymp=1,nSym
             iSymb=mulD2h(iSymp,ISYM)
             if(iSkip(iSymp).ne.0)then
               iWork(ip_Vec-1+8*(jDen-1)+iSymp) = ipVec(iSymp,jDen)
     &                            + nPorb(iSymp,jDen)*nBas(iSymb)*JNUM
             endif
          end do
         End Do

        ELSEIF(iSwap.eq.1)THEN    ! Laq,J

         Do jDen=kDen,nDen
          do iSyma=1,nSym
             iSymq=mulD2h(iSyma,ISYM)
             if(iSkip(iSyma).ne.0)then
               iWork(ip_Vec-1+8*(jDen-1)+iSyma) = ipVec(iSyma,jDen)
     &                           + nBas(iSyma)*nPorb(iSymq,jDen)*JNUM
             endif
          end do
         End Do

        ELSEIF(iSwap.ne.2 .or. iSwap.ne.3)THEN

         write(6,*)SECNAM//': invalid argument. Iswap= ',Iswap
         irc=66
         return

        ENDIF

       End Do  ! end the while loop


      ELSE ! only MO transformation


       JNUM = NUMV

       Call cho_vTra(irc,RedVec,lRedVec,1,IVEC1,JNUM,NUMV,ISYM,IREDC,
     &              iSwap,nDen,kDen,ipMOs,nPorb,iWork(ip_Vec),iSkip)

        if (irc.ne.0) then
           return
        endif


      END IF

      Call GetMem('ip_Vec','Free','Inte',ip_Vec,8*nDen)

      irc=0

      RETURN
      END
