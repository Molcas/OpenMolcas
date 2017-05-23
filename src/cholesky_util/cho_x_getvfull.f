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
      Subroutine Cho_X_getVfull(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                          iSwap,IREDC,ipChoV,iSkip,DoRead)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_getVfull</Name>
*     <Syntax>Call Cho\_X\_getVfull(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,iSwap,IREDC,ipChoV,iSkip,DoRead)</Syntax>
*     <Arguments>
*       \Argument{irc}{return code}{Integer}{out}
*       \Argument{RedVec}{Vectors stored in reduced set(s) [DoRead option off] or scratch space for reading reduced vectors [DoRead option on]}{Real*8}{in}
*       \Argument{lRedVec}{size of the RedVec}{Integer}{in}
*       \Argument{IVEC1}{first vector to read}{Integer}{in}
*       \Argument{NUMV}{number of vectors to read starting from IVEC1}{Integer}{in}
*       \Argument{ISYM}{compund symmetry of the Cholesky vectors}{Integer}{in}
*       \Argument{iSwap}{type of the full storage for the returned Cholesky vectors}{Integer}{in}
*       \Argument{IREDC}{current reduced set in core (location 3)}{Integer}{in}
*       \Argument{ipChoV}{pointers to the target arrays}{Integer}{in}
*       \Argument{iSkip}{skipping parameters for each symmetry block
*       (ab) of compound symmety ISYM.}{Integer}{in}
*       \Argument{DoRead}{flag for reading reduced vectors from disk}{Logical}{in}
*     </Arguments>
*     <Purpose>
*             This routine reorders Cholesky vectors from reduced to
*             full storage. For DoRead=.true. the vectors are read from
*             disk using array RedVec as scratch space, whereas for
*             DoRead=.false. the reduced vectors must be supplied in
*             array RedVec.
*             Given a set of pointers (ipChoV) to target
*             arrays, the routine reorders NUMV Cholesky
*             vectors of compound symmetry ISYM starting with
*             vector JVEC1 and returns them in the target arrays.
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
*       iSwap :   = 0   L(a,b,J) is returned
*                       (in LT-storage if sym(a)=sym(b))
*                 = 1   L(a,J,b) is returned
*                 = 2   L(a,J,b) is returned
*                       (in SQ-storage if sym(a)=sym(b))
*
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

********************************************************
*   Author: F. Aquilante
*
*   Purpose:  This routine reorders Cholesky vectors from reduced to
*             full storage. For DoRead=.true. the vectors are read from
*             disk using array RedVec as scratch space, whereas for
*             DoRead=.false. the reduced vectors must be supplied in
*             array RedVec.
*             Given a set of pointers (ipChoV) to target
*             arrays, the routine reorders NUMV Cholesky
*             vectors of compound symmetry ISYM starting with
*             vector JVEC1 and returns them in the target arrays.
*
*   Input:
*       Ivec1 =  first vector to be reordered
*       NumV  =  # of vectors to be reordered
*
*       iSwap :   = 0   L(a,b,J) is returned
*                       (in LT-storage if sym(a)=sym(b))
*                 = 1   L(a,J,b) is returned
*                 = 2   L(a,J,b) is returned
*                       (in SQ-storage if sym(a)=sym(b))
*
*       iSkip(syma)=0 : skip the symmetry block a.
*                    Any vector L(ab) or L(ba) with syma x symb=ISYM
*                    won't be returned in the target array
*
*       IREDC :  reduced set in core at the moment of
*                the call to the routine.
*                Can be set to -1 by the calling routine
*
********************************************************
      Implicit Real*8 (a-h,o-z)
      Dimension RedVec(lRedVec)
      Integer   ipVec(8),nnBSF(8,8),n2BSF(8,8)
      Integer   ipChoV(*),iSkip(*)
      Logical   DoRead
      Character*14 SECNAM
      Parameter (SECNAM = 'Cho_X_GetVfull')

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      MXUSD = 0
      MUSED = 0

      do i=1,nSym
         ipVec(i) = ipChoV(i)
      end do


      Call set_nnBSF(nSym,Nbas,nnBSF,n2BSF)


      if(iSwap.eq.0)then

       Do iSymq=1,nSym
         iSymp = muld2h(iSym,iSymq)
         if(nnBSF(iSymp,iSymq).gt.0)then
          if (iSymp.ge.iSymq .and. iSkip(iSymp).ne.0 ) then
            Call FZero(Work(ipChoV(iSymp)),nnBSF(iSymp,iSymq)*NUMV)
          endif
         endif
       End Do

      elseif(iSwap.eq.1 .or. iSwap.eq.2)then

       Do iSymq=1,nSym
         iSymp = muld2h(iSym,iSymq)
         if(n2BSF(iSymp,iSymq).gt.0)then
          if (iSymp.ge.iSymq .and. iSkip(iSymp).ne.0 ) then
            Call FZero(Work(ipChoV(iSymp)),n2BSF(iSymp,iSymq)*NUMV)
          endif
         endif
       End Do

      else

         write(6,*)'Wrong parameter! iSwap= ',iSwap
         irc = 66
         Return

      endif


      If (DoRead) Then

       JVEC1 = IVEC1
       IVEC2 = JVEC1 + NUMV - 1

       Do While (jVec1.le.iVec2)
        Call CHO_VECRD(RedVec,lRedVec,JVEC1,IVEC2,ISYM,
     &                   JNUM,IREDC,MUSED)
        MXUSD = MAX(MXUSD,MUSED)

          If (JNUM.le.0 .or. JNUM.gt.(IVEC2-JVEC1+1)) then
             irc=77
             RETURN
          End If

        jVref = JVEC1 - IVEC1 + 1
        Call cho_Reordr(irc,RedVec,lRedVec,jVref,JVEC1,JNUM,NUMV,ISYM,
     &                 IREDC,iSwap,ipVec,iSkip)

        if (irc.ne.0) then
           return
        endif

        jVec1 = jVec1 + JNUM

        do i=1,nSym
           j=mulD2h(i,ISYM)
           IF ( j.ge.i .and. iSkip(j).ne.0 ) THEN
              if(iSwap.eq.0)then
                ipVec(j) = ipVec(j) + nnBSF(j,i)*JNUM
              elseif(iSwap.eq.1)then
                ipVec(j) = ipChoV(j)
              elseif(iSwap.eq.2)then
                ipVec(j) = ipVec(j) + n2BSF(j,i)*JNUM
              endif
           ENDIF
         end do

       End Do  ! end the while loop

      Else ! only reorder

       JNUM=NUMV

       Call cho_Reordr(irc,RedVec,lRedVec,1,IVEC1,JNUM,NUMV,ISYM,IREDC,
     &                iSwap,ipVec,iSkip)

        if (irc.ne.0) then
           return
        endif

      End If

      irc=0

      RETURN
      END
