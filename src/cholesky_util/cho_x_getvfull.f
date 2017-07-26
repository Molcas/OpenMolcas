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
*  Cho_X_getVfull
*
*> @brief
*>   Reorder Cholesky vectors from reduced to full storage
*> @author F. Aquilante
*>
*> @details
*> This routine reorders Cholesky vectors from reduced to
*> full storage. For \p DoRead = ``.true.`` the vectors are read from
*> disk using array \p RedVec as scratch space, whereas for
*> \p DoRead = ``.false.`` the reduced vectors must be supplied in
*> array \p RedVec.
*>
*> Given a set of pointers (\p ipChoV) to target
*> arrays, the routine reorders \p NUMV Cholesky
*> vectors of compound symmetry \p ISYM starting with
*> vector \p JVEC1 and returns them in the target arrays.
*> Each pointer should thereby point to a
*> location where the corresponding Cholesky
*> vector of a given unique symmetry pair
*> of indices has to be stored.
*>
*> - \p iSwap = ``0``: \f$ L(a,b,J) \f$ is returned (in LT-storage if \f$ \mathrm{sym}(a)=\mathrm{sym}(b) \f$)
*> - \p iSwap = ``1``: \f$ L(a,J,b) \f$ is returned
*> - \p iSwap = ``2``: \f$ L(a,J,b) \f$ is returned (in SQ-storage if \f$ \mathrm{sym}(a)=\mathrm{sym}(b) \f$)
*>
*> - \p iSkip(syma) = ``0``: skip the symmetry block \f$ a \f$. Any vector \f$ L_{ab} \f$ or \f$ L_{ba} \f$
*>                           with \c syma &times; \c symb = \p ISYM won't be returned in the target array
*>
*> - \p IREDC: reduced set in core at the moment of the call to the routine.
*>             Can be set to ``-1`` (= unknown or undefined) by the calling routine.
*>
*> @param[out] irc     return code
*> @param[in]  RedVec  vectors stored in reduced set(s) [\p DoRead option off]
*>                     or scratch space for reading reduced vectors [\p DoRead option on]
*> @param[in]  lRedVec size of the \p RedVec
*> @param[in]  IVEC1   first vector to read
*> @param[in]  NUMV    number of vectors to read starting from \p IVEC1
*> @param[in]  ISYM    compound symmetry of the Cholesky vectors
*> @param[in]  iSwap   type of the full storage for the returned Cholesky vectors
*> @param[in]  IREDC   current reduced set in core (location ``3``)
*> @param[in]  ipChoV  pointers to the target arrays
*> @param[in]  iSkip   skipping parameters for each symmetry block \f$ (ab) \f$ of compound symmetry \p ISYM
*> @param[in]  DoRead  flag for reading reduced vectors from disk
************************************************************************
      Subroutine Cho_X_getVfull(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                          iSwap,IREDC,ipChoV,iSkip,DoRead)
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
