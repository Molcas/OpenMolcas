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

      SUBROUTINE CHO_GetShFull(scr,lscr,JNUM,JSYM,
     &                         IREDC,ipChoV,SvShp,iShp_rs)
************************************************************
*   Author: F. Aquilante
*
*
************************************************************
      use ChoArr, only: iSOShl, iShlSO, iBasSh, nBasSh, iRS2F, nDimRS
      use ChoSwp, only: iiBstRSh, IndRSh
      Implicit Real*8 (a-h,o-z)
      Real*8  Scr(lscr),SvShp(*)
      Integer iShp_rs(*)
      Integer cho_isao
      External cho_isao

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Integer   ipChoV
      Parameter (zero = 0.0d0)

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)
****** this is a trick to save memory. Memory in "location 2" is used
******      to store some offset arrays
      iOffShp(i,j) = iiBstRSh(i,j,2)
************************************************************************

**********************************************************
C
C    From Reduced sets to full storage
C    ---------------------------------
C
C     L{a,b,J} ---> L(a,J,b)
C
**********************************************************

      iLoc = 3 ! use scratch location in reduced index arrays


      IF (JSYM.eq.1) THEN

         NREAD = 0

         DO JVEC=1,JNUM

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,IREDC)

            Do jRab=1,nnBstR(jSym,iLoc)

             kRab = iiBstr(jSym,iLoc) + jRab
             iRab = IndRed(kRab,iLoc) ! addr in 1st red set

             iShp = IndRSh(iRab) ! shell pair to which it belongs

             iag  = iRS2F(1,iRab)  !global address
             ibg  = iRS2F(2,iRab)

             iaSh = iSOShl(iag) ! shell to which it belongs
             ibSh = iSOShl(ibg) ! iaSh >= ibSh

             iaSg = iShlSO(iag) !index of SO within its shell
             ibSg = iShlSO(ibg)

             iSyma = cho_isao(iag)  !symmetry block sym(a)=sym(b)

             ias = iaSg - iBasSh(iSyma,iaSh) ! addr within its shell
             ibs = ibSg - iBasSh(isyma,ibSh)

             kscr  = kscr + 1

             kcho1 = nBasSh(iSyma,iaSh)*JNUM*(ibs-1)
     &             + nBasSh(iSyma,iaSh)*(JVEC-1) + ias

             kcho2 = nBasSh(iSyma,ibSh)*JNUM*(ias-1)
     &             + nBasSh(iSyma,ibSh)*(JVEC-1) + ibs

             ioffV = ipChoV + JNUM*iOffShp(iSyma,iShp_rs(iShp)) - 1

             Work(kcho1+ioffV) = Scr(kscr)

             ioffV = ioffV - nBasSh(iSyma,iaSh)*JNUM*nBasSh(iSyma,ibSh)
     &             *Min(0,(ibSh-iaSh))/Max(1,(iaSh-ibSh))

             Work(kcho2+ioffV) = Scr(kscr)

             SvShp(nnShl+iShp_rs(iShp)) = SvShp(nnShl+iShp_rs(iShp))
     &                                  + Scr(kscr)**2

            End Do

            Do jShp=1,nnShl_tot ! Maximize over vectors
               If (iShp_rs(jShp).gt.0) Then
                  SvShp(iShp_rs(jShp)) = Max( SvShp(iShp_rs(jShp)),
     &                                   SvShp(nnShl+iShp_rs(jShp)) )
                  SvShp(nnShl+iShp_rs(jShp)) = zero
               End If
            End Do

         END DO


      ELSE


         NREAD = 0

         DO JVEC=1,JNUM

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,IREDC)

            Do jRab=1,nnBstR(jSym,iLoc)

              kRab = iiBstr(jSym,iLoc) + jRab
              iRab = IndRed(kRab,iLoc) ! addr in 1st red set

              iShp = IndRSh(iRab) ! shell pair to which it belongs

              iag  = iRS2F(1,iRab)  !global address
              ibg  = iRS2F(2,iRab)

              iaSh = iSOShl(iag) ! shell to which it belongs
              ibSh = iSOShl(ibg)

              iaSg = iShlSO(iag) !index of SO within its shell
              ibSg = iShlSO(ibg)

              iSyma = cho_isao(iag)  !symmetry block
              iSymb = muld2h(jSym,iSyma) ! iSyma >= iSymb

c               koff = iOffShp(iSyma,iShp_rs(iShp))

c               If (iaSh.lt.ibSh)   koff = koff
c     &                           + nBasSh(iSyma,ibSh)*nBasSh(iSymb,iaSh)

              koff = iOffShp(iSyma,iShp_rs(iShp)) - nBasSh(iSyma,ibSh)*
     &               nBasSh(iSymb,iaSh)*
     &               Min(0,(iaSh-ibSh))/Max(1,(ibSh-iaSh))

              ias = iaSg - iBasSh(iSyma,iaSh) !addr within its shell
              ibs = ibSg - iBasSh(iSymb,ibSh)

              kchov = nBasSh(iSyma,iaSh)*JNUM*(ibs-1)
     &              + nBasSh(iSyma,iaSh)*(JVEC-1) + ias
     &              + ipChoV + JNUM*kOff - 1

              kscr  = kscr + 1

              Work(kchov) = Scr(kscr)

              SvShp(nnShl+iShp_rs(iShp)) = SvShp(nnShl+iShp_rs(iShp))
     &                                   + Scr(kscr)**2

            End Do

            Do jShp=1,nnShl_tot
               If (iShp_rs(jShp).gt.0) Then
                  SvShp(iShp_rs(jShp)) = Max( SvShp(iShp_rs(jShp)),
     &                                      SvShp(nnShl+iShp_rs(jShp)) )
                  SvShp(nnShl+iShp_rs(jShp)) = zero
               EndIf
            End Do

         END DO


      ENDIF


      Return
      END

**************************************************************
