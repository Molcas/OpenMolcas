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
*               2021, Roland Lindh                                     *
************************************************************************

      SUBROUTINE CHO_GetShFull(LabJ,lLabJ,JNUM,JSYM,IREDC,ChoV,
     &                         SvShp,mmShl,iShp_rs,mmShl_tot)
      use ChoArr, only: iSOShl, iShlSO, iBasSh, iRS2F, nDimRS
      use ChoSwp, only: iiBstRSh, IndRSh, IndRed
      use Data_Structures, only: L_Full_Type
      Implicit Real*8 (a-h,o-z)
      Real*8  LabJ(lLabJ)
      Real*8  SvShp(mmShl , 2)
      Integer iShp_rs( mmShl_tot )
      Integer, External :: cho_isao

#include "cholesky.fh"
#include "choorb.fh"
#include "real.fh"

      Type (L_Full_Type) ChoV

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
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

      ChoV%A0(:)=Zero
      SvShp(:,:)=Zero

      IF (JSYM.eq.1) THEN

         NREAD = 0

         DO JVEC=1,JNUM

            kLabJ = NREAD
            NREAD= NREAD + nDimRS(jSym,IREDC)

            Do jRab=1,nnBstR(jSym,iLoc)

             kRab = iiBstr(jSym,iLoc) + jRab
             iRab = IndRed(kRab,iLoc) ! addr in 1st red set

             iShp = IndRSh(iRab) ! shell pair to which it belongs

             iag  = iRS2F(1,iRab)  !global address
             ibg  = iRS2F(2,iRab)

             iaSh = iSOShl(iag) ! shell to which it belongs
             ibSh = iSOShl(ibg) ! iaSh >= ibSh !!!!!!

             iaSg = iShlSO(iag) !index of SO within its shell
             ibSg = iShlSO(ibg)

             iSyma = cho_isao(iag)  !symmetry block sym(a)=sym(b)

             ias = iaSg - iBasSh(iSyma,iaSh) ! addr within its shell
             ibs = ibSg - iBasSh(isyma,ibSh)

             kLabJ  = kLabJ + 1

             ChoV%SPB(iSyma,iShp_rs(iShp),1)%A3(ias,JVEC,ibs)
     &        = LabJ(kLabJ)

             i1=1
             If (ibSh/=iaSh) i1=2

             ChoV%SPB(iSyma,iShp_rs(iShp),i1)%A3(ibs,JVEC,ias)
     &        = LabJ(kLabJ)

             SvShp(iShp_rs(iShp),2) = SvShp(iShp_rs(iShp),2)
     &                                  + LabJ(kLabJ)**2

            End Do

            Do jShp=1,nnShl_tot ! Maximize over vectors
               If (iShp_rs(jShp).gt.0) Then
                  SvShp(iShp_rs(jShp),1) = Max( SvShp(iShp_rs(jShp),1),
     &                                          SvShp(iShp_rs(jShp),2) )
                  SvShp(iShp_rs(jShp),2) = zero
               End If
            End Do

         END DO


      ELSE

         NREAD = 0

         DO JVEC=1,JNUM

            kLabJ = NREAD
            NREAD= NREAD + nDimRS(jSym,IREDC)

            Do jRab=1,nnBstR(jSym,iLoc)

              kRab = iiBstr(jSym,iLoc) + jRab
              iRab = IndRed(kRab,iLoc) ! addr in 1st red set

              iShp = IndRSh(iRab) ! shell pair to which it belongs

              iag  = iRS2F(1,iRab)  !global address
              ibg  = iRS2F(2,iRab)

              iaSh = iSOShl(iag) ! shell to which it belongs
              ibSh = iSOShl(ibg) ! ibsh<=>iaSh

              iaSg = iShlSO(iag) !index of SO within its shell
              ibSg = iShlSO(ibg)

              iSyma = cho_isao(iag)  !symmetry block
              iSymb = muld2h(jSym,iSyma) ! iSyma >= iSymb

              i1=1
              If (iaSh<ibSh) i1=2

              ias = iaSg - iBasSh(iSyma,iaSh) !addr within its shell
              ibs = ibSg - iBasSh(iSymb,ibSh)

              kLabJ  = kLabJ + 1

              ChoV%SPB(iSyma,iShp_rs(iShp),i1)%A3(ias,JVEC,ibs)
     &        = LabJ(kLabJ)

              SvShp(iShp_rs(iShp),2) = SvShp(iShp_rs(iShp),2)
     &                                   + LabJ(kLabJ)**2

            End Do

            Do jShp=1,nnShl_tot
               If (iShp_rs(jShp).gt.0) Then
                  SvShp(iShp_rs(jShp),1) = Max( SvShp(iShp_rs(jShp),1),
     &                                          SvShp(iShp_rs(jShp),2) )
                  SvShp(iShp_rs(jShp),2) = zero
               EndIf
            End Do

         END DO


      ENDIF


      Return
      END

**************************************************************
