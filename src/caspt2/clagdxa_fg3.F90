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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

      Subroutine CLagDXA_FG3(iSym,nAS,nAshT,NG3,NS,BDER,SDER,           &
     &                       DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,             &
     &                       G2,SC,idxG3)

      use Symmetry_Info, only: Mul
      USE SUPERINDEX, only: KTUV
      use caspt2_module, only: IASYM, EPSA, NTUVES
      use Constants, only: Zero
      use definitions, only: wp, iwp, byte

      implicit none

      integer(kind=iwp), intent(in) :: iSym, nAS, nAshT, NG3, NS
      real(kind=wp), intent(in) :: BDER(nAS,nAS), SDER(nAS,nAS),        &
     &  G2(nAshT,nAshT,nAshT,nAshT), SC(NS)
      real(kind=wp), intent(inout) :: DF1(nAshT,nAshT),                 &
     &  DF2(nAshT,nAshT,nAshT,nAshT), DF3(NG3), DG1(nAshT,nAshT),       &
     &  DG2(nAshT,nAshT,nAshT,nAshT), DG3(NG3), DEPSA(nAshT,nAshT)
      integer(kind=byte), intent(in) :: idxG3(6,NG3)

      integer(kind=iwp) :: iG3, iT, iU, iV, iX, iY, iZ, iST, iSU, iSV,  &
     &                     iSX, iSY, iSZ, ituvs, ixyzs, iTU, iVX, iYZ,  &
     &                     jSYM, ISUP, JSUP, iW, NSEQ
      real(kind=wp) :: F3VAL, G3VAL

      DO iG3=1,NG3
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
        iST=IASYM(iT)
        iSU=IASYM(iU)
        iSV=IASYM(iV)
        iSX=IASYM(iX)
        iSY=IASYM(iY)
        iSZ=IASYM(iZ)
        ituvs=Mul(IST,Mul(ISU,ISV))
        ixyzs=Mul(ISX,Mul(ISY,ISZ))
        F3VAL=Zero
        G3VAL=Zero
        if(ituvs /= ixyzs) cycle
        iTU=iT+NASHT*(iU-1)
        iVX=iV+NASHT*(iX-1)
        iYZ=iY+NASHT*(iZ-1)
!-SVC20100829: 12 equivalent cases, of which the second
!  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
!  - G(tuvxyz) -> SA(xut,vyz)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM == iSYM) THEN
          ISUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
        if (iTU /= iVX .or. iVX /= iYZ) then
          if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
!  - G(vxtuyz) -> SA(uxv,tyz)
            jSYM=Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
            IF (jSYM == iSYM) THEN
              ISUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
              JSUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
              F3VAL = F3VAL + BDER(iSup,jSup)
              G3VAL = G3VAL + SDER(iSup,jSup)
            ENDIF
!  - G(yzvxtu) -> SA(xzy,vtu)
            jSYM=Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
            IF (jSYM == iSYM) THEN
              ISUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
              JSUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
              F3VAL = F3VAL + BDER(iSup,jSup)
              G3VAL = G3VAL + SDER(iSup,jSup)
            ENDIF
!  - G(tuyzvx) -> SA(zut,yvx)
            jSYM=Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
            IF (jSYM == iSYM) THEN
              ISUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
              JSUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
              F3VAL = F3VAL + BDER(iSup,jSup)
              G3VAL = G3VAL + SDER(iSup,jSup)
            ENDIF
          end if
!  - G(yztuvx) -> SA(uzy,tvx)
          jSYM=Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
          IF (jSYM == iSYM) THEN
            ISUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            JSUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
            F3VAL = F3VAL + BDER(iSup,jSup)
            G3VAL = G3VAL + SDER(iSup,jSup)
          ENDIF
!  - G(vxyztu) -> SA(zxv,ytu)
          jSYM=Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
          IF (jSYM == iSYM) THEN
            ISUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            JSUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
            F3VAL = F3VAL + BDER(iSup,jSup)
            G3VAL = G3VAL + SDER(iSup,jSup)
          ENDIF
        end if

        if ((iT /= iU .or. iV /= iX .or. iY /= iZ) .and.                &
     &      (iT /= iU .or. iV /= iZ .or. iX /= iY) .and.                &
     &      (iX /= iV .or. iT /= iZ .or. iU /= iY) .and.                &
     &      (iZ /= iY .or. iV /= iU .or. iX /= iT)) then
!  - G(utxvzy) -> SA(vtu,xzy)
          jSYM=Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
          IF (jSYM == iSYM) THEN
            ISUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
            JSUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            F3VAL = F3VAL + BDER(iSup,jSup)
            G3VAL = G3VAL + SDER(iSup,jSup)
          ENDIF
          if (iTU /= iVX .or. iVX /= iYZ) then
            if (iTU /= iVX .and. iTU /= iYZ .and. iVX /= iYZ) then
!  - G(xvutzy) -> SA(tvx,uzy)
              jSYM=Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
              IF (jSYM == iSYM) THEN
                ISUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
                JSUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
                F3VAL = F3VAL + BDER(iSup,jSup)
                G3VAL = G3VAL + SDER(iSup,jSup)
              ENDIF
!  - G(zyxvut) -> SA(vyz,xut)
              jSYM=Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
              IF (jSYM == iSYM) THEN
                ISUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
                JSUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
                F3VAL = F3VAL + BDER(iSup,jSup)
                G3VAL = G3VAL + SDER(iSup,jSup)
              ENDIF
!  - G(utzyxv) -> SA(ytu,zxv)
              jSYM=Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
              IF (jSYM == iSYM) THEN
                ISUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
                JSUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
                F3VAL = F3VAL + BDER(iSup,jSup)
                G3VAL = G3VAL + SDER(iSup,jSup)
              ENDIF
            end if
!  - G(zyutxv) -> SA(tyz,uxv)
            jSYM=Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
            IF (jSYM == iSYM) THEN
              ISUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
              JSUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
              F3VAL = F3VAL + BDER(iSup,jSup)
              G3VAL = G3VAL + SDER(iSup,jSup)
            ENDIF
!  - G(xvzyut) -> SA(yvx,zut)
            jSYM=Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
            IF (jSYM == iSYM) THEN
              ISUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
              JSUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
              F3VAL = F3VAL + BDER(iSup,jSup)
              G3VAL = G3VAL + SDER(iSup,jSup)
            ENDIF
          end if
        end if

        F3VAL = -F3VAL
        G3VAL = -G3VAL

        !! last line of F3 transformation in mkfg3.f
        G3VAL = G3VAL - (EPSA(iU)+EPSA(iY))*F3VAL
        Do iW = 1, nAshT
          ISUP=KTUV(iV,iW,iT)-nTUVES(iSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(iSYM)
          NSEQ=MAX(iSup,jSup)*(MAX(iSup,jSup)-1)/2 + MIN(iSup,jSup)
          DEPSA(iW,iU) = DEPSA(iW,iU) - F3VAL*SC(NSEQ)

          ISUP=KTUV(iV,iU,iT)-nTUVES(iSYM)
          JSUP=KTUV(iX,iW,iZ)-nTUVES(iSYM)
          NSEQ=MAX(iSup,jSup)*(MAX(iSup,jSup)-1)/2 + MIN(iSup,jSup)
          DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*SC(NSEQ)
        End Do

        !! derivative of <0|EtuEwv,xwEyz|0>*fww
        DF3(iG3) = DF3(iG3) + F3VAL
        !! derivative of <0|EtuEvxEyz|0>
        DG3(iG3) = DG3(iG3) + G3VAL

        !! remaining F3 and G3 transformation in mkfg3.f
        if (iY == iX) then
          DF2(iT,iU,iV,iZ) = DF2(iT,iU,iV,iZ) - F3VAL
          DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - EPSA(iU)*F3VAL
          Do iW = 1, nAshT
            DEPSA(iU,iW) = DEPSA(iU,iW) - F3VAL*G2(iT,iW,iV,iZ)
          End Do
          DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - G3VAL
        end if
        if (iV == iU) then
          DF2(iT,iX,iY,iZ) = DF2(iT,iX,iY,iZ) - F3VAL
          DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - EPSA(iY)*F3VAL
          Do iW = 1, nAshT
            DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*G2(iT,iX,iW,iZ)
          End Do
          DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - G3VAL
        end if
        if (iY == iU) then
          DF2(iV,iX,iT,iZ) = DF2(iV,iX,iT,iZ) - F3VAL
          DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - EPSA(iU)*F3VAL
          DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - G3VAL
        end if
        DEPSA(iY,iU) = DEPSA(iY,iU) - F3VAL*G2(iV,iX,iT,iZ)
        if (iY == iX .and. iV == iU) then
          DF1(iT,iZ) = DF1(iT,iZ) - F3VAL
          DG1(iT,iZ) = DG1(iT,iZ) - G3VAL
        end if
      END DO

      Return

      End Subroutine CLagDXA_FG3
