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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      SUBROUTINE MKBC_F3(ISYM,BC,NBC,NG3,F3,idxG3)
      use Symmetry_Info, only: Mul
      USE SUPERINDEX, only: KTUV
      use caspt2_module, only: NASHT,IASYM,nTUVES
      use definitions, only: iwp, wp, Byte

      IMPLICIT NONE

      INTEGER(KIND=IWP), INTENT(IN):: ISYM, NBC,NG3
      REAL(KIND=WP), INTENT(INOUT):: BC(NBC)
      REAL(KIND=WP), INTENT(IN):: F3(NG3)
      INTEGER(KIND=BYTE), INTENT(IN):: idxG3(6,NG3)

      INTEGER(KIND=IWP) iG3,iT,iU,iV,iX,iY,iZ,iST,iSU,iSV,iSX,iSY,iSZ,  &
     &                  ituvs,ixyzs,iTU,iVX,iYZ,jSYM,ISUP,JSUP,ISADR
      REAL(KIND=WP) F3VAL

!-SVC20100831: determine indices in BC where a certain G3 value will end up
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
        if(ituvs.ne.ixyzs) CYCLE
        iTU=iT+NASHT*(iU-1)
        iVX=iV+NASHT*(iX-1)
        iYZ=iY+NASHT*(iZ-1)
        F3VAL=F3(iG3)
!-SVC20100829: 12 equivalent cases, of which the second
!  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
!  - F(tuvxyz) -> BC(vut,xyz)
        jSYM=Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
        if (.NOT.(iTU.eq.iVX.and.iVX.eq.iYZ)) THEN
        if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - F(vxtuyz) -> BC(txv,uyz)
        jSYM=Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
!  - F(yzvxtu) -> BC(vzy,xtu)
        jSYM=Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
!  - F(tuyzvx) -> BC(yut,zvx)
        jSYM=Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
       ENDIF
!  - F(yztuvx) -> BC(tzy,uvx)
        jSYM=Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
!  - F(vxyztu) -> BC(yxv,ztu)
        jSYM=Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
       ENDIF
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
!  - F(utxvzy) -> BC(xtu,vzy)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
        if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - F(xvutzy) -> BC(uvx,tzy)
        jSYM=Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
!  - F(zyxvut) -> BC(xyz,vut)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
!  - F(utzyxv) -> BC(ztu,yxv)
        jSYM=Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
       ENDIF
!  - F(zyutxv) -> BC(uyz,txv)
        jSYM=Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
!  - F(xvzyut) -> BC(zvx,yut)
        jSYM=Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
      END DO

      END SUBROUTINE MKBC_F3
