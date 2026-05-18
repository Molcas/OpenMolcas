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
      SUBROUTINE MKSC_G3(ISYM,SC,NSC,NG3,G3,idxG3)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp, Byte
      USE SUPERINDEX, only: KTUV
      use caspt2_module, only: NASHT, IASYM, NTUVES
      IMPLICIT None

      integer(kind=iwp), intent(in):: ISYM,NSC,NG3
      real(kind=wp), intent(out):: SC(NSC)
      real(kind=wp), intent(in):: G3(NG3)
      INTEGER(kind=Byte), intent(in):: idxG3(6,NG3)

      integer(kind=iwp) iG3,iT,iU,iV,iX,iY,iZ,iST,iSU,iSV,iSX,iSY,iSZ,  &
     &                  ituvs,ixyzs,iTU,iVX,iYZ,JSYM,ISUP,JSUP,ISADR
      real(kind=wp) G3VAL

!-SVC20100831: determine indices in SC where a certain G3 value will end up
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
        G3VAL=G3(iG3)
!-SVC20100829: 12 equivalent cases, of which the second
!  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
!  - G(tuvxyz) -> SC(vut,xyz)
        jSYM=Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
        if (.NOT.(iTU.eq.iVX.and.iVX.eq.iYZ)) THEN
        if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - G(vxtuyz) -> SC(txv,uyz)
        jSYM=Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(yzvxtu) -> SC(vzy,xtu)
        jSYM=Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(tuyzvx) -> SC(yut,zvx)
        jSYM=Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
       ENDIF
!  - G(yztuvx) -> SC(tzy,uvx)
        jSYM=Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(vxyztu) -> SC(yxv,ztu)
        jSYM=Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
       ENDIF
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
!  - G(utxvzy) -> SC(xtu,vzy)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
        if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - G(xvutzy) -> SC(uvx,tzy)
        jSYM=Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(zyxvut) -> SC(xyz,vut)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(utzyxv) -> SC(ztu,yxv)
        jSYM=Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
       ENDIF
!  - G(zyutxv) -> SC(uyz,txv)
        jSYM=Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(xvzyut) -> SC(zvx,yut)
        jSYM=Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
      END DO

      END SUBROUTINE MKSC_G3
