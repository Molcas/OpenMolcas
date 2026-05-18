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
      SUBROUTINE MKSA_G3(ISYM,SA,NSA,NG3,G3,idxG3)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp, Byte
      USE SUPERINDEX, only: KTUV
      use caspt2_module, only: NASHT, IASYM, NTUVES
      IMPLICIT None

      INTEGER(kind=iwp), intent(in):: ISYM,NSA,NG3
      real(kind=wp), intent(out):: SA(NSA)
      real(kind=wp), intent(in):: G3(NG3)
      INTEGER(kind=Byte), intent(in):: idxG3(6,NG3)

      integer(kind=iwp) iG3,iT,iU,iV,iX,iY,iZ,iST,iSU,iSV,iSX,iSY,iSZ,  &
     &                  ituvs,ixyzs,iTU,iVX,iYZ,JSYM,ISUP,JSUP,ISADR
      real(kind=wp) G3VAL

!-SVC20100831: determine indices in SA where a certain G3 value will end up
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
        G3VAL=-G3(iG3)
!-SVC20100829: 12 equivalent cases, of which the second
!  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
!  - G(tuvxyz) -> SA(xut,vyz)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF

        if (.NOT.(iTU.eq.iVX.and.iVX.eq.iYZ)) THEN

        if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - G(vxtuyz) -> SA(uxv,tyz)
        jSYM=Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(yzvxtu) -> SA(xzy,vtu)
        jSYM=Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(tuyzvx) -> SA(zut,yvx)
        jSYM=Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
       ENDIF

!  - G(yztuvx) -> SA(uzy,tvx)
        jSYM=Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(vxyztu) -> SA(zxv,ytu)
        jSYM=Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF

       ENDIF

        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
!  - G(utxvzy) -> SA(vtu,xzy)
        jSYM=Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE

        if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - G(xvutzy) -> SA(tvx,uzy)
        jSYM=Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(zyxvut) -> SA(vyz,xut)
        jSYM=Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(utzyxv) -> SA(ytu,zxv)
        jSYM=Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
       ENDIF

!  - G(zyutxv) -> SA(tyz,uxv)
        jSYM=Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
!  - G(xvzyut) -> SA(yvx,zut)
        jSYM=Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF

      END DO

      END SUBROUTINE MKSA_G3
