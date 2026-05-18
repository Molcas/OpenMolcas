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

subroutine MKBA_F3(ISYM,BA,MBA,NG3,F3,idxG3)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp, Byte
use SUPERINDEX, only: KTUV
use caspt2_module, only: NASHT, IASYM, nTUVES

implicit none
integer(kind=iwp), intent(in) :: ISYM, MBA, NG3
real(kind=wp), intent(inout) :: BA(MBA)
real(kind=wp), intent(in) :: F3(NG3)
integer(kind=Byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) iG3, iT, iU, iV, iX, iY, iZ, iST, iSU, iSV, iSX, iSY, iSZ, ituvs, ixyzs, iTU, iVX, iYZ, jSYM, ISUP, JSUP, ISADR
real(kind=wp) F3VAL

!-SVC20100831: determine indices in SA where a certain F3 value will end up
do iG3=1,NG3
  iT = idxG3(1,iG3)
  iU = idxG3(2,iG3)
  iV = idxG3(3,iG3)
  iX = idxG3(4,iG3)
  iY = idxG3(5,iG3)
  iZ = idxG3(6,iG3)
  iST = IASYM(iT)
  iSU = IASYM(iU)
  iSV = IASYM(iV)
  iSX = IASYM(iX)
  iSY = IASYM(iY)
  iSZ = IASYM(iZ)
  ituvs = Mul(IST,Mul(ISU,ISV))
  ixyzs = Mul(ISX,Mul(ISY,ISZ))
  if (ituvs /= ixyzs) cycle
  iTU = iT+NASHT*(iU-1)
  iVX = iV+NASHT*(iX-1)
  iYZ = iY+NASHT*(iZ-1)
  F3VAL = -F3(iG3)
  !-SVC20100829: 12 equivalent cases, of which the second
  !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
  !  - F(tuvxyz) -> BA(xut,vyz)
  jSYM = Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iX,iU,iT)-nTUVES(jSYM)
    JSUP = KTUV(iV,iY,iZ)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      BA(ISADR) = BA(ISADR)+F3VAL
    end if
  end if
  if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - F(vxtuyz) -> BA(uxv,tyz)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iU,iX,iV)-nTUVES(jSYM)
        JSUP = KTUV(iT,iY,iZ)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = (ISUP*(ISUP-1))/2+JSUP
          BA(ISADR) = BA(ISADR)+F3VAL
        end if
      end if
      !  - F(yzvxtu) -> BA(xzy,vtu)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iX,iZ,iY)-nTUVES(jSYM)
        JSUP = KTUV(iV,iT,iU)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = (ISUP*(ISUP-1))/2+JSUP
          BA(ISADR) = BA(ISADR)+F3VAL
        end if
      end if
      !  - F(tuyzvx) -> BA(zut,yvx)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iZ,iU,iT)-nTUVES(jSYM)
        JSUP = KTUV(iY,iV,iX)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = (ISUP*(ISUP-1))/2+JSUP
          BA(ISADR) = BA(ISADR)+F3VAL
        end if
      end if
    end if
    !  - F(yztuvx) -> BA(uzy,tvx)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iU,iZ,iY)-nTUVES(jSYM)
      JSUP = KTUV(iT,iV,iX)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        BA(ISADR) = BA(ISADR)+F3VAL
      end if
    end if
    !  - F(vxyztu) -> BA(zxv,ytu)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iZ,iX,iV)-nTUVES(jSYM)
      JSUP = KTUV(iY,iT,iU)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        BA(ISADR) = BA(ISADR)+F3VAL
      end if
    end if
  end if
  if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
  if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
  if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
  if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
  !  - F(utxvzy) -> BA(vtu,xzy)
  jSYM = Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iV,iT,iU)-nTUVES(jSYM)
    JSUP = KTUV(iX,iZ,iY)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      BA(ISADR) = BA(ISADR)+F3VAL
    end if
  end if
  if ((iTU == iVX) .and. (iVX == iYZ)) cycle
  if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
    !  - F(xvutzy) -> BA(tvx,uzy)
    jSYM = Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iT,iV,iX)-nTUVES(jSYM)
      JSUP = KTUV(iU,iZ,iY)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        BA(ISADR) = BA(ISADR)+F3VAL
      end if
    end if
    !  - F(zyxvut) -> BA(vyz,xut)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iV,iY,iZ)-nTUVES(jSYM)
      JSUP = KTUV(iX,iU,iT)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        BA(ISADR) = BA(ISADR)+F3VAL
      end if
    end if
    !  - F(utzyxv) -> BA(ytu,zxv)
    jSYM = Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iY,iT,iU)-nTUVES(jSYM)
      JSUP = KTUV(iZ,iX,iV)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = (ISUP*(ISUP-1))/2+JSUP
        BA(ISADR) = BA(ISADR)+F3VAL
      end if
    end if
  end if
  !  - F(zyutxv) -> BA(tyz,uxv)
  jSYM = Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iT,iY,iZ)-nTUVES(jSYM)
    JSUP = KTUV(iU,iX,iV)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      BA(ISADR) = BA(ISADR)+F3VAL
    end if
  end if
  !  - F(xvzyut) -> BA(yvx,zut)
  jSYM = Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iY,iV,iX)-nTUVES(jSYM)
    JSUP = KTUV(iZ,iU,iT)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = (ISUP*(ISUP-1))/2+JSUP
      BA(ISADR) = BA(ISADR)+F3VAL
    end if
  end if
end do

end subroutine MKBA_F3
