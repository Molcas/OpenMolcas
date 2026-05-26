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

subroutine MKBC_F3(ISYM,BC,NBC,NG3,F3,idxG3)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use SUPERINDEX, only: KTUV
use caspt2_module, only: IASYM, NASHT, nTUVES
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: ISYM, NBC, NG3
real(kind=wp), intent(inout) :: BC(NBC)
real(kind=wp), intent(in) :: F3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: iG3, ISADR, iST, iSU, ISUP, iSV, iSX, iSY, iSZ, iT, iTU, ituvs, iU, iV, iVX, iX, ixyzs, iY, iYZ, iZ, JSUP, jSYM
real(kind=wp) :: F3VAL

!-SVC20100831: determine indices in BC where a certain G3 value will end up
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
  F3VAL = F3(iG3)
  !-SVC20100829: 12 equivalent cases, of which the second
  !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
  !  - F(tuvxyz) -> BC(vut,xyz)
  jSYM = Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iV,iU,iT)-nTUVES(jSYM)
    JSUP = KTUV(iX,iY,iZ)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = iTri(ISUP,JSUP)
      BC(ISADR) = BC(ISADR)+F3VAL
    end if
  end if
  if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - F(vxtuyz) -> BC(txv,uyz)
      jSYM = Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iT,iX,iV)-nTUVES(jSYM)
        JSUP = KTUV(iU,iY,iZ)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = iTri(ISUP,JSUP)
          BC(ISADR) = BC(ISADR)+F3VAL
        end if
      end if
      !  - F(yzvxtu) -> BC(vzy,xtu)
      jSYM = Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iV,iZ,iY)-nTUVES(jSYM)
        JSUP = KTUV(iX,iT,iU)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = iTri(ISUP,JSUP)
          BC(ISADR) = BC(ISADR)+F3VAL
        end if
      end if
      !  - F(tuyzvx) -> BC(yut,zvx)
      jSYM = Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
      if (jSYM == iSYM) then
        ISUP = KTUV(iY,iU,iT)-nTUVES(jSYM)
        JSUP = KTUV(iZ,iV,iX)-nTUVES(jSYM)
        if (JSUP <= ISUP) then
          ISADR = iTri(ISUP,JSUP)
          BC(ISADR) = BC(ISADR)+F3VAL
        end if
      end if
    end if
    !  - F(yztuvx) -> BC(tzy,uvx)
    jSYM = Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iT,iZ,iY)-nTUVES(jSYM)
      JSUP = KTUV(iU,iV,iX)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = iTri(ISUP,JSUP)
        BC(ISADR) = BC(ISADR)+F3VAL
      end if
    end if
    !  - F(vxyztu) -> BC(yxv,ztu)
    jSYM = Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iY,iX,iV)-nTUVES(jSYM)
      JSUP = KTUV(iZ,iT,iU)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = iTri(ISUP,JSUP)
        BC(ISADR) = BC(ISADR)+F3VAL
      end if
    end if
  end if
  if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
  if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
  if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
  if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
  !  - F(utxvzy) -> BC(xtu,vzy)
  jSYM = Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iX,iT,iU)-nTUVES(jSYM)
    JSUP = KTUV(iV,iZ,iY)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = iTri(ISUP,JSUP)
      BC(ISADR) = BC(ISADR)+F3VAL
    end if
  end if
  if ((iTU == iVX) .and. (iVX == iYZ)) cycle
  if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
    !  - F(xvutzy) -> BC(uvx,tzy)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iU,iV,iX)-nTUVES(jSYM)
      JSUP = KTUV(iT,iZ,iY)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = iTri(ISUP,JSUP)
        BC(ISADR) = BC(ISADR)+F3VAL
      end if
    end if
    !  - F(zyxvut) -> BC(xyz,vut)
    jSYM = Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iX,iY,iZ)-nTUVES(jSYM)
      JSUP = KTUV(iV,iU,iT)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = iTri(ISUP,JSUP)
        BC(ISADR) = BC(ISADR)+F3VAL
      end if
    end if
    !  - F(utzyxv) -> BC(ztu,yxv)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      ISUP = KTUV(iZ,iT,iU)-nTUVES(jSYM)
      JSUP = KTUV(iY,iX,iV)-nTUVES(jSYM)
      if (JSUP <= ISUP) then
        ISADR = iTri(ISUP,JSUP)
        BC(ISADR) = BC(ISADR)+F3VAL
      end if
    end if
  end if
  !  - F(zyutxv) -> BC(uyz,txv)
  jSYM = Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iU,iY,iZ)-nTUVES(jSYM)
    JSUP = KTUV(iT,iX,iV)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = iTri(ISUP,JSUP)
      BC(ISADR) = BC(ISADR)+F3VAL
    end if
  end if
  !  - F(xvzyut) -> BC(zvx,yut)
  jSYM = Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
  if (jSYM == iSYM) then
    ISUP = KTUV(iZ,iV,iX)-nTUVES(jSYM)
    JSUP = KTUV(iY,iU,iT)-nTUVES(jSYM)
    if (JSUP <= ISUP) then
      ISADR = iTri(ISUP,JSUP)
      BC(ISADR) = BC(ISADR)+F3VAL
    end if
  end if
end do

end subroutine MKBC_F3
