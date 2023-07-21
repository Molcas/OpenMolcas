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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_DecChk_Int(irc,lUnit,Col,Nai,Nbj,ibj1,NumVec,Work,lWork,Fac)
!
! Thomas Bondo Pedersen, Jan. 2005.
!
! Purpose: compute consecutive columns of (ai|bj) matrix from
!          vectors on file (unit: lUnit)
!          vectors from MP2 decomposition).

implicit real*8(a-h,o-z)
real*8 Col(Nai,Nbj), Work(lWork)

irc = 0

! Check dimensions.
! -----------------

if ((Nai < 1) .or. (Nbj < 1) .or. (Nbj > Nai)) then
  irc = -1
  return
end if
ibj2 = ibj1+Nbj-1
if ((ibj1 < 1) .or. (ibj2 > Nai)) then
  irc = -2
  return
end if

! Scale result array.
! -------------------

call dScal_(Nai*Nbj,Fac,Col,1)
if (NumVec < 1) return

! Set up batch.
! -------------
nVec = min(lWork/Nai,NumVec)
if (nVec < 1) then
  irc = 1
  return
end if
nBat = (NumVec-1)/nVec+1

! Start batch loop.
! -----------------

do iBat=1,nBat

  ! Set batch info.
  ! ---------------

  if (iBat == nBat) then
    NumV = NumVec-nVec*(nBat-1)
  else
    NumV = nVec
  end if
  iVec1 = nVec*(iBat-1)+1

  ! Read vectors.
  ! -------------

  iOpt = 2
  lTot = Nai*NumV
  iAdr = Nai*(iVec1-1)+1
  call ddaFile(lUnit,iOpt,Work,lTot,iAdr)

  ! Compute integrals.
  ! ------------------

  call DGEMM_('N','T',Nai,Nbj,NumV,1.0d0,Work,Nai,Work(ibj1),Nai,1.0d0,Col,Nai)

end do

end subroutine ChoMP2_DecChk_Int
