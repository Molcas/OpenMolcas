!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Picky(mDij,DeDe,nDeDe,nSD,iSD4,i,j)

use Symmetry_Info, only: nIrrep
use Definitions, only: wp, iwp, u6
use Dens_Stuff, only: mDCR12=>mDCRij,mDCR34=>mDCRkl,mDCR13=>mDCRik,mDCR14=>mDCRil,mDCR23=>mDCRjk,mDCR24=>mDCRjl
use Dens_Stuff, only:  ipD12=> ipDij, ipD34=> ipDkl, ipD13=> ipDik, ipD14=> ipDil, ipD23=> ipDjk, ipD24=> ipDjl
use Dens_Stuff, only: ipDD12=>ipDDij,ipDD34=>ipDDkl,ipDD13=>ipDDik,ipDD14=>ipDDil,ipDD23=>ipDDjk,ipDD24=>ipDDjl

implicit none
integer(kind=iwp), intent(in) :: nDeDe, i, j, nSD, iSD4(0:nSD,4)
integer(kind=iwp), intent(out) :: mDij
real(kind=wp), intent(inout) :: DeDe(nDeDe)
integer(kind=iwp) :: ii1, ii2, ii3, jj1, jj2, jj3, i1, i2, i3, j1, j2, j3
integer(kind=iwp) :: iCmpi,jCmpj,iBasi,jBasj,iPrimi,jPrimj,iShell,jShell,iBsInc,jBsInc
integer(kind=iwp) :: iBasAO, jBasAO, iBasn, jBasn
integer(kind=iwp), pointer :: mDCRij=>Null(), ipDij=>Null(), ipDDij=>Null()

iCmpi =iSD4( 2,i)
iBasi =iSD4( 3,i)
iPrimi=iSD4( 5,i)
iShell=iSD4(11,i)
iBsInc=iSD4( 4,i)
iBasAO=iSD4( 8,i)+1
iBasn =iSD4(19,i)

jCmpj =iSD4( 2,j)
jBasj =iSD4( 3,j)
jPrimj=iSD4( 5,j)
jShell=iSD4(11,j)
jBsInc=iSD4( 4,j)
jBasAO=iSD4( 8,j)+1
jBasn =iSD4(19,j)

If (i==1 .and. j==2) Then
   mDCRij=>mDCR12
    ipDij=> ipD12
   ipDDij=>ipDD12
Else If (i==1 .and. j==3) Then
   mDCRij=>mDCR13
    ipDij=> ipD13
   ipDDij=>ipDD13
Else If (i==1 .and. j==4) Then
   mDCRij=>mDCR14
    ipDij=> ipD14
   ipDDij=>ipDD14
Else If (i==2 .and. j==3) Then
   mDCRij=>mDCR23
    ipDij=> ipD23
   ipDDij=>ipDD23
Else If (i==2 .and. j==4) Then
   mDCRij=>mDCR24
    ipDij=> ipD24
   ipDDij=>ipDD24
Else If (i==3 .and. j==4) Then
   mDCRij=>mDCR34
    ipDij=> ipD34
   ipDDij=>ipDD34
Else
   Write (u6,*) 'Picky: illegal i and j combination'
   Write (u6,*) 'i,j=',i,j
   Call Abend()
End If


if (nIrrep == 1) then
  ii1 = 0
  ii2 = 1
  ii3 = 0
  jj1 = 0
  jj2 = 1
  jj3 = 0
else
  ii1 = iBasi
  ii2 = iBasAO
  ii3 = iBasn
  jj1 = jBasj
  jj2 = jBasAO
  jj3 = jBasn
end if
if (mDCRij /= 0) then
  if (iShell >= jShell) then
    i1 = ii1
    i2 = ii2
    i3 = ii3
    j1 = jj1
    j2 = jj2
    j3 = jj3
  else
    i1 = jj1
    i2 = jj2
    i3 = jj3
    j1 = ii1
    j2 = ii2
    j3 = ii3
  end if
  if ((iBasi == iBsInc) .and. (jBasj == jBsInc)) then
    ipDDij = ipDij
  else
    call Picky_inner(DeDe(ipDij),i1,j1,iPrimi*jPrimj,iCmpi*jCmpj,mDCRij,i2,i2+i3-1,j2,j2+j3-1,DeDe(ipDDij))
  end if
end if
mDij = (ii3*jj3+1)*iCmpi*jCmpj+iPrimi*jPrimj+1

return

end subroutine Picky
