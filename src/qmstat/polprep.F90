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

subroutine PolPrep(Dist,DistIm,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,iCNum,nSize)

use qmstat_global, only: CordIm, Cordst, nCent, nPart, nPol
use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCNum, nSize
real(kind=wp), intent(in) :: Dist(nCent,nCent,nTri_Elem(nPart-iCNum-1)), DistIm(nCent,nPart-iCNum,nCent,nPart-iCNum)
real(kind=wp), intent(out) :: xx(nSize,nSize), yy(nSize,nSize), zz(nSize,nSize), rr3(nSize,nSize), xxi(nSize,nSize), &
                              yyi(nSize,nSize), zzi(nSize,nSize), Gri(nSize,nSize)
integer(kind=iwp) :: i, i1, ild, imd, Indco1, IndCo2, Indp1, IndP2, IndTr, IndTr1, IndTri, j, j1, jnd, k, l, ncParm

!----------------------------------------------------------------------*
! Simply compute some vectors etc. for the ensuing polarization        *
! calculation.                                                         *
!----------------------------------------------------------------------*
ncParm = nCent*nPart-nCent*icNum
! Loop over solvent polarization sites.
rr3(nPol*iCNum+1:,nPol*iCNum+1:) = Zero
IndTr1 = 0
! If (iCNum+1) == nPart this loop will not be run, but that is
! okey since then we only have one solvent and thus it can not
! experience any field from other solvent molecules.
do i=1,nPol
  Indp1 = i+iCnum*nPol
  Indco1 = i+iCNum*nCent
  do j=iCnum+2,nPart
    IndP1 = IndP1+nPol
    IndCo1 = IndCo1+nCent
    IndTr = nTri_Elem(j-(iCnum+2))
    do k=1,nPol
      IndP2 = k+(iCnum-1)*nPol
      IndCo2 = k+(iCNum-1)*nCent
      do l=iCnum+1,j-1
        IndTr1 = Indtr1+1
        Indp2 = Indp2+nPol
        Indco2 = Indco2+nCent
        IndTri = IndTr+l-iCnum
        xx(Indp1,Indp2) = (Cordst(1,indco1)-Cordst(1,indco2))*Dist(k,i,IndTri)
        yy(Indp1,Indp2) = (Cordst(2,indco1)-Cordst(2,indco2))*Dist(k,i,IndTri)
        zz(Indp1,Indp2) = (Cordst(3,indco1)-Cordst(3,indco2))*Dist(k,i,IndTri)
        rr3(IndP1,IndP2) = Dist(k,i,IndTri)**3
        ! Why should xx(indp2,indp1)=xx(indp1,indp2), you wonder, should
        ! they not be of different sign? The answer is, it does not
        ! matter. Recall that the formula for the field from an ideal
        ! dipole changes sign twice, thus no time in effect, when
        ! the sign of the r-vector is changed.
        xx(Indp2,Indp1) = xx(Indp1,Indp2)
        yy(Indp2,Indp1) = yy(Indp1,Indp2)
        zz(Indp2,Indp1) = zz(Indp1,Indp2)
        rr3(Indp2,Indp1) = rr3(Indp1,Indp2)
      end do
    end do
  end do
end do
Gri(:,:) = Zero
do i=1,nPol
  k = i+(iCnum-1)*nCent
  do i1=iCnum+1,nPart
    k = k+nCent
    imd = (i1-1)*nPol+i
    do j=1,nPol
      l = j+(iCnum-1)*nCent
      jnd = ((i1-(iCnum+1))*nCent+i-1)*ncparm+j-nCent
      do j1=iCnum+1,nPart
        l = l+nCent
        ild = (j1-1)*nPol+j
        jnd = jnd+nCent
        Gri(imd,ild) = DistIm(j,j1-iCnum,i,i1-iCnum)
        xxi(imd,ild) = (Cordim(1,k)-Cordst(1,l))*DistIm(j,j1-iCnum,i,i1-iCnum)
        yyi(imd,ild) = (Cordim(2,k)-Cordst(2,l))*DistIm(j,j1-iCnum,i,i1-iCnum)
        zzi(imd,ild) = (Cordim(3,k)-Cordst(3,l))*DistIm(j,j1-iCnum,i,i1-iCnum)
      end do
    end do
  end do
end do

return

end subroutine PolPrep
