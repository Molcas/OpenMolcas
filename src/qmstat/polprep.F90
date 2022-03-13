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

subroutine PolPrep(iDist,iDistIm,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,iCNum,nSize)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iDist, iDistIm, iCNum, nSize
real(kind=wp) :: xx(nSize,nSize), yy(nSize,nSize), zz(nSize,nSize), rr3(nSize,nSize), xxi(nSize,nSize), yyi(nSize,nSize), &
                 zzi(nSize,nSize), Gri(nSize,nSize)
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1, ii, ild, imd, Indco1, IndCo2, Indp1, IndP2, IndTr, IndTr1, IndTri, j, j1, jj, jnd, k, l, ncParm

!----------------------------------------------------------------------*
! Simply compute some vectors etc. for the ensuing polarization        *
! calculation.                                                         *
!----------------------------------------------------------------------*
ncParm = ncent*npart-ncent*icNum
do i=nPol*iCNum+1,nPol*nPart !Loop over solvent polarization sites.
  do j=nPol*iCNum+1,nPol*nPart
    rr3(i,j) = Zero
  end do
end do
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
    IndTr = ((j-(iCnum+2))*(j-(iCNum+1)))/2*nCent**2+(i-1)*nCent
    do k=1,nPol
      IndP2 = k+(iCnum-1)*nPol
      IndCo2 = k+(iCNum-1)*nCent
      do l=iCnum+1,j-1
        IndTr1 = Indtr1+1
        Indp2 = Indp2+nPol
        Indco2 = Indco2+nCent
        IndTri = IndTr+(l-(iCnum+1))*nCent**2+k
        xx(Indp1,Indp2) = (Cordst(indco1,1)-Cordst(indco2,1))*Work(iDist+indtri-1)
        yy(Indp1,Indp2) = (Cordst(indco1,2)-Cordst(indco2,2))*Work(iDist+indtri-1)
        zz(Indp1,Indp2) = (Cordst(indco1,3)-Cordst(indco2,3))*Work(iDist+indtri-1)
        rr3(IndP1,IndP2) = Work(iDist+Indtri-1)**3
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
do ii=1,nSize
  do jj=1,nSize
    Gri(ii,jj) = Zero
  end do
end do
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
        Gri(imd,ild) = Work(iDistIm+jnd-1)
        xxi(imd,ild) = (Cordim(k,1)-Cordst(l,1))*Work(iDistIm+jnd-1)
        yyi(imd,ild) = (Cordim(k,2)-Cordst(l,2))*Work(iDistIm+jnd-1)
        zzi(imd,ild) = (Cordim(k,3)-Cordst(l,3))*Work(iDistIm+jnd-1)
      end do
    end do
  end do
end do

return

end subroutine PolPrep
