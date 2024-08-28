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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Trnsps_Seward(ijCmp,iCmp,jCmp,iAng,jAng,iShll,jShll,kOp,ijkl,ij,AOInt,Scrtch)
!***********************************************************************
!                                                                      *
!  Object: to transpose the integrals in order to resolve the          *
!          redundancy (faA,fbB)=(fcC,fdD). In this case both sides will*
!          have the same DCR, i.e. (R)=(S). In this case we will only  *
!          need the unique combinations. For the off diagonal comb-    *
!          inations (R=/=S) we will pick up two terms. The terms are   *
!          (faA,fbR(B)|faT(A),fbTS(B)) and                             *
!          (faA,fbS(B)|faT(A),fbTR(B)). Since T and T-1 are the same   *
!          in D2h it is simple to see that after applying T on the     *
!          second integral we will end up with the first one.          *
!                                                                      *
!          However, since we compute the integrals in batches there    *
!          will not be a simple one to one correspondes between the    *
!          integrals in batch one and two. But after transposing the   *
!          pair arguments we will achive that one to one correspondens.*
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************

use Basis_Info, only: Shells
use Real_Spherical, only: iSphCr
use Symmetry_Info, only: iChBas
use Constants, only: One
use Definitions, only: wp

implicit none
integer ijCmp, iCmp, jCmp, iAng, jAng, iShll, jShll, kOp, ijkl, ij
real*8 AOInt(ijkl,ijCmp,ijCmp), Scrtch(ijkl,ijCmp,ijCmp)
integer, external :: iPrmt
integer ixyz, iOff
integer ii, jj, i1, i2, i3, i4, iChBs, jChBs, kChBs, lChBs, ij2, i12, i34, ij1
real*8 pa1T, pb1T, pa2T, pb2T, Factor
! Statement Function
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

!call RecPrt(' In Trnsps: AOInt ',' ',AOInt,ijkl,ijCmp*ijCmp)

! Change phase factor. This is only necessary if T=/=E.

if ((kOp == 0) .or. (ijCmp == 0)) Go To 14
ii = iOff(iAng)
jj = iOff(jAng)
do i1=1,iCmp
  iChBs = iChBas(ii+i1)
  if (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pa1T = real(iPrmt(kOp,iChBs),kind=wp)
  do i2=1,jCmp
    jChBs = iChBas(jj+i2)
    if (Shells(jShll)%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pb1T = real(iPrmt(kOp,jChBs),kind=wp)
    ij1 = iCmp*(i2-1)+i1

    do i3=1,iCmp
      kChBs = iChBas(ii+i3)
      if (Shells(iShll)%Transf) kChBs = iChBas(iSphCr(ii+i3))
      pa2T = real(iPrmt(kOp,kChBs),kind=wp)
      do i4=1,jCmp
        lChBs = iChBas(jj+i4)
        if (Shells(jShll)%Transf) lChBs = iChBas(iSphCr(jj+i4))
        pb2T = real(iPrmt(kOp,lChBs),kind=wp)
        ij2 = iCmp*(i4-1)+i3
        Factor = pa1T*pb1T*pa2T*pb2T
        if (Factor /= One) call DScal_(ijkl,Factor,AOInt(1,ij1,ij2),1)
      end do
    end do
  end do
end do
14 continue

! Transpose ijkl,abcd to klij,cdab

if ((ijCmp == 1) .or. (ij == 1)) then
  call DGeTMI(AOInt,ijCmp*ij,ijCmp*ij)
else
  do i12=1,ijCmp
    do i34=1,ijCmp
      call DGeTMO(AOInt(1,i12,i34),ij,ij,ij,Scrtch(1,i34,i12),ij)

    end do
  end do
  call dcopy_(ijkl*ijCmp*ijCmp,Scrtch,1,AOInt,1)
end if

!call RecPrt(' Exit Trnsps: AOInt ',' ',AOInt,ijkl,ijCmp*ijCmp)

end subroutine Trnsps_Seward
