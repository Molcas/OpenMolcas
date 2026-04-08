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
! Copyright (C) 1990,2005, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine PLF_RI_2(AOint,ijkl,jCmp,lCmp,iAO,iAOst,jBas,lBas,kOp,TInt,nTInt,iSO2Ind,iOffA,nSOs)
!***********************************************************************
!                                                                      *
!  object: to sift and index the petite list format integrals.         *
!                                                                      *
!          the indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          May '90                                                     *
!          Modified to 2-center RI June '05                            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Constants, only: One
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ijkl, jCmp, lCmp, iAO(4), iAOst(4), jBas, lBas, kOp(4), nTInt, nSOs, iSO2Ind(nSOs), iOffA(4)
real(kind=wp), intent(in) :: AOint(ijkl,jCmp,lCmp)
real(kind=wp), intent(_OUT_) :: TInt(nTInt)
integer(kind=iwp) :: i2, i4, iAOj, iAOl, iAOstj, iAOstl, ij, iOff, iOffA_, iSO, jSO, jSOj, lSO, lSOl, mm_, mx, nijkl, nn
#ifdef _DEBUGPRINT_
real(kind=wp) :: r1, r2
real(kind=wp), external :: ddot_
#endif

! Note on the ordering of the basis functions (valence or auxiliary).
!
! 1) The basis functions are first ordered according to the order in which the basis sets appear in the input.
!    Hence, for example, for a heterodiatomic molecule the order of the basis sets for the centers A and B occurs
!    according to which center is first in the input.
! 2) For each basis set the basis functions are ordered in shells with fixed total angular momentum.
!    These shells are given running shell indices, which are used to identify a particular angular shell, iS.
! 3) For a given total angular momentum, iAng, the basis sets are ordered with respect to the elements of the
!    shell. For example, for l=2, we have the basis functions ordered as m=2,1,0,-1,-2.
! 4) For each component the basis functions finally runs over all contracted basis functions of this particular shell.

! For the 2-center integrals shell 1 and 3 are dummy shells.

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
r1 = DDot_(ijkl*jCmp*lCmp,AOInt,1,[One],0)
r2 = DDot_(ijkl*jCmp*lCmp,AOInt,1,AOInt,1)
write(u6,*) ' Sum=',r1
write(u6,*) ' Dot=',r2
call RecPrt(' In Plf_RI_2: AOInt',' ',AOInt,ijkl,jCmp*lCmp)
#endif

! If the code cannot handle all the contracted functions at the same time the loop over them is partitioned.
! iAOst(:) list the starting index -- if no partitioning the value is 0.
iAOstj = iAOst(2)
iAOstl = iAOst(4)

! This is an index to be used in iAOtSO(:,:). This table lists the canonical index of the first basis function with a
! particular angular momentum element. As the basis is ordered the next contracted basis function with the same m
! value simply increment this value by one.
!
! NOTE: the auxiliary functions have canonical indices starting after the valence functions!

iAOj = iAO(2)
iAOl = iAO(4)

iOff = nBas(0) ! Number of valence functions.

iOffA_ = iOffA(1)  ! Offset to where the subsection of the A matrix for the jS shell starts
mm_ = iOffA(4)     ! mm_ (iOffA(4)): is the number of basis functions up to and including shell jS
nn = mm_-iOffA(2)  ! nn:  the number of basis functions before shell jS
                   ! iOffA(2) is the number of basis functions of shell jS
mx = nTri_Elem(nn) ! mx: is the number of elements in A before the block belonging to shell jS

#ifdef _DEBUGPRINT_
write(u6,*) 'nn,mx,mm_=',nn,mx,mm_
write(u6,*) 'mx-iOffA_=',mx-iOffA_
write(u6,*) 'iOff=',iOff
write(u6,*) 'lBas,jBas=',lBas,jBas
write(u6,*) 'lCmp,jCmp=',lCmp,jCmp
#endif

do i2=1,jCmp
  jSO = iAOtSO(iAOj+i2,kOp(2))+iAOstj-iOff ! starting canonical index for shell jS
  do i4=1,lCmp
    lSO = iAOtSO(iAOl+i4,kOp(4))+iAOstl-iOff ! starting canonical index for shell lS

    nijkl = 0
    do lSOl=lSO,lSO+lBas-1

      do jSOj=jSO,jSO+jBas-1

        iSO = iSO2Ind(jSOj)+nn
        ij = iTri(iSO,lSOl)-mx+iOffA_

        nijkl = nijkl+1
        TInt(ij) = AOint(nijkl,i2,i4)

      end do
    end do

  end do
end do

end subroutine PLF_RI_2
