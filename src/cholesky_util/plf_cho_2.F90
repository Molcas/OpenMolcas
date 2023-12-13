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

subroutine PLF_Cho_2(TInt,lInt,AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
!***********************************************************************
!                                                                      *
!  object: to sift and index the petite list format integrals.         *
!                                                                      *
!          the indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          May '90                                                     *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use Cholesky, only: iShlSO, iShP2Q, iShP2RS, iSOShl, LuPri, nBstSh, nnBstR, ShA, ShAB, ShB, ShC, ShCD, ShD
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp), intent(in) :: lInt, ijkl, iCmp, jCmp, kCmp, lCmp, iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4)
real(kind=wp), intent(inout) :: TInt(lInt)
real(kind=wp), intent(in) :: AOint(ijkl,iCmp,jCmp,kCmp,lCmp)
integer(kind=iwp) :: A, AB, ABCD, B, C, CD, CDAB, D, i1, i2, i3, i4, IAB, iAOi, iAOj, iAOk, iAOl, iAOsti, iAOstj, iAOstk, iAOstl, &
                     ICD, irout, ISHLAB, ISHLCD, ISHLI, ISHLJ, ISHLK, ISHLL, iSO, iSOi, iSOs(4), jprint, jSO, jSOj, kSO, kSOk, &
                     lSO, lSOl, nijkl, NTELM, NUMA, NUMB, NUMC, NUMD
real(kind=wp) :: r1, r2
real(kind=wp), external :: ddot_

irout = 109
jprint = nprint(irout)
if (jPrint >= 49) then
  r1 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
  r2 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
  write(u6,*) ' Sum=',r1
  write(u6,*) ' Dot=',r2
end if
if (jPrint >= 99) call RecPrt(' In Plf_Cho_2: AOInt',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)

NUMC = NBSTSH(SHC)
NUMD = NBSTSH(SHD)
NUMA = NBSTSH(SHA)
NUMB = NBSTSH(SHB)

ISHLCD = SHCD
ISHLAB = SHAB

! to avoid stupid compiler warnings:

C = 0
D = 0
A = 0
B = 0

NTELM = 0

! Allocate space to store integrals to gether with their
! Symmetry batch and sequence number.
! To avoid conflicts in using memory this is done in the
! subroutine PSOAO

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

iAOsti = iAOst(1)
iAOstj = iAOst(2)
iAOstk = iAOst(3)
iAOstl = iAOst(4)
iAOi = iAO(1)
iAOj = iAO(2)
iAOk = iAO(3)
iAOl = iAO(4)

do i1=1,iCmp
  iSOs(1) = iAOtSO(iAOi+i1,kOp(1))+iAOsti
  do i2=1,jCmp
    iSOs(2) = iAOtSO(iAOj+i2,kOp(2))+iAOstj
    do i3=1,kCmp
      iSOs(3) = iAOtSO(iAOk+i3,kOp(3))+iAOstk
      do i4=1,lCmp
        iSOs(4) = iAOtSO(iAOl+i4,kOp(4))+iAOstl

        iSO = iSOs(1)
        jSO = iSOs(2)
        kSO = iSOs(3)
        lSO = iSOs(4)

        nijkl = 0
        do lSOl=lSO,lSO+lBas-1
          do kSOk=kSO,kSO+kBas-1
            do jSOj=jSO,jSO+jBas-1
              do iSOi=iSO,iSO+iBas-1

                nijkl = nijkl+1
                NTELM = NTELM+1

                ISHLI = ISOSHL(ISOI)
                ISHLJ = ISOSHL(JSOJ)
                ISHLK = ISOSHL(KSOK)
                ISHLL = ISOSHL(LSOL)

                if ((ISHLI == SHC) .and. (ISHLJ == SHD) .and. (ISHLK == SHA) .and. (ISHLL == SHB)) then
                  C = ISHLSO(ISOI)
                  D = ISHLSO(JSOJ)
                  A = ISHLSO(KSOK)
                  B = ISHLSO(LSOL)
                else if ((ISHLJ == SHC) .and. (ISHLI == SHD) .and. (ISHLK == SHA) .and. (ISHLL == SHB)) then
                  C = ISHLSO(JSOJ)
                  D = ISHLSO(ISOI)
                  A = ISHLSO(KSOK)
                  B = ISHLSO(LSOL)
                else if ((ISHLI == SHC) .and. (ISHLJ == SHD) .and. (ISHLL == SHA) .and. (ISHLK == SHB)) then
                  C = ISHLSO(ISOI)
                  D = ISHLSO(JSOJ)
                  A = ISHLSO(LSOL)
                  B = ISHLSO(KSOK)
                else if ((ISHLJ == SHC) .and. (ISHLI == SHD) .and. (ISHLL == SHA) .and. (ISHLK == SHB)) then
                  C = ISHLSO(JSOJ)
                  D = ISHLSO(ISOI)
                  A = ISHLSO(LSOL)
                  B = ISHLSO(KSOK)
                else if ((ISHLK == SHC) .and. (ISHLL == SHD) .and. (ISHLI == SHA) .and. (ISHLJ == SHB)) then
                  C = ISHLSO(KSOK)
                  D = ISHLSO(LSOL)
                  A = ISHLSO(ISOI)
                  B = ISHLSO(JSOJ)
                else if ((ISHLL == SHC) .and. (ISHLK == SHD) .and. (ISHLI == SHA) .and. (ISHLJ == SHB)) then
                  C = ISHLSO(LSOL)
                  D = ISHLSO(KSOK)
                  A = ISHLSO(ISOI)
                  B = ISHLSO(JSOJ)
                else if ((ISHLK == SHC) .and. (ISHLL == SHD) .and. (ISHLJ == SHA) .and. (ISHLI == SHB)) then
                  C = ISHLSO(KSOK)
                  D = ISHLSO(LSOL)
                  A = ISHLSO(JSOJ)
                  B = ISHLSO(ISOI)
                else if ((ISHLL == SHC) .and. (ISHLK == SHD) .and. (ISHLJ == SHA) .and. (ISHLI == SHB)) then
                  C = ISHLSO(LSOL)
                  D = ISHLSO(KSOK)
                  A = ISHLSO(JSOJ)
                  B = ISHLSO(ISOI)
                else
                  write(LUPRI,*) 'Shell quadruple requested: ',SHC,SHD,SHA,SHB
                  write(LUPRI,*) 'Shell quadruple of element ',NTELM,':',ISHLI,ISHLJ,ISHLK,ISHLL
                  call CHO_QUIT('Logical error in PLF_Cho_2',103)
                end if

                if (SHA == SHB) then
                  AB = ITRI(A,B)
                else
                  AB = NUMA*(B-1)+A
                end if
                if (SHC == SHD) then
                  CD = ITRI(C,D)
                else
                  CD = NUMC*(D-1)+C
                end if

                ICD = ISHP2RS(1,CD)
                IAB = ISHP2Q(1,AB)
                if ((ICD > 0) .and. (IAB > 0)) then
                  CDAB = nnBstR(1,2)*(IAB-1)+ICD
                  TINT(CDAB) = AOint(nijkl,i1,i2,i3,i4)
                end if

                if (ISHLCD == ISHLAB) then
                  if ((SHC == SHD) .or. (SHC == SHA)) then
                    IAB = ISHP2RS(1,AB)
                    ICD = ISHP2Q(1,CD)
                    if ((ICD > 0) .and. (IAB > 0)) then
                      ABCD = nnBstR(1,2)*(ICD-1)+IAB
                      TINT(ABCD) = AOint(nijkl,i1,i2,i3,i4)
                    end if
                  else if (SHC == SHB) then
                    AB = NUMB*(A-1)+B
                    CD = NUMD*(C-1)+D
                    IAB = ISHP2RS(1,AB)
                    ICD = ISHP2Q(1,CD)
                    if ((ICD > 0) .and. (IAB > 0)) then
                      ABCD = nnBstR(1,2)*(ICD-1)+IAB
                      TINT(ABCD) = AOint(nijkl,i1,i2,i3,i4)
                    end if
                  end if
                end if

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine PLF_Cho_2
