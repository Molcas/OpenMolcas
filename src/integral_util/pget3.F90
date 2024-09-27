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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PGet3(PAO,ijkl,nPAO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,PAOPam,n1,n2,n3,n4,iPam,MapPam,mDim,Cred,nCred,Scr1,nScr1, &
                 Scr2,nScr2,PMax)
!***********************************************************************
!                                                                      *
!  Object: to assemble the index list of the batch of the 2nd order    *
!          density matrix.                                             *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!             Modified from PGet1, June '92                            *
!***********************************************************************

use SOAO_Info, only: iAOtSO, iOffSO
use pso_stuff, only: G_ToC, Gamma_On, lSA
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4), n1, n2, n3, n4, mDim, &
                                 nCred, nScr1, nScr2
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PAOPam(n1,n2,n3,n4), iPam(n1+n2+n3+n4), MapPam(4,mDim), Cred(nCred), Scr1(nScr1,2), &
                              Scr2(nScr2), PMax
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, iiBas(4), in1, in2, ipAO, iSO, iSOi, jAOj, jPAM, jSO, jSOj, k1, k2, k3, k4, kAOk, kSO, &
                     kSOk, lAOl, loc, lSO, lSOl, nijkl, nPam(4), nPSOPam
real(kind=wp) :: Val

#ifdef _DEBUGPRINT_
write(u6,*) ' nBases..=',iBas,jBas,kBas,lBas
#endif

! Prepare some data for Pam

iiBas(1) = iBas
iiBas(2) = jBas
iiBas(3) = kBas
iiBas(4) = lBas
nPSOPam = n1*n2*n3*n4

! Set up table with SO indices in iPam and a table
! with the number of basis functions in nPam.

nPam(:) = 0
in1 = 0
do jPam=1,4
  in2 = 0
  do i1=1,iCmp(jPam)
    iSO = iAOtSO(iAO(jPam)+i1,0)+iAOst(jPam)
    nPam(jPam) = nPam(jPam)+iiBas(jPam)
    do iAOi=0,iiBas(jPam)-1
      iSOi = iSO+iAOi
      in2 = in2+1
      iPam(in1+in2) = real(iSOi,kind=wp)
      MapPam(jPam,iSOi) = real(in2,kind=wp)
    end do
  end do
  in1 = in1+in2
end do

! Get the scrambled 2nd order density matrix

if (LSA) then
  call PTrans_sa(nPam,iPam,n1+n2+n3+n4,PAOPam,nPSOPam,Cred,nCred/2,Scr1(:,1),nScr1,Scr2,nScr2,Scr1(:,2),nScr1)
else
  call PTrans(nPam,iPam,n1+n2+n3+n4,PAOPam,nPSOPam,Cred,nCred,Scr1,nScr1,Scr2,nScr2)
end if

! Quadruple loop over elements of the basis functions angular
! description.
! Observe that we will walk through the memory in PAO in a
! sequential way.

PMax = Zero
iPAO = 0
do i1=1,iCmp(1)
  do i2=1,iCmp(2)
    do i3=1,iCmp(3)
      do i4=1,iCmp(4)

        ! Unfold the way the eight indices have been reordered.
        iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)+iOffSO(kOp(1))
        jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)+iOffSO(kOp(2))
        kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)+iOffSO(kOp(3))
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)+iOffSO(kOp(4))

        iPAO = iPAO+1
        nijkl = 0
        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          k4 = int(MapPam(4,lSOl))
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk
            k3 = int(MapPam(3,kSOk))
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj
              k2 = int(MapPam(2,jSOj))
              do iAOi=0,iBas-1
                iSOi = iSO+iAOi
                k1 = int(MapPam(1,iSOi))
                nijkl = nijkl+1

                PMax = max(PMax,abs(PAOPam(k1,k2,k3,k4)))
                PAO(nijkl,iPAO) = PAOPam(k1,k2,k3,k4)
                if (Gamma_On) then
                  Loc = k2-1+n2*(k4-1+n4*(k1-1+n1*(k3-1)))
                  Val = G_Toc(Loc+1)
                  PAO(nijkl,iPAO) = PAO(nijkl,iPAO)*1+Val
                end if

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do
if (iPAO /= nPAO) then
  call WarningMessage(2,' Error in PGet3!')
  call Abend()
end if

return

end subroutine PGet3
