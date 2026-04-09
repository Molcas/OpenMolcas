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
subroutine PGet4(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,PSOPam,n1,n2,n3,n4,iPam,MapPam,mDim,Cred,nCred,Scr1,nScr1,Scr2, &
                 nScr2,PMax)
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
!             Modified from PGet2, October '92.                        *
!***********************************************************************

use SOAO_Info, only: iAOtSO, iOffSO
use pso_stuff, only: lSA
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), ijkl, nPSO, n1, n2, n3, n4, mDim, nCred, &
                                 nScr1, nScr2
real(kind=wp), intent(out) :: PSO(ijkl,nPSO), PSOPam(n1,n2,n3,n4), iPam(n1+n2+n3+n4), MapPam(4,mDim), Cred(nCred), Scr1(nScr1,2), &
                              Scr2(nScr2), PMax
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, iiBas(4), in1, in2, iS, iSO, iSOi, iSym(0:7), j, j1, j12, j123, j2, j3, j4, jAOj, jPam, &
                     jS, jSO, jSOj, jSym(0:7), k1, k2, k3, k4, kAOk, kS, kSO, kSOk, kSym(0:7), lAOl, lS, lSO, lSOl, lSym(0:7), &
                     MemSO2, nijkl, niSym, njSym, nkSym, nlSym, nPam(4,0:7), nPSOPam

! Prepare some data for Pam

iiBas(1) = iBas
iiBas(2) = jBas
iiBas(3) = kBas
iiBas(4) = lBas
nPSOPam = n1*n2*n3*n4

! Set up table with SO indices in iPam and a table
! with number of basis functions in each irrep in nPam.
! Observe that the SO index is only within a given irrep.

nPam(:,:) = 0
in1 = 0
do jPam=1,4
  in2 = 0
  do j=0,nIrrep-1
    do i1=1,iCmp(jPam)
      if (iAOtSO(iAO(jPam)+i1,j) > 0) then
        iSO = iAOtSO(iAO(jPam)+i1,j)+iAOst(jPam)
        nPam(jPam,j) = nPam(jPam,j)+iiBas(jPam)
        do iAOi=0,iiBas(jPam)-1
          iSOi = iSO+iAOi
          in2 = in2+1
          iPam(in1+in2) = real(iSOi,kind=wp)
          MapPam(jPam,iSOi+iOffSO(j)) = real(in2,kind=wp)
        end do
      end if
    end do
  end do
  in1 = in1+in2
end do

! Get the scrambled 2nd order density matrix

if (LSA) then
  call PTrans_sa(nPam,iPam,n1+n2+n3+n4,PSOPam,nPSOPam,Cred,nCred/2,Scr1(:,1),nScr1,Scr2,nScr2,Scr1(:,2),nScr1)
else
  call PTrans(nPam,iPam,n1+n2+n3+n4,PSOPam,nPSOPam,Cred,nCred,Scr1,nScr1,Scr2,nScr2)
end if

! Quadruple loop over elements of the basis functions angular description.
! Observe that we will walk through the memory in AOInt in a sequential way.

PMax = Zero
MemSO2 = 0
do i1=1,iCmp(1)
  niSym = 0
  do j=0,nIrrep-1
    if (iAOtSO(iAO(1)+i1,j) > 0) then
      iSym(niSym) = j
      niSym = niSym+1
    end if
  end do
  do i2=1,iCmp(2)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do
    do i3=1,iCmp(3)
      nkSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(3)+i3,j) > 0) then
          kSym(nkSym) = j
          nkSym = nkSym+1
        end if
      end do
      do i4=1,iCmp(4)
        nlSym = 0
        do j=0,nIrrep-1
          if (iAOtSO(iAO(4)+i4,j) > 0) then
            lSym(nlSym) = j
            nlSym = nlSym+1
          end if
        end do

        ! Loop over irreps which are spanned by the basis function.
        ! Again, the loop structure is restricted to ensure unique
        ! integrals.

        do is=0,niSym-1
          j1 = iSym(is)

          do js=0,njSym-1
            j2 = jSym(js)
            j12 = ieor(j1,j2)

            do ks=0,nkSym-1
              j3 = kSym(ks)
              j123 = ieor(j12,j3)
              do ls=0,nlSym-1
                j4 = lSym(ls)
                if (j123 /= j4) cycle

                MemSO2 = MemSO2+1

                ! Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)

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
                        nijkl = nijkl+1
                        k1 = int(MapPam(1,iSOi))

                        ! Pick up the contribution.

                        PMax = max(PMax,abs(PSOPam(k1,k2,k3,k4)))
                        PSO(nijkl,MemSO2) = PSOPam(k1,k2,k3,k4)

                      end do
                    end do
                  end do
                end do

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do
if (nPSO /= MemSO2) then
  call WarningMessage(2,'PGet4: nPSO /= MemSO2')
  call Abend()
end if

return

end subroutine PGet4
