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
! Copyright (C) 1992,2007, Roland Lindh                                *
!***********************************************************************

subroutine PGet2_CD3(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,nijkl,PSO,nPSO,DSO,nDSO,CoulFac,PMax,V_k,mV_k)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density matrix.                 *
!                                                                      *
!          The indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!             Modified for Cholesky 1-center gradients May 2007 by     *
!             R. Lindh                                                 *
!***********************************************************************

use Index_Functions, only: iTri
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use Symmetry_Info, only: Mul, nIrrep
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), nijkl, nPSO, nDSO, mV_k
real(kind=wp), intent(out) :: PSO(nijkl,nPSO), PMax
real(kind=wp), intent(in) :: DSO(nDSO), CoulFac, V_k(mV_k)
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, Indij, Indkl, iPntij, iPntkl, is, iSO, iSOi, iSym(0:7), j, j1, j12, j123, j2, j3, j4, &
                     jAOj, js, jSO, jSOj, jSym(0:7), kAOk, ks, kSO, kSOk, kSym(0:7), lAOl, lOper, ls, lSO, lSOl, lSym(0:7), &
                     MemSO2, mijkl, niSym, njSym, nkSym, nlSym
real(kind=wp) :: Fac, temp
integer(kind=iwp), external :: iPntSO

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iComp = 1
call PrMtrx(' In PGet2_CD3:DSO ',[iD0Lbl],iComp,1,D0)
call RecPrt('V_K',' ',V_K,1,mV_K)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!Fac = Quart
Fac = Half
lOper = 1
PMax = Zero

! Quadruple loop over elements of the basis functions angular description.
! Observe that we will walk through the memory in AOInt in a sequential way.

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

        do is=0,niSym-1
          j1 = iSym(is)

          do js=0,njSym-1
            j2 = jSym(js)
            j12 = Mul(j1+1,j2+1)-1

            do ks=0,nkSym-1
              j3 = kSym(ks)
              j123 = Mul(j12+1,j3+1)-1
              do ls=0,nlSym-1
                j4 = lSym(ls)
                if (j123 /= j4) cycle

                MemSO2 = MemSO2+1

                ! Unfold the way the eight indicies have been reordered.
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

                if ((j1 /= j2) .and. (j1 /= j3) .and. (j1 /= j4)) then
                  ! all irreps are different and the 2nd order density
                  ! matrix will be identical to zero for a SCF type wave
                  ! function.
                  PSO(:,MemSO2) = Zero
                  exit
                end if

                mijkl = 0
                do lAOl=0,lBas-1
                  lSOl = lSO+lAOl
                  do kAOk=0,kBas-1
                    kSOk = kSO+kAOk
                    do jAOj=0,jBas-1
                      jSOj = jSO+jAOj
                      do iAOi=0,iBas-1
                        iSOi = iSO+iAOi
                        mijkl = mijkl+1

                        ! Contribution V_k(ij)*D(kl) to P(ijkl)
                        if (j1 == j2) then
                          ! j3 == j4 also
                          iPntij = iPntSO(j1,j2,lOper,nbas)
                          iPntkl = iPntSO(j3,j4,lOper,nbas)
                          Indij = iPntij+iTri(iSOi,jSOj)
                          Indkl = iPntkl+iTri(kSOk,lSOl)
                          temp = V_k(Indij)*DSO(Indkl)*CoulFac
                        else
                          temp = Zero
                        end if

                        ! Contribution -1/4*D(ik)*D(jl) to P(ijkl)
                        !if (j1 == j3) then
                        !  ! j2 == j4 also
                        !end if

                        ! Contribution -1/4*D(il)*D(jk) to P(ijkl)
                        !if (j1 == j4) then
                        !  ! j2 == j3 also
                        !end if

                        PMax = max(PMax,abs(Temp))
                        PSO(mijkl,MemSO2) = Fac*temp

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
  write(u6,*) ' PGet2_CD3: nPSO /= MemSO2'
  write(u6,*) nPSO,MemSO2
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet2_CD3:PSO ',' ',PSO,nijkl,nPSO)
#endif

return

end subroutine PGet2_CD3
