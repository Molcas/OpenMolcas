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
subroutine PGet2(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,nijkl,PSO,nPSO,DSO,DSSO,nDSO,ExFac,CoulFac,PMax)
!***********************************************************************
!                                                                      *
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density matrix.                 *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!***********************************************************************

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, Quart
use Definitions, only: wp, iwp, u6
#ifdef _DEBUGPRINT_
use pso_stuff, only: D0, iD0Lbl
#endif

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), nijkl, nPSO, nDSO
real(kind=wp), intent(out) :: PSO(nijkl,nPSO), PMax
real(kind=wp), intent(in) :: DSO(nDSO), DSSO(nDSO), ExFac, CoulFac
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, IndIJ, IndIK, IndIL, IndJK, IndJL, IndKL, iS, iSO, iSOi, iSym(0:7), j, j1, j12, j123, &
                     j2, j3, j4, jAOj, jS, jSO, jSOj, jSym(0:7), kAOk, kS, kSO, kSOk, kSym(0:7), lAOl, lOper, lS, lSO, lSOl, &
                     lSym(0:7), MemSO2, mijkl, niSym, njSym, nkSym, nlSym
real(kind=wp) :: t14, Temp
integer(kind=iwp), external :: iPntSO
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iComp
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iComp = 1
call PrMtrx(' In PGet2:DSO ',[iD0Lbl],iComp,1,D0)
#endif
lOper = 1
t14 = Quart*ExFac
PMax = Zero

! Quadruple loop over elements of the basis functions angular
! description.
! Observe that we will walk through the memory in AOInt in a
! sequential way.

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
            j12 = ieor(j1,j2)

            do ks=0,nkSym-1
              j3 = kSym(ks)
              j123 = ieor(j12,j3)
              do ls=0,nlSym-1
                j4 = lSym(ls)
                if (j123 /= j4) cycle

                MemSO2 = MemSO2+1

                ! Unfold the way the eight indices have been reordered.
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

                        ! Contribution D(ij)*D(kl) to P(ijkl)
                        if (j1 == j2) then
                          ! j3 == j4 also
                          Indij = iPntSO(j1,j2,lOper,nbas)+iTri(iSOi,jSOj)
                          Indkl = iPntSO(j3,j4,lOper,nbas)+iTri(kSOk,lSOl)
                          temp = DSO(Indij)*DSO(Indkl)*CoulFac
                        else
                          temp = Zero
                        end if

                        ! Contribution -1/4*D(ik)*D(jl) to P(ijkl)
                        if (j1 == j3) then
                          ! j2 == j4 also
                          Indik = iPntSO(j1,j3,lOper,nbas)+iTri(iSOi,kSOk)
                          Indjl = iPntSO(j2,j4,lOper,nbas)+iTri(jSOj,lSOl)
                          temp = temp-t14*(DSO(Indik)*DSO(Indjl)+DSSO(Indik)*DSSO(Indjl))
                        end if

                        ! Contribution -1/4*D(il)*D(jk) to P(ijkl)
                        if (j1 == j4) then
                          ! j2 == j3 also
                          Indil = iPntSO(j1,j4,lOper,nbas)+iTri(iSOi,lSOl)
                          Indjk = iPntSO(j2,j3,lOper,nbas)+iTri(jSOj,kSOk)
                          temp = temp-t14*(DSO(Indil)*DSO(Indjk)+DSSO(Indil)*DSSO(Indjk))
                        end if

                        PMax = max(PMax,abs(Temp))
                        PSO(mijkl,MemSO2) = temp

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
  call WarningMessage(2,' PGet2: nPSO /= MemSO2')
  write(u6,*) nPSO,MemSO2
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet2:PSO ',' ',PSO,nijkl,nPSO)
#endif

return

end subroutine PGet2
