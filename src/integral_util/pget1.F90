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
subroutine PGet1(PAO,ijkl,nPAO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,DSO,DSSO,nDSO,ExFac,CoulFac,PMax)
!***********************************************************************
!                                                                      *
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density.                        *
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
use Constants, only: Zero, Quart
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use pso_stuff, only: D0, iD0Lbl
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4), nDSO
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PMax
real(kind=wp), intent(in) :: DSO(nDSO), DSSO(nDSO), ExFac, CoulFac
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, IndIJ, IndIK, IndIL, IndJK, IndJL, IndKL, iPAO, iSO, iSOi, jAOj, jSO, jSOj, kAOk, kSO, &
                     kSOk, lAOl, lSO, lSOl, nijkl
real(kind=wp) :: t14, Temp
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, iComp
real(kind=wp), external :: DDot_
#endif

#ifdef _DEBUGPRINT_
iComp = 1
call PrMtrx('DSO     ',[iD0Lbl],iComp,1,D0)
write(u6,*) ' nBases..=',iBas,jBas,kBas,lBas
#endif

! Quadruple loop over elements of the basis functions angular
! description.
! Observe that we will walk through the memory in PAO in a
! sequential way.

PMax = Zero
iPAO = 0
t14 = Quart*ExFac
do i1=1,iCmp(1)
  do i2=1,iCmp(2)
    do i3=1,iCmp(3)
      do i4=1,iCmp(4)

        ! Unfold the way the eight indices have been reordered.
        iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
        jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
        kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0
        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj
              do iAOi=0,iBas-1
                iSOi = iSO+iAOi
                nijkl = nijkl+1

                ! D(ij)*D(kl)

                Indij = iTri(iSOi,jSOj)
                Indkl = iTri(kSOk,lSOl)
                temp = DSO(Indij)*DSO(Indkl)*CoulFac

                ! -0.25*D(ik)*D(jl)

                Indik = iTri(iSOi,kSOk)
                Indjl = iTri(jSOj,lSOl)
                temp = temp-t14*(DSO(Indik)*DSO(Indjl)+DSSO(Indik)*DSSO(Indjl))

                ! -0.25*D(il)*D(jk)

                Indil = iTri(iSOi,lSOl)
                Indjk = iTri(jSOj,kSOk)
                temp = temp-t14*(DSO(Indil)*DSO(Indjk)+DSSO(Indil)*DSSO(Indjk))

                PMax = max(PMax,abs(Temp))
                PAO(nijkl,iPAO) = temp

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do
if (iPAO /= nPAO) then
  call WarningMessage(2,' Error in PGet1!')
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet1:PAO ',' ',PAO,ijkl,nPAO)
do i=1,ijkl
  write(u6,*) DDot_(nPAO,PAO(i,1),ijkl,PAO(i,1),ijkl)
end do
#endif

return

end subroutine PGet1
