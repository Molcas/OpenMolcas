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
! Copyright (C) 1992,2000, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PGet1_Aces(PAO,ijkl,nPAO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,DSO,DSO_Var,DSSO,DSSO_Var,nDSO,Gmma,nGamma,iSO2cI,nSOs, &
                      iSO2Sh,PMax)
!***********************************************************************
!                                                                      *
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density.                        *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          DSO: HF 1st order density                                   *
!          DSO_Var: 1st order density of correlated wf.                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!     Modified to Aces 2 by RL, July 2000, Gainesville, FL, USA        *
!***********************************************************************

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use pso_stuff, only: Gamma_MRCISD
use Constants, only: Zero, One, Quart, Four
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use pso_stuff, only: iD0Lbl, DVar, D0
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4), nDSO, nGamma, nSOs, &
                                 iSO2cI(2,nSOs), iSO2Sh(nSOs)
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PMax
real(kind=wp), intent(in) :: DSO(nDSO), DSO_Var(nDSO), DSSO(nDSO), DSSO_Var(nDSO), Gmma(nGamma)
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, Index_A, Index_AB, Index_ABCD, Index_B, Index_C, Index_CD, Index_D, Indij, Indik, &
                     Indil, Indjk, Indjl, Indkl, iPAO, iShell_A, iShell_AB, iShell_B, iShell_C, iShell_CD, iShell_D, iSO, iSOi, &
                     jAOj, jSO, jSOj, kAOk, kSO, kSOk, lAOl, lSO, lSOl, nDim_A, nDim_AB, nDim_B, nDim_C, nDim_CD, nDim_D, nijkl
real(kind=wp) :: t14, Temp
real(kind=wp), parameter :: exfac = One
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, iComp
real(kind=wp), external :: DDOt_
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iComp = 1
call PrMtrx('DSO     ',[iD0Lbl],iComp,1,D0)
call PrMtrx('DSO_Var ',[iD0Lbl],iComp,1,DVar)
write(u6,*) ' nBases..=',iBas,jBas,kBas,lBas
write(u6,*) 'iSO2Sh=',iSO2Sh
write(u6,*) 'iSO2cI(1)',(iSO2cI(1,i),i=1,nSOs)
write(u6,*) 'iSO2cI(2)',(iSO2cI(2,i),i=1,nSOs)
call RecPrt('PGet1: Gmma',' ',Gmma,1,nGamma)
#endif

! Quadruple loop over elements of the basis functions angular
! description.
! Observe that we will walk through the memory in PAO in a
! sequential way.

PMax = Zero
iPAO = 0
t14 = Quart*exfac
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
          iShell_D = iSO2Sh(lSOl)
          Index_D = iSO2cI(1,lSOl)
          nDim_D = iSO2cI(2,lSOl)
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk
            iShell_C = iSO2Sh(kSOk)
            Index_C = iSO2cI(1,kSOk)
            nDim_C = iSO2cI(2,kSOk)
            nDim_CD = nDim_C*nDim_D
            iShell_CD = iTri(iShell_C,iShell_D)
            if (iShell_C > iShell_D) then
              Index_CD = (Index_D-1)*nDim_C+Index_C
            else if (iShell_C == iShell_D) then
              Index_CD = iTri(Index_C,Index_D)
            else
              Index_CD = (Index_C-1)*nDim_D+Index_D
            end if
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj
              iShell_B = iSO2Sh(jSOj)
              Index_B = iSO2cI(1,jSOj)
              nDim_B = iSO2cI(2,jSOj)
              do iAOi=0,iBas-1
                iSOi = iSO+iAOi
                iShell_A = iSO2Sh(iSOi)
                Index_A = iSO2cI(1,iSOi)
                nDim_A = iSO2cI(2,iSOi)
                nDim_AB = nDim_A*nDim_B
                iShell_AB = iTri(iShell_A,iShell_B)
                if (iShell_A > iShell_B) then
                  Index_AB = (Index_B-1)*nDim_A+Index_A
                else if (iShell_A == iShell_B) then
                  Index_AB = iTri(Index_A,Index_B)
                else
                  Index_AB = (Index_A-1)*nDim_B+Index_B
                end if
                if (iShell_AB > iShell_CD) then
                  Index_ABCD = (Index_CD-1)*nDim_AB+Index_AB
                else if (iShell_AB == iShell_CD) then
                  Index_ABCD = iTri(Index_AB,Index_CD)
                else
                  Index_ABCD = (Index_AB-1)*nDim_CD+Index_CD
                end if
                nijkl = nijkl+1

                !*********** columbus interface ****************************************
                !do not reconstruct the two-particle density from the one-particle
                !density or partial two-particle densities but simply read them from
                !file
                if (gamma_mrcisd) then
                  temp = Gmma(Index_ABCD)
                else

                  ! D(ij)*D(kl)

                  Indij = iTri(iSOi,jSOj)
                  Indkl = iTri(kSOk,lSOl)
                  temp = DSO(Indij)*DSO(Indkl)+(DSO_Var(Indij)-DSO(Indij))*DSO(Indkl)+DSO(Indij)*(DSO_Var(Indkl)-DSO(Indkl))

                  ! -0.25*D(ik)*D(jl)

                  Indik = iTri(iSOi,kSOk)
                  Indjl = iTri(jSOj,lSOl)
                  temp = temp-t14*(DSO(Indik)*DSO(Indjl)+(DSO_Var(Indik)-DSO(Indik))*DSO(Indjl)+ &
                                   DSO(Indik)*(DSO_Var(Indjl)-DSO(Indjl))+DSSO(Indik)*DSSO(Indjl)+ &
                                   (DSSO_Var(Indik)-DSSO(Indik))*DSSO(Indjl)+DSSO(Indik)*(DSSO_Var(Indjl)-DSSO(Indjl)))

                  ! -0.25*D(il)*D(jk)

                  Indil = iTri(iSOi,lSOl)
                  Indjk = iTri(jSOj,kSOk)
                  temp = temp-t14*(DSO(Indil)*DSO(Indjk)+(DSO_Var(Indil)-DSO(Indil))*DSO(Indjk)+ &
                                   DSO(Indil)*(DSO_Var(Indjk)-DSO(Indjk))+DSSO(Indil)*DSSO(Indjk)+ &
                                   (DSSO_Var(Indil)-DSSO(Indil))*DSSO(Indjk)+DSSO(Indil)*(DSSO_Var(Indjk)-DSSO(Indjk)))

                  temp = temp+Four*Gmma(Index_ABCD)
                end if

                PMax = max(PMax,abs(temp))
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
  call WarningMessage(2,' Error in PGet1_Aces!')
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet1:PAO ',' ',PAO,ijkl,nPAO)
do i=1,ijkl
  write(u6,*) DDot_(nPAO,PAO(i,1),ijkl,PAO(i,1),ijkl)
end do
#endif

return

end subroutine PGet1_Aces
