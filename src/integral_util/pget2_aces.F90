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
subroutine PGet2_Aces(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,nijkl,PSO,nPSO,DSO,DSO_Var,DSSO,DSSO_Var,nDSO,Gmma,nGamma,iSO2cI,nSOs, &
                      iSO2Sh,PMax)
!***********************************************************************
!                                                                      *
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density matrix.                 *
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
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO, iOffSO
use pso_stuff, only: Gamma_MRCISD
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One, Four, Quart
use Definitions, only: wp, iwp, u6
#ifdef _DEBUGPRINT_
use pso_stuff, only: D0, DVar, iD0Lbl
#endif

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), nijkl, nPSO, nDSO, nGamma, nSOs, &
                                 iSO2cI(2,nSOs), iSO2Sh(nSOs)
real(kind=wp), intent(out) :: PSO(nijkl,nPSO), PMax
real(kind=wp), intent(in) :: DSO(nDSO), DSO_Var(nDSO), DSSO(nDSO), DSSO_Var(nDSO), Gmma(nGamma)
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, Index_A, Index_AB, Index_ABCD, Index_B, Index_C, Index_CD, Index_D, Indij, Indik, &
                     Indil, Indjk, Indjl, Indkl, iS, iShell_A, iShell_AB, iShell_B, iShell_C, iShell_CD, iShell_D, iSO_A, iSO_R, &
                     iSOi, iSOi_A, iSym(0:7), j, j1, j12, j123, j2, j3, j4, jAOj, jS, jSO_A, jSO_R, jSOj, jSOj_A, jSym(0:7), kAOk, &
                     kS, kSO_A, kSO_R, kSOk, kSOk_A, kSym(0:7), lAOl, lOper, lS, lSO_A, lSO_R, lSOl, lSOl_A, lSym(0:7), MemSO2, &
                     mijkl, nDim_A, nDim_AB, nDim_B, nDim_C, nDim_CD, nDim_D, niSym, njSym, nkSym, nlSym
real(kind=wp) :: T14, Temp
real(kind=wp), parameter :: exfac = One
integer(kind=iwp), external :: iPntSO
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iComp
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'nSOs=',nSOs
write(u6,*) 'iSO2Sh=',iSO2Sh
iComp = 1
call PrMtrx(' In PGet2:DSO ',[iD0Lbl],iComp,1,D0)
call PrMtrx(' In PGet2:DSO_Var ',[iD0Lbl],iComp,1,DVar)
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
                iSO_r = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO_r = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO_r = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO_r = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
                iSO_a = iSO_r+iOffSO(j1)
                jSO_a = jSO_r+iOffSO(j2)
                kSO_a = kSO_r+iOffSO(j3)
                lSO_a = lSO_r+iOffSO(j4)

                mijkl = 0
                do lAOl=0,lBas-1
                  lSOl = lSO_r+lAOl
                  lSOl_a = lSO_a+lAOl
                  iShell_D = iSO2Sh(lSOl_a)
                  Index_D = iSO2cI(1,lSOl_a)
                  nDim_D = iSO2cI(2,lSOl_a)
                  do kAOk=0,kBas-1
                    kSOk = kSO_r+kAOk
                    kSOk_a = kSO_a+kAOk
                    iShell_C = iSO2Sh(kSOk_a)
                    Index_C = iSO2cI(1,kSOk_a)
                    nDim_C = iSO2cI(2,kSOk_a)
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
                      jSOj = jSO_r+jAOj
                      jSOj_a = jSO_a+jAOj
                      iShell_B = iSO2Sh(jSOj_a)
                      Index_B = iSO2cI(1,jSOj_a)
                      nDim_B = iSO2cI(2,jSOj_a)
                      do iAOi=0,iBas-1
                        iSOi = iSO_r+iAOi
                        iSOi_a = iSO_a+iAOi
                        iShell_A = iSO2Sh(iSOi_a)
                        Index_A = iSO2cI(1,iSOi_a)
                        nDim_A = iSO2cI(2,iSOi_a)
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
                        mijkl = mijkl+1

                        !*********** columbus interface ****************************************
                        !do not reconstruct the two-particle density from the one-particle
                        !density or partial two-particle densities but simply read them from
                        !file

                        if (gamma_mrcisd) then
                          temp = Gmma(Index_ABCD)
                        else
                          ! Contribution D(ij)*D(kl) to P(ijkl)
                          if (j1 == j2) then
                            ! j3 == j4 also
                            Indij = iPntSO(j1,j2,lOper,nbas)+iTri(iSOi,jSOj)
                            Indkl = iPntSO(j3,j4,lOper,nbas)+iTri(kSOk,lSOl)
                            temp = DSO(Indij)*DSO(Indkl)+(DSO_Var(Indij)-DSO(Indij))*DSO(Indkl)+ &
                                   DSO(Indij)*(DSO_Var(Indkl)-DSO(Indkl))
                          else
                            temp = Zero
                          end if

                          ! Contribution -1/4*D(ik)*D(jl) to P(ijkl)
                          if (j1 == j3) then
                            ! j2 == j4 also
                            Indik = iPntSO(j1,j3,lOper,nbas)+iTri(iSOi,kSOk)
                            Indjl = iPntSO(j2,j4,lOper,nbas)+iTri(jSOj,lSOl)
                            temp = temp-t14*(DSO(Indik)*DSO(Indjl)+(DSO_Var(Indik)-DSO(Indik))*DSO(Indjl)+ &
                                             DSO(Indik)*(DSO_Var(Indjl)-DSO(Indjl))+DSSO(Indik)*DSSO(Indjl)+ &
                                             (DSSO_Var(Indik)-DSSO(Indik))*DSSO(Indjl)+DSSO(Indik)*(DSSO_Var(Indjl)-DSSO(Indjl)))
                          end if

                          ! Contribution -1/4*D(il)*D(jk) to P(ijkl)
                          if (j1 == j4) then
                            ! j2 == j3 also
                            Indil = iPntSO(j1,j4,lOper,nbas)+iTri(iSOi,lSOl)
                            Indjk = iPntSO(j2,j3,lOper,nbas)+iTri(jSOj,kSOk)
                            temp = temp-t14*(DSO(Indil)*DSO(Indjk)+(DSO_Var(Indil)-DSO(Indil))*DSO(Indjk)+ &
                                             DSO(Indil)*(DSO_Var(Indjk)-DSO(Indjk))+DSSO(Indil)*DSSO(Indjk)+ &
                                             (DSSO_Var(Indil)-DSSO(Indil))*DSSO(Indjk)+DSSO(Indil)*(DSSO_Var(Indjk)-DSSO(Indjk)))
                          end if

                          !write(u6,*) 'iSO:',iSOi_a,jSOj_a,kSOk_a,lSOl_a
                          !write(u6,*) 'iShell:',iShell_A,iShell_B,iShell_C,iShell_D
                          !write(u6,*) 'nDim:',nDim_A,nDim_B,nDim_C,nDim_D
                          !write(u6,*)  temp , Gmma(Index_ABCD),Index_ABCD
                          temp = temp+Four*Gmma(Index_ABCD)
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
  call WarningMessage(2,'PGet2_Aces: nPSO /= MemSO2')
  write(u6,*) nPSO,MemSO2
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet2:PSO ',' ',PSO,nijkl,nPSO)
#endif

return

end subroutine PGet2_Aces
