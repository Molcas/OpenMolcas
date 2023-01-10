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
!               2009, Francesco Aquilante                              *
!***********************************************************************

subroutine PGet1_CD2(PAO,ijkl,nPAO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,ExFac,CoulFac,PMax,V_k,U_k,mV_k,Z_p_K,nnP1)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density.                        *
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
!                                                                      *
!             Modified for RI-HF/CAS, Dec 2009 (F. Aquilante)          *
!***********************************************************************

use Index_Functions, only: iTri
use pso_stuff, only: ij2K, iOff_ij2K
use SOAO_Info, only: iAOtSO
use RI_glob, only: CijK, iAdrCVec, iMP2prpt, LuCVector, nAuxVe, nIJ1, tavec
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4), mV_k, nnP1
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PMax
real(kind=wp), intent(in) :: ExFac, CoulFac, V_k(mV_k), U_K(mV_K), Z_p_K(nnP1,mV_K)
integer(kind=iwp) :: i1, i2, i3, i4, iAdrJ, iAdrL, iAOi, iE, ijVec, Indij, Indkl, iPAO, iS, iSO, iSOi, iSym, jAOj, jp, jSO, jSOj, &
                     jSym, kAOk, klVec, kSO, kSOk, lAOl, lSO, lSOl, lSym, nijkl, NumIK
real(kind=wp) :: Cpu, Cpu1, Cpu2, Fac, Fac_ij, Fac_kl, temp, tempJ_mp2, tempK, tempK_mp2, Wall, Wall1, Wall2
real(kind=wp), pointer :: CiKj(:,:), V2(:)
real(kind=wp), external :: dDot_

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('PGet1_CD2: V_k',' ',V_k,1,mV_k)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Quadruple loop over elements of the basis functions angular description.
! Observe that we will walk through the memory in PAO in a sequential way.
!
!Fac = Quart

call CWTime(Cpu1,Wall1)

Fac = One
PMax = Zero
iPAO = 0

iSym = 1
jSym = 1
lSym = 1
iSO = 1

if (ExFac == Zero) then

  do i1=1,iCmp(1)
    do i2=1,iCmp(2)
      do i3=1,iCmp(3)
        do i4=1,iCmp(4)

          ! Unfold the way the eight indicies have been reordered.
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

                  ! V_k(ij)*V_k(kl)

                  Indij = iTri(iSOi,jSOj)
                  Indkl = iTri(kSOk,lSOl)
                  temp = V_k(Indij)*V_k(Indkl)*CoulFac
                  ! Active space contribution (any factor?)
                  ijVec = ij2K(iOff_ij2K(1)+Indij)
                  klVec = ij2K(iOff_ij2K(1)+Indkl)
                  if ((ijVec /= 0) .and. (klVec /= 0)) then
                    do jp=1,nnP1
                      temp = temp+Z_p_K(jp,ijVec)*Z_p_K(jp,klVec)
                    end do
                  end if

                  PMax = max(PMax,abs(temp))
                  PAO(nijkl,iPAO) = Fac*temp

                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
else if (iMP2prpt /= 2) then
  NumIK = nIJ1(iSym,lSym,iSO)
  if (NumIK == 0) return

  iS = 1
  iE = NumIK*2
  CiKj(1:NumIK,1:2) => CijK(iS:iE)

  do i1=1,iCmp(1)
    do i2=1,iCmp(2)
      iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
      jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

      do i3=1,iCmp(3)
        do i4=1,iCmp(4)
          kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
          lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

          iPAO = iPAO+1
          nijkl = 0

          do lAOl=0,lBas-1
            lSOl = lSO+lAOl
            do kAOk=0,kBas-1
              kSOk = kSO+kAOk
              Indkl = iTri(kSOk,lSOl)
              klVec = ij2K(iOff_ij2K(1)+Indkl)

              if (klvec /= 0) then
                iAdrL = NumIK*(klVec-1)+iAdrCVec(jSym,iSym,iSO)
                call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,1),NumIK,iAdrL)
              end if

              do jAOj=0,jBas-1
                jSOj = jSO+jAOj
                do iAOi=0,iBas-1
                  iSOi = iSO+iAOi
                  nijkl = nijkl+1

                  Indij = iTri(iSOi,jSOj)
                  ijVec = ij2K(iOff_ij2K(1)+Indij)

                  if ((ijVec /= klVec) .and. (ijvec /= 0)) then
                    iAdrJ = NumIK*(ijVec-1)+iAdrCVec(jSym,iSym,iSO)
                    call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,2),NumIK,iAdrJ)
                    V2(1:) => CiKj(1:,2)
                  else
                    V2(1:) => CiKj(1:,1)
                  end if

                  temp = V_k(Indij)*V_k(Indkl)*CoulFac

                  if ((ijVec /= 0) .and. (klVec /= 0)) then
                    if (iSOi == jSOj) then
                      Fac_ij = One
                    else
                      Fac_ij = Half
                    end if
                    if (kSOk == lSOl) then
                      Fac_kl = One
                    else
                      Fac_kl = Half
                    end if

                    ! Exchange contribution

                    temp = temp-ExFac*Fac_ij*Fac_kl*dDot_(NumIK,CiKJ(:,1),1,V2,1)
                    ! Active space contribution (any factor?)
                    do jp=1,nnP1
                      temp = temp+Z_p_K(jp,ijVec)*Z_p_K(jp,klVec)
                    end do
                  end if

                  PMax = max(PMax,abs(temp))
                  PAO(nijkl,iPAO) = Fac*temp

                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
else
  NumIK = nIJ1(iSym,lSym,iSO)
  if (NumIK == 0) return

  iS = 1
  iE = NumIK*2
  CiKj(1:NumIK,1:2) => CijK(iS:iE)

  do i1=1,iCmp(1)
    do i2=1,iCmp(2)
      iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
      jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
      do i3=1,iCmp(3)
        do i4=1,iCmp(4)
          kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
          lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

          iPAO = iPAO+1
          nijkl = 0

          do lAOl=0,lBas-1
            lSOl = lSO+lAOl
            do kAOk=0,kBas-1
              kSOk = kSO+kAOk
              Indkl = iTri(kSOk,lSOl)
              klVec = ij2K(iOff_ij2K(1)+Indkl)
              if (klvec /= 0) then
                iAdrL = NumIK*(klVec-1)+iAdrCVec(jSym,iSym,iSO)
                call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,1),NumIK,iAdrL)
              end if

              do jAOj=0,jBas-1
                jSOj = jSO+jAOj
                do iAOi=0,iBas-1
                  iSOi = iSO+iAOi
                  nijkl = nijkl+1

                  Indij = iTri(iSOi,jSOj)
                  ijVec = ij2K(iOff_ij2K(1)+Indij)

                  if ((ijVec /= klVec) .and. (ijvec /= 0)) then
                    iAdrJ = NumIK*(ijVec-1)+iAdrCVec(jSym,iSym,iSO)
                    call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,2),NumIK,iAdrJ)
                    V2(1:) => CiKj(:,2)
                  else
                    V2(1:) => CiKj(:,1)
                  end if

                  temp = V_k(Indij)*V_k(Indkl)*CoulFac+V_K(Indij)*U_K(Indkl)*CoulFac+U_K(Indij)*V_K(Indkl)*CoulFac

                  if ((ijVec /= 0) .and. (klVec /= 0)) then
                    if (iSOi == jSOj) then
                      Fac_ij = One
                    else
                      Fac_ij = Half
                    end if
                    if (kSOk == lSOl) then
                      Fac_kl = One
                    else
                      Fac_kl = Half
                    end if

                    ! MP2 contribution
                    call Compute_A_jk_Mp2(ijVec,klVec,tempJ_mp2,Fac_ij,Fac_kl,nAuxVe,2)
                    temp = temp+tempJ_mp2*CoulFac

                    ! Exchange contribution

                    tempK = Two*Fac_ij*Fac_kl*dDot_(NumIK,CiKj(:,1),1,V2,1)
                    call compute_A_jk_Mp2(ijVec,klVec,tempK_mp2,fac_ij,fac_kl,nAuxVe,1)

                    tempK = tempK+tempK_mp2

                    temp = temp-ExFac*tempK*Half
                    ! Active space contribution (any factor?)
                    do jp=1,nnP1
                      temp = temp+Z_p_K(jp,ijVec)*Z_p_K(jp,klVec)
                    end do
                  end if

                  PMax = max(PMax,abs(temp))
                  PAO(nijkl,iPAO) = Fac*temp

                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end if

nullify(CiKj,V2)

if (iPAO /= nPAO) then
  write(u6,*) ' Error in PGet1_CD2!'
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet1_CD2:PAO ',' ',PAO,ijkl,nPAO)
#endif

call CWTime(Cpu2,Wall2)
Cpu = Cpu2-Cpu1
Wall = Wall2-Wall1
tavec(1) = tavec(1)+Cpu
tavec(2) = tavec(2)+Wall

return

end subroutine PGet1_CD2
