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

subroutine PGet1_CD2(PAO,ijkl,nPAO,iCmp,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,ExFac,CoulFac,PMax,V_k,U_k,mV_k,Z_p_K,nnP1)
!***********************************************************************
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
!                                                                      *
!             Modified for Cholesky 1-center gradients May 2007 by     *
!             R. Lindh                                                 *
!                                                                      *
!             Modified for RI-HF/CAS, Dec 2009 (F. Aquilante)          *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use ExTerm, only: CijK, iMP2prpt, nAuxVe

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "exterm.fh"
real*8 PAO(ijkl,nPAO), V_k(mV_k), U_K(mV_K), Z_p_K(nnP1,mV_K), Fac_ij, Fac_kl
integer iAO(4), kOp(4), iAOst(4), iCmp(4)
logical Shijij
external mn2K
real*8, pointer :: CiKj(:,:) => null(), V2(:) => null()

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

                  Indi = max(iSOi,jSOj)
                  Indj = iSOi+jSOj-Indi
                  Indk = max(kSOk,lSOl)
                  Indl = kSOk+lSOl-Indk
                  Indij = (Indi-1)*Indi/2+Indj
                  Indkl = (Indk-1)*Indk/2+Indl
                  temp = V_k(Indij)*V_k(Indkl)*CoulFac
                  ! Active space contribution (any factor?)
                  ijVec = mn2K(Indij,1)
                  klVec = mn2K(Indkl,1)
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
              Indk = max(kSOk,lSOl)
              Indl = kSOk+lSOl-Indk
              Indkl = (Indk-1)*Indk/2+Indl
              klVec = mn2K(Indkl,1)

              if (klvec /= 0) then
                iAdrL = NumIK*(klVec-1)+iAdrCVec(jSym,iSym,iSO)
                call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,1),NumIK,iAdrL)
              end if

              do jAOj=0,jBas-1
                jSOj = jSO+jAOj
                do iAOi=0,iBas-1
                  iSOi = iSO+iAOi
                  nijkl = nijkl+1

                  Indi = max(iSOi,jSOj)
                  Indj = iSOi+jSOj-Indi
                  Indij = (Indi-1)*Indi/2+Indj
                  ijVec = mn2K(Indij,1)

                  if ((ijVec /= klVec) .and. (ijvec /= 0)) then
                    iAdrJ = NumIK*(ijVec-1)+iAdrCVec(jSym,iSym,iSO)
                    call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,2),NumIK,iAdrJ)
                    V2(1:) => CiKj(1:,2)
                  else
                    V2(1:) => CiKj(1:,1)
                  end if

                  temp = V_k(Indij)*V_k(Indkl)*CoulFac

                  if ((ijVec /= 0) .and. (klVec /= 0)) then
                    if (Indi == Indj) then
                      Fac_ij = 1.0d0
                    else
                      Fac_ij = 0.5d0
                    end if
                    if (Indk == Indl) then
                      Fac_kl = 1.0d0
                    else
                      Fac_kl = 0.5d0
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
              Indk = max(kSOk,lSOl)
              Indl = kSOk+lSOl-Indk
              Indkl = (Indk-1)*Indk/2+Indl
              klVec = mn2K(Indkl,1)
              if (klvec /= 0) then
                iAdrL = NumIK*(klVec-1)+iAdrCVec(jSym,iSym,iSO)
                call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,1),NumIK,iAdrL)
              end if

              do jAOj=0,jBas-1
                jSOj = jSO+jAOj
                do iAOi=0,iBas-1
                  iSOi = iSO+iAOi
                  nijkl = nijkl+1

                  Indi = max(iSOi,jSOj)
                  Indj = iSOi+jSOj-Indi
                  Indij = (Indi-1)*Indi/2+Indj
                  ijVec = mn2K(Indij,1)

                  if ((ijVec /= klVec) .and. (ijvec /= 0)) then
                    iAdrJ = NumIK*(ijVec-1)+iAdrCVec(jSym,iSym,iSO)
                    call dDaFile(LuCVector(jSym,iSO),2,CiKj(:,2),NumIK,iAdrJ)
                    V2(1:) => CiKj(:,2)
                  else
                    V2(1:) => CiKj(:,1)
                  end if

                  temp = V_k(Indij)*V_k(Indkl)*CoulFac+V_K(Indij)*U_K(Indkl)*CoulFac+U_K(Indij)*V_K(Indkl)*CoulFac

                  if ((ijVec /= 0) .and. (klVec /= 0)) then
                    if (Indi == Indj) then
                      Fac_ij = 1.0d0
                    else
                      Fac_ij = 0.5d0
                    end if
                    if (Indk == Indl) then
                      Fac_kl = 1.0d0
                    else
                      Fac_kl = 0.5d0
                    end if

                    ! MP2 contribution
                    call Compute_A_jk_Mp2(1,ijVec,klVec,tempJ_mp2,Fac_ij,Fac_kl,nAuxVe,2)
                    temp = temp+tempJ_mp2*CoulFac

                    ! Exchange contribution

                    tempK = 2.0d0*Fac_ij*Fac_kl*dDot_(NumIK,CiKj(:,1),1,V2,1)
                    call compute_A_jk_Mp2(1,ijVec,klVec,tempK_mp2,fac_ij,fac_kl,nAuxVe,1)

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

CiKj => null()
V2 => null()

if (iPAO /= nPAO) then
  write(6,*) ' Error in PGet1_CD2!'
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
! Avoid unused argument warnings
if (.false.) then
  call Unused_logical(Shijij)
end if

end subroutine PGet1_CD2
