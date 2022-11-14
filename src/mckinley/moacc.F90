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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine MOAcc(AOInt,Temp1,Temp2,nTemp,ishell,Ck,nCk,Cl,nCl,moip,nACO,pert,nOp,iBasa,iCmpa,icar,icnt,indgrd,fact,iaost,Buffer, &
                 nij,nkl,nbasi,nbasj,icmp,jcmp)
!***********************************************************************
!                                                                      *
!     Transforms a batch of unsymmetrized integrals to                 *
!     active integral batches and FM                                   *
!     All MO combinations are constructed                              *
!     They will be needed with unsymmetric perurbations                *
!                                                                      *
!     Author: Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
!             University of Lund, Sweden. Januar '96                   *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: iChTbl, iOper, nIrrep, Prmt
use Gateway_Info, only: CutInt
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTemp, ishell(4), nCk, nCl, moip(0:7), nACO, nOp(4), ibasa(4), icmpa(4), icar, icnt, &
                                 indgrd(3,4,0:7), iAOST(4), nij, nkl, nbasi, nbasj, icmp, jcmp
real(kind=wp), intent(in) :: AOInt(nkl,nij), fact
real(kind=wp), intent(out) :: Temp1(nTemp), Temp2(naco,naco)
real(kind=wp), intent(inout) :: Ck(nCk), Cl(nCl), Buffer(nbasi,icmp,nbasj,jcmp,0:nirrep-1,nTri_Elem(naco),*)
logical(kind=iwp), intent(in) :: pert(0:7)
#include "etwas.fh"
integer(kind=iwp) :: i, ib, iBas, ic, iCB, iirr, ij, il, ipC, ipM, irest, iSPert, j, jb, jBas, jc, jIrr, k, kAsh, kBas, kCmp, &
                     kIrr, kIrrep, kk, kMax, l, lAsh, lBas, lCmp, lIrr, ll, nt
real(kind=wp) :: rFact, rFact2, rk, rl, rPij, rPj, sfact, vij
integer(kind=iwp), external :: NrOpr
real(kind=wp), external :: DNrm2_

iCB = 2**(icar-1)
rFact = Prmt(ioper(nOp(icnt)),icb)*fact

iBas = iBasa(1)
jBas = iBasa(2)
kBas = iBasa(3)
lBas = iBasa(4)

kCmp = iCmpa(3)
lCmp = iCmpa(4)
kk = 0
do kIrrep=0,nIrrep-1
  sfact = real(ichtbl(kirrep,nop(3)),kind=wp)
  Ck(kk+1:kk+nAsh(kIrrep)*kcmp*kbas) = Ck(kk+1:kk+nAsh(kIrrep)*kcmp*kbas)*sFact
  kk = kk+nAsh(kIrrep)*kcmp*kbas
end do
kk = 0
do kIrrep=0,nIrrep-1
  sfact = real(ichtbl(kirrep,nop(4)),kind=wp)
  Cl(kk+1:kk+nAsh(kIrrep)*lcmp*lbas) = Cl(kk+1:kk+nAsh(kIrrep)*lcmp*lbas)*sFact
  kk = kk+nAsh(kIrrep)*lcmp*lbas
end do
rk = DNrm2_(nck,ck,1)
rl = DNrm2_(ncl,cl,1)
ij = 0
nt = lBas*lCmp*kcmp*kbas
do jc=1,jcmp
  do jb=1+iaost(2),iaost(2)+jbas
    do ic=1,icmp
      do ib=1+iaost(1),iaost(1)+ibas

        ij = ij+1
        vij = DNrm2_(nt,AOInt(1,ij),1)
        if (abs(vij*rk*rl) < cutint) cycle
        ipC = 0
        Temp1(1:nAco*lbas*lcmp) = Zero
        do kAsh=1,nAco
          ipM = (kAsh-1)*lbas*lcmp
          il = 0
          do i=1,lbas*lcmp
            do k=1,kCmp*kBas
              il = il+1
              Temp1(ipM+i) = Temp1(ipM+i)+Ck(ipC+k)*AOINT(il,ij)
            end do
          end do
          ipC = ipC+kBas*kCmp
        end do
        ipC = 0
        Temp2(:,:) = Zero
        do lAsh=1,naco
          il = 0
          do kash=1,naco
            do l=1,lbas*lcmp
              il = il+1
              Temp2(kash,lash) = Temp2(kash,lash)+Cl(ipC+l)*Temp1(il)
            end do
          end do
          ipC = ipC+lBas*lCmp
        end do

        if (iShell(3) /= ishell(4)) then

          do iSPert=0,nIrrep-1
            if (pert(isPert)) then
              rFact2 = rFact*real(iChtbl(ispert,nop(icnt)),kind=wp)
              k = abs(indgrd(icar,icnt,ispert))
              j = 0
              do lIrr=0,nIrrep-1
                do lAsh=1,nAsh(lIrr)
                  do kIrr=0,lIrr
                    irest = ieor(ieor(ioper(ispert),ioper(kirr)),ioper(lirr))
                    kMax = nAsh(kIrr)
                    if (kIrr == lIrr) kMax = lAsh
                    do kAsh=1,kMax
                      kk = kash+moip(kirr)
                      ll = lash+moip(lirr)
                      j = j+1
                      do jIrr=0,nIrrep-1
                        rPj = real(iChTbl(jIrr,nop(2)),kind=wp)
                        iirr = nropr(ieor(iOPER(jirr),irest))
                        rPij = rPj*real(iChTbl(iIrr,nop(1)),kind=wp)*rfact2
                        buffer(ib,ic,jb,jc,iirr,j,k) = buffer(ib,ic,jb,jc,iirr,j,k)+rpij*Temp2(kk,ll)+rpij*Temp2(ll,kk)
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end do
        else

          do iSPert=0,nIrrep-1
            if (pert(isPert)) then
              rFact2 = rFact*real(iChtbl(ispert,nop(icnt)),kind=wp)
              k = abs(indgrd(icar,icnt,ispert))
              j = 0
              do lIrr=0,nIrrep-1
                do lAsh=1,nAsh(lIrr)
                  do kIrr=0,lIrr
                    irest = ieor(ieor(ioper(ispert),ioper(kirr)),ioper(lirr))
                    kMax = nAsh(kIrr)
                    if (kIrr == lIrr) kMax = lAsh
                    do kAsh=1,kMax
                      kk = kash+moip(kirr)
                      ll = lash+moip(lirr)
                      j = j+1
                      do jIrr=0,nIrrep-1
                        rPj = real(iChTbl(jIrr,nop(2)),kind=wp)
                        iirr = nropr(ieor(iOPER(jirr),irest))
                        rPij = rPj*real(iChTbl(iIrr,nop(1)),kind=wp)*rfact2
                        buffer(ib,ic,jb,jc,iirr,j,k) = buffer(ib,ic,jb,jc,iirr,j,k)+rpij*Temp2(kk,ll)
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end do
        end if
      end do
    end do
  end do
end do

kk = 0
do kIrrep=0,nIrrep-1
  sfact = real(ichtbl(kirrep,nop(3)),kind=wp)
  Ck(kk+1:kk+nAsh(kIrrep)*kcmp*kbas) = Ck(kk+1:kk+nAsh(kIrrep)*kcmp*kbas)*sFact
  kk = kk+nAsh(kIrrep)*kcmp*kbas
end do
kk = 0
do kIrrep=0,nIrrep-1
  sfact = real(ichtbl(kirrep,nop(4)),kind=wp)
  Cl(kk+1:kk+nAsh(kIrrep)*lcmp*lbas) = Cl(kk+1:kk+nAsh(kIrrep)*lcmp*lbas)*sFact
  kk = kk+nAsh(kIrrep)*lcmp*lbas
end do

return

end subroutine MOAcc
