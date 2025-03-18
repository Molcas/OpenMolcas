!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine FockGen_sp(d_0,rDens1,rdens2,Fock,fockout,idsym)
!***********************************************************************
!                                                                      *
!   Constructs active fockmatrix and Q matrix                          *
!                                                                      *
!   Input: rkappa : Rotation matrix                                    *
!          idsym  : symmetry of perturbation                           *
!                                                                      *
!                                                                      *
!   Output:MO     :MO integrals                                        *
!          Fock   :Fock matrix (one index transformed integrals)       *
!          MOtilde:MO (one index transformed integrals)                *
!                                                                      *
!***********************************************************************

use Arrays, only: FIMO
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nNA, ipCM, ipMat, nA, nDens2
use input_mclr, only: nSym, nAsh, nIsh, nBas

implicit none
real*8 d_0
integer idSym
real*8 Fock(*), fockout(*), rdens2(*), rDens1(*)
!real*8 Fock(nDens2), fockout(*), rdens2(*), rDens1(nna*nna)
real*8, allocatable :: MO(:), Scr(:)
integer n1, iS, n2, ipS, kS, jS, iB, jA, jAA, kA, kAA, ipM, ipF, iA, ip1, ip2
real*8 rd

!                                                                      *
!***********************************************************************
!                                                                      *
!  Coulomb term: F  = 2(pk|ji)d
!                 kp           ij

call dcopy_(nDens2,[0.0d0],0,Fock,1)

n1 = 0
do iS=1,nSym
  n1 = max(n1,nBas(iS))
end do
n2 = n1**2
call mma_allocate(MO,n2,Label='MO')
call mma_allocate(Scr,n2,Label='Scr')
!
do ips=1,nSym
  do ks=1,nSym
    do is=1,nSym
      jS = ieor(ieor(ips-1,ks-1),is-1)+1

      ! Exchange term   F = -(pk|ji)d     (i>j)
      !                  pl          kj

      if ((ieor(ips-1,iS-1)+1 == iDsym) .and. (nBas(ipS) > 0)) then
        do iB=1,nIsh(iS)
          do jA=1,nAsh(jS)
            jAA = jA+nIsh(js)

            call Coul(ipS,kS,iS,jS,iB,jAA,MO,Scr)

            do kA=1,nAsh(ks)
              kAA = kA+nIsh(kS)

              ipM = 1+(kAA-1)*nBas(ipS)
              ipF = ipMat(ipS,iS)+nBas(ipS)*(iB-1)
              rd = rDens1(jA+nA(jS)+(kA+nA(ks)-1)*nna)
              call DaXpY_(nBas(ipS),-rd,MO(ipM),1,Fock(ipF),1)
            end do
          end do
        end do
      end if

    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  if (nBas(iS) > 0) then
    jS = ieor(is-1,iDSym-1)+1
    do iA=1,nAsh(is)
      do jA=1,nAsh(js)
        rd = rDens1(iA+nA(iS)+(jA+nA(js)-1)*nna)
        ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        call DaxPy_(nBAs(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
      end do
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call CreQADD_sp(Fock,rdens2,idsym,MO,Scr,n2)

call mma_deallocate(Scr)
call mma_deallocate(MO)

do iS=1,nSym
  js = ieor(is-1,idsym-1)+1
  if (nbas(is)*nBas(js) /= 0) &
    call DGESUB(Fock(ipMat(is,js)),nBas(is),'N',Fock(ipMat(js,is)),nBas(js),'T',FockOut(ipMat(is,js)),nBas(is),nBas(is),nBas(js))
end do
call DScal_(ndens2,2.0d0,FockOut,1)
if (idsym == 1) call Add2(Fockout,d_0)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FockGen_sp
