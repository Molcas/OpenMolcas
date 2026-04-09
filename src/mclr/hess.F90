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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine Hess(FockC,FockX,rCon,Temp1,Temp2,Temp3,Temp4,idsym,jdisp,idisp)
! Constructs the connection parts that is dependend on the first
! derivative of the connection.

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, DspVec, F0SQMO, Hss, ipCM, ipMat, lDisp, nDens
use input_mclr, only: nBas, nOrb, nSym, nTPert
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: FockC(nDens), FockX(nDens), rcon(nDens)
real(kind=wp), intent(out) :: Temp1(nDens), Temp2(nDens), Temp3(nDens)
real(kind=wp), intent(_OUT_) :: temp4(*)
integer(kind=iwp), intent(in) :: idSym, jDisp, iDisp
integer(kind=iwp) :: IndX, iOp, iOpt, iP, iRC, iS, jS, kDisp, mDisp, nIn, nNJ
real(kind=wp) :: Fact
character(len=8) :: Label
real(kind=wp), external :: DDot_

Temp3(:) = FockX(:)
if (btest(ntpert(idisp),3)) then
  do iS=1,nSym
    js = Mul(is,idSym)
    nnj = nOrb(js)!nash(js)+nash(js)
    if (nOrb(is)*nOrb(js) /= 0) call DGEMM_('N','N',nOrb(is),nnj,nnj,-Half,rCon(ipMat(is,js)),nOrb(is),F0SQMO(ipCM(jS)),nOrb(js), &
                                            One,Temp3(ipMat(is,js)),nOrb(is))

  end do
  Temp3(:) = Temp3(:)+Half*FockC(:)
end if

!              xa     ca      xa
! Temp3 = Y = F + 1/2F  -1/2 S  F

nIn = 0
mdisp = sum(lDisp(1:idsym-1))
do iS=1,idsym-1
  nIn = nIn+nTri_Elem(lDisp(is))
end do

do kDisp=1,ldisp(idsym)
  mDisp = mdisp+1
  if (.not. btest(ntpert(mdisp),3)) cycle
  iRC = -1
  iOpt = 0
  Label = 'OvrGrd'
  iOp = ibset(0,idSym)
  call dRdMck(iRC,iOpt,Label,DspVec(mDisp),Temp1,iop)
  if (iRc /= 0) then
    write(u6,*) 'Hess: Error reading MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  ip = 1
  do iS=1,nSym
    do jS=1,iS
      if (nOrb(is)*nOrb(jS) /= 0) then
        if (Mul(iS,jS) == idsym) then
          if (is == js) then
            call Square(Temp1(ip),Temp2(ipMat(iS,jS)),1,nBas(is),nBas(is))
            ip = ip+nTri_Elem(nBas(is))
          else
            if (nBas(is)*nBas(js) /= 0) Temp2(ipMat(iS,jS):ipMat(iS,jS)+nBas(iS)*nBas(jS)-1) = Temp1(ip:ip+nBas(iS)*nBas(jS)-1)
            ip = ip+nBas(is)*nBas(js)
          end if
          if (nBas(is)*nBas(js) /= 0) then
            call DGEMM_('T','N',nOrb(iS),nBAs(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),Temp2(ipMat(iS,jS)),nBas(iS),Zero,Temp4, &
                        nOrb(is))
            Temp2(ipMat(is,js):ipMat(is,js)+nBas(is)*nBas(js)-1) = Zero
            call DGEMM_('N','N',nOrb(iS),nOrb(jS),nBas(jS),One,Temp4,nOrb(iS),CMO(ipCM(jS)),nBas(jS),Zero,Temp2(ipMat(iS,jS)), &
                        nOrb(iS))
            if (is /= js) then
              Temp2(ipMat(js,is):ipMat(js,is)+nBas(is)*nBas(js)-1) = Zero
              call DGEMM_('T','T',nOrb(js),nOrb(iS),nBas(js),One,CMO(ipCM(js)),nBas(js),Temp4,nOrb(is),Zero,Temp2(ipMat(js,is)), &
                          nOrb(js))
            end if
          end if
        end if
      end if
    end do
  end do
  Fact = One
  if (kDisp == jDisp) Fact = Two
  Indx = nIn+iTri(kDisp,jDisp)
  Hss(Indx) = Hss(Indx)-fact*ddot_(nDens,Temp2,1,Temp3,1)
end do

end subroutine Hess
