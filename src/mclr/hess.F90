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
use MCLR_Data, only: Hss, CMO, F0SQMO
use MCLR_Data, only: nDens2, ipCM, ipMat
use MCLR_Data, only: DspVec, lDisp
use input_mclr, only: nSym, nBas, nOrb, nTPert
use Constants, only: Zero, One, Two, Half
use Definitions, only: u6

implicit none
real*8 Temp1(nDens2), Temp2(nDens2), Temp3(nDens2), FockC(nDens2), FockX(nDens2), rcon(nDens2), temp4(*)
integer idSym, jDisp, iDisp
character(len=8) Label
integer iS, jS, nNJ, Len, iSym, nIn, mDisp, kDisp, iRC, iOpt, iOp, iP, IndX
real*8 Fact
real*8, external :: DDot_

Temp3(:) = Zero
if (iand(ntpert(idisp),2**3) == 8) then
  do iS=1,nSym
    js = ieor(is-1,idSym-1)+1
    nnj = nOrb(js)!nash(js)+nash(js)
    if (nOrb(is)*nOrb(js) /= 0) call DGEMM_('N','N',nOrb(is),nnj,nnj,One,rCon(ipMat(is,js)),nOrb(is),F0SQMO(ipCM(jS)),nOrb(js), &
                                            Zero,Temp3(ipMat(is,js)),nOrb(is))

  end do
  !Temp3(:) = -Half*Temp3(:)+Half*FockC(:)+FockX(:)
  call DScal_(ndens2,-Half,Temp3,1)
  call DaXpY_(nDens2,Half,FockC,1,Temp3,1)
  call DaXpY_(nDens2,One,FockX,1,Temp3,1)
else
  Temp3(:) = FockX(:)
end if

!              xa     ca      xa
! Temp3 = Y = F + 1/2F  -1/2 S  F

Len = 0
do iSym=1,nSym
  Len = Len+nTri_Elem(lDisp(iSym))
end do

nIn = 0
mdisp = 0
do iS=1,idsym-1
  mdisp = mdisp+ldisp(is)
  nIn = nIn+nTri_Elem(lDisp(is))
end do

do kDisp=1,ldisp(idsym)
  mDisp = mdisp+1
  if (iand(ntpert(mdisp),2**3) == 0) cycle
  iRC = -1
  iOpt = 0
  Label = 'OvrGrd'
  iOp = 2**idSym
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
        if (ieor(iS-1,jS-1) == idsym-1) then
          if (is == js) then
            call Square(Temp1(ip),Temp2(ipMat(iS,jS)),1,nBas(is),nBas(is))
            ip = ip+nTri_Elem(nBas(is))
          else
            if (nBas(is)*nBas(js) /= 0) call dcopy_(nBas(iS)*nBas(jS),Temp1(ip),1,Temp2(ipMat(iS,jS)),1)
            ip = ip+nBas(is)*nBas(js)
          end if
          if (nBas(is)*nBas(js) /= 0) then
            call DGEMM_('T','N',nOrb(iS),nBAs(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),Temp2(ipMat(iS,jS)),nBas(iS),Zero,Temp4, &
                        nOrb(is))
            call dcopy_(nBas(is)*nBas(js),[Zero],0,temp2(ipMat(is,js)),1)
            call DGEMM_('N','N',nOrb(iS),nOrb(jS),nBas(jS),One,Temp4,nOrb(iS),CMO(ipCM(jS)),nBas(jS),Zero,Temp2(ipMat(iS,jS)), &
                        nOrb(iS))
            if (is /= js) then
              call dcopy_(nBas(is)*nBas(js),[Zero],0,temp2(ipMat(js,is)),1)
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
  Hss(Indx) = Hss(Indx)-fact*ddot_(nDens2,Temp2,1,Temp3,1)
end do

end subroutine Hess
