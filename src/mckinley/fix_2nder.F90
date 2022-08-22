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

subroutine Fix_2nder(Fa1,Fa2,Fb1,Fb2,rFinal,nalpha,nbeta,ishll,la,lb,iang,jfhess,nfun,fact)
! Fa1 includes  (  grad   < a | c > , < a | c >  )
! Fb1 includes  (  grad   < c | b > , < c | b >  )

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nalpha, nbeta, ishll, la, lb, iang, nfun
real(kind=wp) :: FA1((2*iAng+1)*((la+1)*(la+2)/2)*nAlpha*nFun*4), FA2((2*iAng+1)*((la+1)*(la+2)/2)*nAlpha*nFun*6), &
                 FB1((2*iAng+1)*((lb+1)*(lb+2)/2)*nBeta*nFun*4), FB2((2*iAng+1)*((lb+1)*(lb+2)/2)*nBeta*nFun*6), &
                 rFinal(nAlpha*nbeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,21), fact
!BS              rFinal(nAlpha*nbeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6)
logical(kind=iwp) :: jfhess(3,3,4,4)
!BS                  jfhess(3,3,2,2)  ! now same dimensions as in prjhss.f
integer(kind=iwp) :: ia, iaC, ib, iC, iCar, iCb, ipaC, ipCb, ipFA1a, ipFB1a, jCar, mVec
! Statement functions
integer(kind=iwp) :: nElem, ixyz, iTri, i, j
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

! OK NOW EVERYTHING SHOULD BE MERGED TOGETHER

!BS write(u6,*) 'nFun ',nfun
!BS write(u6,*) 'FA1 ',FA1
!BS write(u6,*) 'FA2 ',FA2
!BS write(u6,*) 'FB1 ',FB1
!BS write(u6,*) 'FB2 ',FB2
do iCar=1,3
  do jCar=1,3
    mvec = itri(iCar+3,jcar)
    if (jfHess(iCar,jcar,2,1)) then
      ipFA1a = 1+iCar*nAlpha*nFun*nElem(la)*(2*iAng+1)
      ipFB1a = 1+jcar*nFun*nBeta*(2*iAng+1)*nElem(lb)

      do ib=1,nElem(lb)
        do ia=1,nElem(la)

          do iC=1,(2*iAng+1)
            iaC = (iC-1)*nElem(la)+ia
            ipaC = (iaC-1)*nAlpha*nFun+ipFA1a
            iCb = (ib-1)*(2*iAng+1)+iC
            ipCb = (iCb-1)*nFun*nBeta+ipFB1a

            call DGEMM_('N','N',nAlpha,nBeta,nFun,Fact,FA1(ipaC),nAlpha,FB1(ipCb),nFun,One,rFinal(1,ia,ib,mVec),nAlpha)

          end do
        end do
      end do
    end if
  end do
end do

!BS if (.false.) then
do iCar=1,3
  do jCar=1,icar
    mvec = itri(iCar,jcar)
    if (jfHess(iCar,jcar,1,1)) then
      ipFA1a = 1+(itri(iCar,jCar)-1)*nAlpha*nFun*nElem(la)*(2*iAng+1)

      do ib=1,nElem(lb)
        do ia=1,nElem(la)

          do iC=1,(2*iAng+1)
            iaC = (iC-1)*nElem(la)+ia
            ipaC = (iaC-1)*nAlpha*nFun+ipFA1a
            iCb = (ib-1)*(2*iAng+1)+iC
            ipCb = (iCb-1)*nFun*nBeta+1

            call DGEMM_('N','N',nAlpha,nBeta,nFun,Fact,FA2(ipaC),nAlpha,FB1(ipCb),nFun,One,rFinal(1,ia,ib,mVec),nAlpha)

          end do
        end do
      end do
    end if
  end do
end do
!BS end if

!BS if (.false) then ! for testing the different contributions
do iCar=1,3
  do jCar=1,icar
    mvec = itri(3+icar,3+jcar)
    if (jfHess(iCar,jcar,2,2)) then
      ipFB1a = 1+(itri(icar,jcar)-1)*nBeta*nFun*nElem(lb)*(2*iAng+1)
      !BS nAlpha*nFun*nElem(la)*(2*iAng+1)
      !BS looks much better ....  (Oh Anders ...)

      do ib=1,nElem(lb)
        do ia=1,nElem(la)

          do iC=1,(2*iAng+1)
            iaC = (iC-1)*nElem(la)+ia
            ipaC = (iaC-1)*nAlpha*nFun+1
            iCb = (ib-1)*(2*iAng+1)+iC
            ipCb = (iCb-1)*nFun*nBeta+ipFB1a

            call DGEMM_('N','N',nAlpha,nBeta,nFun,Fact,FA1(ipaC),nAlpha,FB2(ipCb),nFun,One,rFinal(1,ia,ib,mVec),nAlpha)

          end do
        end do
      end do
    end if
  end do
end do
!BS end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(ishll)

end subroutine Fix_2nder
