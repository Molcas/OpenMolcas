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

subroutine m1kernel(rFinal,Hess,nHess,DAO,nDAO,iAng,nRys,nZeta,Alpha,Beta,Zeta,rKappa,P,TC,Coor,CoorAc,Array,nArray,ifgrd,indgrd, &
                    ifhss,indhss,ifg,tr,nop,iuvwx,kCnttp,fact,loper,idcar)

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc
use Symmetry_Info, only: nIrrep
use Constants, only: One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nHess, nDAO, iAng(4), nRys, nZeta, nArray, indgrd(3,4,0:7), indhss(3,4,3,4,0:7), nop(4), &
                                 iuvwx(4), kCnttp, loper, idcar
real(kind=wp), intent(inout) :: rFinal(*), Hess(nHess)
real(kind=wp), intent(in) :: DAO(nZeta,nDAO), Alpha(nZeta), Beta(nZeta), Zeta(nZeta), rKappa(nZeta), P(nZeta,3), TC(3), Coor(3,4), &
                             CoorAC(3,2), fact
real(kind=wp), intent(out) :: Array(nArray)
logical(kind=iwp), intent(in) :: ifgrd(3,4), ifhss(3,4,3,4), ifg(4)
logical(kind=iwp), intent(inout) :: tr(4)
integer(kind=iwp) :: iDAO, iElem, iM1xp, indi, Indx(3,4), ip, ipDAO, ipDAOt, ipK, ipPx, ipPy, ipPz, ipZ, ipZI, iZeta, &
                     jndgrd(3,4,0:7), jndhss(3,4,3,4,0:7), nb, nGr
real(kind=wp) :: coori(3,4), FactECP, Gmma, PTC2, Tmp0, Tmp1
logical(kind=iwp) :: jfg(4), jfgrd(3,4), jfhss(3,4,3,4), lGrad, lHess
logical(kind=iwp), external :: EQ
external :: Cff2D, Fake, TNAI1

lGrad = idcar /= 0
lHess = nHess /= 0
coori(:,:) = coor
if ((.not. EQ(coor(1,1),coor(1,2))) .or. (.not. EQ(coor(1,1),coor(1,3)))) then
  Coori(1,1) = Coori(1,1)+One
end if

ip = 1
ipK = ip
ip = ip+nZeta
ipZ = ip
ip = ip+nZeta
ipZI = ip
ip = ip+nZeta
ipPx = ip
ip = ip+nZeta
ipPy = ip
ip = ip+nZeta
ipPz = ip
ip = ip+nZeta
ipDAO = ip
ip = ip+nDAO*nZeta
if (ip >= narray) then
  write(u6,*) 'Out of memory in m1kernel (',narray,',',ip,')'
  call Abend()
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'nM1=',dbsc(kCnttp)%nM1,'kCnttp=',kCnttp
#endif

do iM1xp=1,dbsc(kCnttp)%nM1
  Gmma = dbsc(kCnttp)%M1xp(iM1xp)
  FactECP = dbsc(kCnttp)%M1cf(iM1xp)*Fact

# ifdef _DEBUGPRINT_
  write(u6,*) 'Fact=',FactECP,Fact
  write(u6,*) 'im1xp=',iM1xp
  write(u6,*) 'Gamma=',Gmma
# endif

  ! Modify the original basis. Observe that
  ! simplification due to A=B are not valid for the
  ! exponent index, eq. P-A=/=0.

  do iZeta=1,nZeta
    PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
    Tmp0 = Zeta(iZeta)+Gmma
    Tmp1 = exp(-Zeta(iZeta)*Gmma*PTC2/Tmp0)
    Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
    Array(ipZ+iZeta-1) = Tmp0
    Array(ipZI+iZeta-1) = One/Tmp0
    Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gmma*TC(1))/Tmp0
    Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gmma*TC(2))/Tmp0
    Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gmma*TC(3))/Tmp0
  end do

  do iDAO=1,nDAO
    ipDAOt = ipDAO+nZeta*(iDAO-1)
    Array(ipDAOt:ipDAOt+nZeta-1) = FactECP*Array(ipK:ipK+nZeta-1)*Array(ipZI:ipZI+nZeta-1)*Two*Pi*DAO(:,iDAO)
  end do

  ! Compute integrals with the Rys quadrature.

  jfg(:) = ifg

  jfgrd(:,:) = ifgrd
  jfhss(:,:,:,:) = ifhss

  jndgrd(:,:,0:nirrep-1) = indgrd(:,:,0:nirrep-1)
  jndhss(:,:,:,:,0:nirrep-1) = indhss(:,:,:,:,0:nirrep-1)

  call Rysg2(iAng,nRys,nZeta,Alpha,Beta,[One],[One],Array(ipZ),Array(ipZI),nZeta,[One],[One],1,Array(ipPx),nZeta,TC,1,Coori,Coor, &
             CoorAC,Array(ip),nArray-ip+1,TNAI1,Fake,Cff2D,Array(ipDAO),nDAO,Hess,nHess,jfGrd,jndGrd,jfHss,jndHss,nOp,iuvwx,jfg, &
             nGr,Indx,lgrad,lhess,tr)
  if (lGrad) then
    nb = nzeta*nTri_Elem1(iAng(1))*nTri_Elem1(iAng(2))
    do iElem=1,nTri_Elem1(iAng(1))*nTri_Elem1(iAng(2))*ngr
      indi = ip+(iElem-1)*nZeta
      Array(indi:indi+nZeta-1) = FactECP*Two*PI*Array(ipK:ipK+nZeta-1)*Array(ipZI:ipZI+nZeta-1)*Array(indi:indi+nZeta-1)
    end do

    call SmAdNa(Array(ip),nb,rFinal,nop(1:3),loper,jndGrd,iuvwx(1:3),Indx,idcar,One,tr)
  end if

end do

return

end subroutine m1kernel
