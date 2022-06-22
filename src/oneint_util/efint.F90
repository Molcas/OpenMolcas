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
! Copyright (C) 1991,1995, Roland Lindh                                *
!***********************************************************************

subroutine EFInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of electric field         *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!                                                                      *
! Modified for explicit code, R. Lindh, February '95.                  *
!***********************************************************************

implicit real*8(A-H,O-Z)
external TNAI, Fake, XCff2D, XRys2D
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
integer iDCRT(0:7), iStabO(0:7)
real*8 TC(3), Coori(3,4), CoorAC(3,2)
logical EQ, NoSpecial
integer iAnga(4)
character*80 Label
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 200
iPrint = nPrint(iRout)

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1
mcdMin = nabSz(nOrdOp-1)+1
mcdMax = nabSz(nOrdop)
lab = (mabMax-mabMin+1)
kab = nElem(la)*nElem(lb)
lcd = (mcdMax-mcdMin+1)
labcd = lab*lcd

! Compute Flop's and size of work array which HRR will Use.

call mHRR(la,lb,nFLOP,nMem)

! Distribute the work array

ip2 = 1
ip1 = ip2+nZeta*max(labcd,lcd*nMem)
mArr = nArr-max(labcd,lcd*nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CCoor,TC)
  call dcopy_(3,TC,1,CoorAC(1,2),1)
  call dcopy_(3,TC,1,Coori(1,3),1)
  call dcopy_(3,TC,1,Coori(1,4),1)

  ! Compute integrals with the Rys-Gauss quadrature.

  nT = nZeta
  NoSpecial = .true.
  call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coori,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
           Array(ip1),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

  ! The integrals are now ordered as ijkl,e,f

  ! a) Change the order to f,ijkl,e
  ! b) Unfold e to ab, f,ijkl,ab
  ! c) Change the order back to ijkl,ab,f

  ! a)

  call DGetMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array(ip2),lcd)

  ! b) Use the HRR to unfold e to ab

  call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
  ip3 = ip2-1+ipIn

  ! c)

  call DGetMO(Array(ip3),lcd,lcd,nZeta*kab,Array(ip1),nZeta*kab)

  ! Modify to traceless form, the sixth element contains r*r and

  if (nOrdOp == 2) then
    !if (.false.) then
    nzab = nZeta*kab
    iOffxx = ip1
    iOffyy = ip1+nzab*3
    iOffzz = ip1+nzab*5
    ThreeI = One/Three
    do i=0,nzab-1
      RR = Array(iOffxx+i)+Array(iOffyy+i)+Array(iOffzz+i)
      XX = Two*Array(iOffxx+i)-Array(iOffyy+i)-Array(iOffzz+i)
      YY = Two*Array(iOffyy+i)-Array(iOffxx+i)-Array(iOffzz+i)
      Array(iOffxx+i) = XX*ThreeI
      Array(iOffyy+i) = YY*ThreeI
      Array(iOffzz+i) = RR
    end do
  end if

  ! Stored as nZeta,iElem,jElem,iComp

  if (iPrint >= 49) then
    write(6,*) ' In EFInt la,lb=',la,lb
    nzab = nZeta*kab
    do iElem=1,nElem(la)
      do jElem=1,nElem(lb)
        ij = (jElem-1)*nElem(la)+iElem
        ip = ip1+nZeta*(ij-1)
        do iComp=1,nComp
          write(Label,'(A,I2,A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',',iComp,') '
          call RecPrt(Label,' ',Array(ip),nZeta,1)
          ip = ip+nzab
        end do
      end do
    end do
  end if

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ip1),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine EFInt
