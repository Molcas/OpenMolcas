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

subroutine Wr_Prop(nAtoms,nCenters,nBas,nMltPl,NOCOB,NOCOB_b,orbe,orbe_b,iPol,LAllCenters)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nCenters, nBas, nMltPl, NOCOB, NOCOB_b, iPol
real(kind=wp), intent(in) :: orbe(NOCOB), orbe_b(NOCOB_b)
logical(kind=iwp), intent(in) :: LAllCenters
#include "MpData.fh"
#include "WrkSpc.fh"
#include "MolProp.fh"
integer(kind=iwp) :: i, iComp, iCompMat(0:nMltPl,0:nMltPl,0:nMltPl), il, iMltpl, ip, iq, j, nA, nB, nComp, nl, np, nq, nTotCen
real(kind=wp) :: fac, rnloveril, rnPoveriP, rnqoveriq, xfac, yfac, zfac
character(len=8) :: MemLabel

nTotCen = 0
do i=1,nAtoms
  nTotCen = nTotCen+1
  write(CEN_LAB(i*(i+1)/2),'(A)') Labe(i)
  do j=1,i
    if (BondMat(i,j)) then
      nTotCen = nTotCen+1
      write(CEN_LAB(i*(i-1)/2+j),'(3A)') LABE(i),'- ',LABE(j)
    end if
  end do
end do

do iMltpl=0,nMltPl
  iComp = 0
  nComp = (iMltPl+1)*(iMltPl+2)/2
  write(MemLabel,'(A4,I4.4)') 'AtTo',iMltPl
  call GetMem(MemLabel,'Allo','Real',iAtMltPlTotAd(iMltPl),nComp)
  do i=0,nComp-1
    Work(iAtMltPlTotAd(iMltPl)+i) = Zero
  end do
  write(MemLabel,'(A4,I4.4)') 'BoTo',iMltPl
  call GetMem(MemLabel,'Allo','Real',iAtBoMltPlTotAd(iMltPl),nComp)
  do i=0,nComp-1
    Work(iAtBoMltPlTotAd(iMltPl)+i) = Zero
  end do
  do np=iMltpl,0,-1
    do nq=iMltpl-np,0,-1
      nl = iMltpl-np-nq
      iComp = iComp+1
      iCompMat(np,nq,nl) = iComp
      do nA=1,nAtoms
        do ip=0,np
          call NoverP(np,ip,rnPoveriP)
          if (np == ip) then
            xfac = rnpoverip
          else
            xfac = rnPoveriP*(Cor(1,nA,nA))**(np-ip)
          end if
          do iq=0,nq
            call NoverP(nq,iq,rnqoveriq)
            if (nq == iq) then
              yfac = rnqoveriq
            else
              yfac = rnqoveriq*(Cor(2,nA,nA))**(nq-iq)
            end if
            do il=0,nl
              call NoverP(nl,il,rnloveril)
              if (nl == il) then
                zfac = rnloveril
              else
                zfac = rnloveril*(Cor(3,nA,nA))**(nl-il)
              end if
              fac = xfac*yfac*zfac*Work(iAtMltPlAd(ip+iq+il)+nAtoms*(iCompMat(ip,iq,il)-1)+nA-1)
              Work(iAtMltPlTotAd(iMltpl)+iComp-1) = Work(iAtMltPlTotAd(iMltpl)+iComp-1)+fac
            end do
          end do
        end do
        do nB=1,nA
          if ((nA == nB) .or. BondMat(nA,nB)) then
            do ip=0,np
              call NoverP(np,ip,rnPoveriP)
              if (np == ip) then
                xfac = rnpoverip
              else
                xfac = rnPoveriP*(Cor(1,nA,nB))**(np-ip)
              end if
              do iq=0,nq
                call NoverP(nq,iq,rnqoveriq)
                if (nq == iq) then
                  yfac = rnqoveriq
                else
                  yfac = rnqoveriq*(Cor(2,nA,nB))**(nq-iq)
                end if
                do il=0,nl
                  call NoverP(nl,il,rnloveril)
                  if (nl == il) then
                    zfac = rnloveril
                  else
                    zfac = rnloveril*(Cor(3,nA,nB))**(nl-il)
                  end if
                  fac = xfac*yfac*zfac*Work(iAtBoMltPlAd(ip+iq+il)+nCenters*(iCompMat(ip,iq,il)-1)+nA*(nA-1)/2+nB-1)
                  Work(iAtBoMltPlTotAd(iMltpl)+iComp-1) = Work(iAtBoMltPlTotAd(iMltpl)+iComp-1)+fac
                end do
              end do
            end do
          end if
        end do
      end do
    end do
  end do
end do

call Wr_MpProp(nAtoms,nCenters,nMltPl,iPol)
!EB call Wr_Files(nAtoms,nCenters,nMltPl,nBas,NOCOB,orbe,iBond,
call Wr_Files(nAtoms,nCenters,nMltPl,nBas,NOCOB,NOCOB_b,orbe,orbe_b,LAllCenters)

do iMltpl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  write(MemLabel,'(A4,I4.4)') 'AtTo',iMltPl
  call GetMem(MemLabel,'Free','Real',iAtMltPlTotAd(iMltPl),iMltPl*nComp)
  write(MemLabel,'(A4,I4.4)') 'BoTo',iMltPl
  call GetMem(MemLabel,'Free','Real',iAtBoMltPlTotAd(iMltPl),iMltPl*nComp)
end do

return

!EB 111 format(A,3F15.5,6F10.3,10I5)

end subroutine Wr_Prop
