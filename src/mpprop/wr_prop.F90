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

use MPProp_globals, only: AtBoMltPl, AtBoMltPlTot, AtMltPl, AtMltPlTot, BondMat, Cen_Lab, Cor, Labe
use Data_Structures, only: Allocate_DT, Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nCenters, nBas, nMltPl, NOCOB, NOCOB_b, iPol
real(kind=wp), intent(in) :: orbe(NOCOB), orbe_b(NOCOB_b)
logical(kind=iwp), intent(in) :: LAllCenters
integer(kind=iwp) :: i, iComp, il, iMltpl, ip, iq, j, nA, nB, nComp, nl, np, nq, nTotCen
real(kind=wp) :: fac, rnloveril, rnPoveriP, rnqoveriq, xfac, yfac, zfac
character(len=8) :: MemLabel
integer(kind=iwp), allocatable :: iCompMat(:,:,:)

call mma_allocate(Cen_Lab,nAtoms*(nAtoms+1)/2,label='Cen_Lab')
call mma_allocate(iCompMat,[0,nMltPl],[0,nMltPl],[0,nMltPl],label='iCompMat')

nTotCen = 0
do i=1,nAtoms
  nTotCen = nTotCen+1
  write(Cen_Lab(i*(i+1)/2),'(A)') Labe(i)
  do j=1,i-1
    if (BondMat(i,j)) then
      nTotCen = nTotCen+1
      write(Cen_Lab(i*(i-1)/2+j),'(3A)') LABE(i),'- ',LABE(j)
    end if
  end do
end do

call Allocate_DT(AtMltPlTot,[0,nMltPl],'AtMltPlTot')
call Allocate_DT(AtBoMltPlTot,[0,nMltPl],'AtBoMltPlTot')
do iMltpl=0,nMltPl
  iComp = 0
  nComp = (iMltPl+1)*(iMltPl+2)/2
  write(MemLabel,'(a4,i4.4)') 'AtTo',iMltPl
  call mma_allocate(AtMltPlTot(iMltPl)%A,nComp,label=MemLabel)
  AtMltPlTot(iMltPl)%A(:) = Zero
  write(MemLabel,'(a4,i4.4)') 'BoTo',iMltPl
  call mma_allocate(AtBoMltPlTot(iMltPl)%A,nComp,label=MemLabel)
  AtBoMltPlTot(iMltPl)%A(:) = Zero
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
              fac = xfac*yfac*zfac*AtMltPl(ip+iq+il)%A(iCompMat(ip,iq,il),nA)
              AtMltPlTot(iMltpl)%A(iComp) = AtMltPlTot(iMltpl)%A(iComp)+fac
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
                  fac = xfac*yfac*zfac*AtBoMltPl(ip+iq+il)%A(iCompMat(ip,iq,il),nA*(nA-1)/2+nB)
                  AtBoMltPlTot(iMltpl)%A(iComp) = AtBoMltPlTot(iMltpl)%A(iComp)+fac
                end do
              end do
            end do
          end if
        end do
      end do
    end do
  end do
end do

call mma_deallocate(iCompMat)

call Wr_MpProp(nAtoms,nCenters,nMltPl,iPol)
!EB call Wr_Files(nAtoms,nCenters,nMltPl,nBas,NOCOB,orbe,iBond,
call Wr_Files(nAtoms,nCenters,nMltPl,nBas,NOCOB,NOCOB_b,orbe,orbe_b,LAllCenters)

call Deallocate_DT(AtMltPlTot)
call Deallocate_DT(AtBoMltPlTot)

call mma_deallocate(Cen_Lab)

return

!EB 111 format(A,3F15.5,6F10.3,10I5)

end subroutine Wr_Prop
