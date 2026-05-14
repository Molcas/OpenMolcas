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

subroutine SpinDens(LS,RS,iL,iR,rP1,rp2,rp3,rp4,rp5,rDe1,rde2,itype)
! Input:
!
! LS : CI Coeff for left state
! RS : CI Coeff for right state
! iL : Symmetry of left state
! iR : Symmetry of right state
!
! Output:
!
!  rP1 : Two Electron -- Density
!  rP2 : Two Electron ++ Density
!  rP3 : Two Electron spin Density
!  rD1 : One Electron Spin adapted Density
!  rD1 : One Electron spin density Density
!
! itype=1
! Densities for:
!          0
! <L|[{Q:S} ,H]|R>
!          0
! itype 2
! Densities for
!           0       0
! <L||[{Q:S} ,[{Q:S} ,H]]|R>
!           0       0

use Index_Functions, only: iTri
use MCLR_Data, only: n1Dens, n2Dens, nNA, XISPSM
use CandS, only: ICSM, ISSM
use input_mclr, only: nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: LS(*), RS(*)
integer(kind=iwp), intent(in) :: iL, iR, iType
real(kind=wp), intent(out) :: rP1(nna,nna,nna,nna), rP2(nna,nna,nna,nna), rP3(nna,nna,nna,nna), rP4(nna,nna,nna,nna), &
                              rP5(nna,nna,nna,nna), rDe1(nna,nna), rde2(nna,nna)
integer(kind=iwp) :: iA, ijkl, ijklAB, ijklBA, jA, jilk, jilkAB, jilkBA, kA, lA, n2, nConfL, nConfR
real(kind=wp), allocatable :: CIL(:), CIR(:), Dens(:,:), Pens(:)

n2 = 2*n2Dens+nnA**4
call mma_allocate(Dens,n1Dens,2,Label='Dens')
call mma_allocate(Pens,n2,Label='Pens')
Dens(:,:) = Zero
Pens(:) = Zero
nConfL = max(ncsf(il),nint(xispsm(il,1)))
nConfR = max(ncsf(iR),nint(xispsm(iR,1)))

call mma_allocate(CIL,nConfL,Label='CIL')
call mma_allocate(CIR,nConfR,Label='CIR')
call CSF2SD(LS,CIL,iL)
call CSF2SD(RS,CIR,iR)
icsm = iR
issm = iL
call Densi2_mclr(2,Dens,Pens,CIL,CIR,0,0,1,n1Dens,n2Dens)

if (itype == 1) then

  !      0
  ! <0|[Q  ,H]|0>
  !      0pq

  do iA=1,nna
    do jA=1,nnA
      rde1(ia,ja) = Dens(nna*(ia-1)+ja,1)-Dens(nna*(ia-1)+ja,2)+Dens(nna*(ja-1)+ia,1)-Dens(nna*(ja-1)+ia,2)
    end do
  end do
  do iA=1,NnA
    do jA=1,nnA
      do kA=1,nnA
        do lA=1,nnA
          ijklab = (ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
          jilkab = (ja-1)*nna**3+(ia-1)*nna**2+(la-1)*nna+ka
          ijklba = (ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
          jilkba = (la-1)*nna**3+(ka-1)*nna**2+(ja-1)*nna+ia
          ijkl = iTri((ia-1)*nna+ja,(ka-1)*nna+la)
          jilk = iTri((ja-1)*nna+ia,(la-1)*nna+ka)
          rP1(ia,ja,ka,la) = Pens(ijkl)-Pens(ijkl+n2Dens)+Pens(n2Dens*2+ijklab)-Pens(ijklba+2*n2Dens)+Pens(jilk)- &
                             Pens(jilk+n2Dens)-Pens(n2Dens*2+jilkab)+Pens(jilkba+2*n2Dens)
        end do
      end do
    end do
  end do

else if (itype == 2) then
  ! OK CONSTRUCT

  ! --
  ! ++
  ! -+
  ! spindensity
  ! spin adapted density

  rde1(:,:) = reshape(Dens(:,2)-Dens(:,1),[nna,nna])
  rde2(:,:) = reshape(Dens(:,1)+Dens(:,2),[nna,nna])

  do iA=1,nnA
    do jA=1,nnA
      do kA=1,nnA
        do lA=1,nnA
          ijklab = (ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
          ijklba = (ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
          ijkl = iTri((ia-1)*nna+ja,(ka-1)*nna+la)
          rP1(ia,ja,ka,la) = Pens(ijkl)+Pens(n2Dens+ijkl)-Pens(ijklab+2*n2Dens)-Pens(ijklba+2*n2Dens)
        end do
      end do
    end do
  end do
  do iA=1,nnA
    do jA=1,nnA
      do kA=1,nnA
        do lA=1,nnA
          ijklab = (ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
          ijklba = (ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
          ijkl = iTri((ia-1)*nna+ja,(ka-1)*nna+la)
          rP2(ia,ja,ka,la) = Pens(ijkl)-Pens(n2Dens+ijkl)-Pens(ijklab+2*n2Dens)+Pens(ijklba+2*n2Dens)
        end do
      end do
    end do
  end do

  do iA=1,NnA
    do jA=1,nnA
      do kA=1,nnA
        do lA=1,nnA
          ijklab = (ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
          ijklba = (ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
          ijkl = iTri((ia-1)*nna+ja,(ka-1)*nna+la)
          rP3(ia,ja,ka,la) = Pens(ijkl)+Pens(ijkl+n2Dens)+Pens(n2Dens*2+ijklab)+Pens(ijklba+2*n2Dens)
        end do
      end do
    end do
  end do
  do iA=1,NnA
    do jA=1,nnA
      do kA=1,nnA
        do lA=1,nnA
          ijklab = (ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
          ijklba = (ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
          ijkl = iTri((ia-1)*nna+ja,(ka-1)*nna+la)
          rP4(ia,ja,ka,la) = -Pens(n2Dens*2+ijklab)+Pens(ijklba+2*n2Dens)
        end do
      end do
    end do
  end do
  do iA=1,NnA
    do jA=1,nnA
      do kA=1,nnA
        do lA=1,nnA
          ijklab = (ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
          ijklba = (ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
          ijkl = iTri((ia-1)*nna+ja,(ka-1)*nna+la)
          rP5(ia,ja,ka,la) = Pens(n2Dens*2+ijklab)+Pens(ijklba+2*n2Dens)
        end do
      end do
    end do
  end do
end if
call mma_deallocate(Dens)
call mma_deallocate(Pens)
call mma_deallocate(CIL)
call mma_deallocate(CIR)

end subroutine SpinDens
