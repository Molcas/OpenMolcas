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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine projects out specified vectors from the density matrix. *
! The procedure is to read the projection orbitals, generate a         *
! complimentary space, and orthonormalize it. Then the density matrix  *
! is transformed into this representation, the corresponding matrix    *
! are zeroed out, and transformed back to AO representaion.            *
!                                                                      *
!======================================================================*
!                                                                      *
! NOTE: The transformations and orthonormalizations are done as        *
!       (On^4) processes and not as O(n^3) processed due to a lazy     *
!       programmer. I will rewrite this when I get the time.           *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine Proj1()

use Genano_globals, only: MxLqn, nPrim, iSymBk, Ssym, tDsym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: err, i, iB, iBlk, iLqn, ind, indT, iO, iShell, j, jB, jO, Lu, nB, nO
real(kind=wp) :: s
real(kind=wp), allocatable :: pOrb(:,:), DAO(:,:), DMO(:,:), Sc(:,:), Ovlp(:,:)
integer(kind=iwp), external :: IsFreeUnit

Lu = IsFreeUnit(17)
call molcas_open(Lu,'PROJ')
iBlk = 0
!--- Loop over l quantum number ---*
do iLqn=0,MxLqn
  read(Lu,*,iostat=err) nB,nO
  if (err /= 0) exit
  if (nB*nO <= 0) cycle
  if (nB /= nPrim(iLqn)) then
    write(u6,*) 'Project: inconsistency in number of functions'
    call Abend()
  end if
  call mma_allocate(pOrb,nB,nB,label='pOrb')
  call mma_allocate(DAO,nB,nB,label='DAO')
  call mma_allocate(DMO,nB,nB,label='DMO')
  call mma_allocate(Sc,nB,nB,label='Sc')
  call mma_allocate(Ovlp,nB,nB,label='Ovlp')
  do iB=1,nB
    read(Lu,*) (pOrb(iB,iO),iO=1,nO)
  end do
  !--- Generate complementary space ---*
  do iB=1,nB
    do iO=nO+1,nB
      pOrb(iB,iO) = cos(Pi*iO*iB/nB)
    end do
  end do
  !write(u6,'(a,i5)') ' Projection orbitals',iLqn
  !do iB=1,nB
  !  write(u6,'(1x,6f12.6)') (pOrb(iB,iO)),iO=1,nO)
  !end do
  !--- Loop over m quantum number ---*
  do iShell=-iLqn,iLqn
    iBlk = iBlk+1
    !write(u6,'(a,2i5)') ' Density matrix block',iLqn,iShell
    !call Triprt(' ','(6f12.6)',tDsym(iSymBk(iBlk)),nPrim(iLqn))
    !write(u6,'(a,2i5)') ' Overlap matrix block',iLqn,iShell
    !call Triprt(' ','(6f12.6)',Ssym(iSymBk(iBlk)),nPrim(iLqn))
    !--- Copy and square D and S ---*
    do i=1,nB
      do j=1,nB
        indT = min(i,j)+max(i,j)*(max(i,j)-1)/2
        DAO(i,j) = tDsym(iSymBk(iBlk)-1+indT)
        Ovlp(i,j) = Ssym(iSymBk(iBlk)-1+indT)
      end do
    end do
    !--- Orthonormalize orbitals ---*
    do iO=1,nB
      s = Zero
      do iB=1,nB
        do jB=1,nB
          s = s+pOrb(iB,iO)*ovlp(iB,jB)*pOrb(jB,iO)
        end do
      end do
      s = One/sqrt(s)
      do iB=1,nB
        pOrb(iB,iO) = s*pOrb(iB,iO)
      end do
      do jO=iO+1,nB
        s = Zero
        do iB=1,nB
          do jB=1,nB
            s = s+pOrb(iB,iO)*ovlp(iB,jB)*pOrb(jB,jO)
          end do
        end do
        do iB=1,nB
          pOrb(iB,jO) = pOrb(iB,jO)-s*pOrb(iB,iO)
        end do
      end do
    end do
    !--- Transform to MO basis ---*
    do iB=1,nB
      do iO=1,nB
        s = Zero
        do jB=1,nB
          s = s+ovlp(iB,jB)*pOrb(jB,iO)
        end do
        Sc(iB,iO) = s
      end do
    end do
    do iO=1,nB
      do jO=1,nB
        s = Zero
        do iB=1,nB
          do jB=1,nB
            s = s+Sc(iB,iO)*DAO(iB,jB)*Sc(jB,jO)
          end do
        end do
        DMO(iO,jO) = s
      end do
    end do
    !--- Project ---*
    do iO=1,nO
      do jO=1,nB
        DMO(iO,jO) = Zero
        DMO(jO,iO) = Zero
      end do
    end do
    !--- Transform back to AO ---*
    ind=0
    do iB=1,nB
      do jB=1,iB
        ind = ind+1
        s = Zero
        do iO=1,nB
          do jO=1,nB
            s = s+pOrb(iB,iO)*DMO(iO,jO)*pOrb(jB,jO)
          end do
        end do
        tDsym(iSymBk(iBlk)-1+ind) = s
      end do
    end do
    call mma_deallocate(pOrb)
    call mma_deallocate(DAO)
    call mma_deallocate(DMO)
    call mma_deallocate(Sc)
    call mma_deallocate(Ovlp)
    !--- ---*
    !write(u6,'(a,2i5)') ' Density matrix block',iLqn,iShell
    !call Triprt(' ','(6f12.6)',tDsym(iSymBk(iBlk)),nPrim(iLqn))
  end do
end do
!--- Finished ---*
close(Lu)

return

end subroutine Proj1
