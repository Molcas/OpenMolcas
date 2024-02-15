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

subroutine CRYSTALFIELD_1(nDIMcf,nlanth,MM,ESOJ,GRAD,iprint)
! This soubrutine calculates the crystal field parameters on the basis
! of the given fron RASSI - J multiplet.
! In a second step, the first largest 27 parameters will be used to
! recalculate the S-O energies, eigenfunctions, g- and D- tensors.

!  Employed parameters:

!  IPRINT = the print level of the calculation

!  IReturn = the error value.
!        0 = no error, happy landing

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDIMcf, nlanth, iprint
complex(kind=wp), intent(in) :: MM(3,nDIMcf,nDIMcf)
real(kind=wp), intent(in) :: ESOJ(nDIMcf)
logical(kind=iwp), intent(in) :: GRAD
integer(kind=iwp) :: i, info, j, k, LuCF, q
real(kind=wp), allocatable :: BNC(:,:), BNS(:,:), Bstev(:,:), Eloc(:), Winit(:)
complex(kind=wp), allocatable :: Akq(:,:), HCF(:,:), Z(:,:), Zinit(:,:)
integer(kind=iwp), external :: IsFreeUnit

call mma_allocate(Winit,nDIMcf,'Winit')
call mma_allocate(Eloc,nDIMcf,'Eloc')
call mma_allocate(Zinit,nDIMcf,nDIMcf,'Zinit')
call mma_allocate(Z,nDIMcf,nDIMcf,'Z')
call mma_allocate(HCF,nDIMcf,nDIMcf,'HCF')
call mma_allocate(BNC,[1,nDIMcf],[0,nDIMcf],label='BNC')
call mma_allocate(BNS,[1,nDIMcf],[0,nDIMcf],label='BNS')
call mma_allocate(Bstev,[1,nDIMcf],[-nDIMcf,nDIMcf],label='Bstev')
call mma_allocate(Akq,[1,nDIMcf-1],[-(nDIMcf-1),nDIMcf-1],label='Akq')

! find the J-pseudospin:
!iDir = 3
call pseudospin(MM,nDIMcf,Z,3,1,iprint)
call rtrace(nDIMcf,ESOJ,ELOC)
! re-write the CF matrix in J-pseudospin basis:
! energy units = cm-1
HCF(:,:) = cZero
do j=1,nDIMcf
  do k=1,nDIMcf
    HCF(:,j) = HCF(:,j)+ELOC(k)*conjg(Z(k,:))*Z(k,j)
  end do
end do
! diagonalize the initial CF matrix:
call DIAG_C2(HCF,nDIMcf,INFO,Winit,Zinit)
call print_ZFS('Ab Initio Calculated Crystal-Field Splitting Matrix written in the basis of Pseudospin Eigenfunctions',HCF,nDIMCF)

if (IPRINT > 2) then
  write(u6,*)
  write(u6,'(5X,A)') 'MAIN VALUES OF THE INITIAL CRYSTAL-FIELD HAMILTONIAN:'
  write(u6,*)
  if (mod(nDIMcf,2) == 1) then
    do I=1,nDIMcf
      write(u6,'(3X,A,I3,A,F25.16)') '|',(nDIMcf-1)/2+(1-I),'> = ',Winit(I)-Winit(1)
    end do
  else
    do I=1,nDIMcf
      write(u6,'(3X,A,I3,A,F25.16)') '|',(nDIMcf-1)-2*(I-1),'/2 > = ',Winit(i)-Winit(1)
    end do
  end if
  write(u6,*)

  write(u6,'(5X,A)') 'EIGENVECTORS OF THE INITIAL CRYSTAL-FIELD HAMILTONIAN:'
  write(u6,*)
  call print_ZFS_naoya('J',Zinit,nDIMcf)
  ! End the checking of the main values of the initial crystal-field
end if

! calculating the coeficients of the crystal filed operators Bnm
!  Akq=(2k+1)/(2J+1) * 1/|< J || O || J >|^2 * Trace{HCF*O(k,-q)}
call NEWCF(HCF,nDIMcf,Akq,BNC,BNS,Bstev)
!#ifdef _DEBUGPRINT_
!call recover_CF(nDIMCF,HCF,Akq,BNC,BNS,Bstev)
!#endif

call print_CFP_alpha(nlanth,nDIMCF,BNC,BNS)
if (iprint >= 4) then
  call print_CFP_LCLU(nDIMCF,BNC,BNS,.true.)
  call print_CFP_stev(nDIMCF,Bstev,.true.)
  call print_CFP_naoya(nDIMcf,Akq,.true.)
else
  call print_CFP_LCLU(nDIMCF,BNC,BNS,.false.)
  call print_CFP_stev(nDIMCF,Bstev,.false.)
  call print_CFP_naoya(nDIMcf,Akq,.false.)
end if

!-----------------------------------------------------------------------
write(u6,'(/)')
if (mod(nDIMcf,2) == 1) then
  write(u6,'(A,I0)') 'DECOMPOSITION OF THE RASSI WAVE FUNCTIONS CORRESPONDING TO THE LOWEST ATOMIC MULTIPLET J =',(nDIMcf-1)/2
  write(u6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE PROJECTION OF THE TOTAL MOMENT ON THE QUANTIZATION AXIS'
else ! mod(nDIMcf,2) == 0
  write(u6,'(A,I0,A)') 'DECOMPOSITION OF THE RASSI WAVE FUNCTIONS CORRESPONDING TO THE LOWEST ATOMIC MULTIPLET J = ',(nDIMcf-1),'/2'
  write(u6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE PROJECTION OF THE TOTAL MOMENT ON THE QUANTIZATION AXIS'
end if

call print_ZFS_naoya('J',Zinit,nDIMcf)
call individual_ranks(nDIMCF,BNC,BNS,HCF,'J',iprint)
!-----------------------------------------------------------------------
! saving some information for tests:
call Add_Info('CRYS_BNMC_20',BNC(2,0),1,4)
call Add_Info('CRYS_BNMC_40',BNC(4,0),1,4)
call Add_Info('CRYS_BNMC_60',BNC(6,0),1,4)
!-----------------------------------------------------------------------
! for the interface related to CF gradient calculation:
if (GRAD) then
  LuCF = IsFreeUnit(81)
  call molcas_open(LuCF,'CFMAT')
  do k=2,nDIMcf-1,2
    do q=0,k
      write(LuCF,'(I3,I3,1x,2ES25.15)') k,q,BNC(k,q),BNS(k,q)
    end do
  end do
  close(LuCF)
end if
!-----------------------------------------------------------------------
call mma_deallocate(Winit)
call mma_deallocate(Eloc)
call mma_deallocate(Zinit)
call mma_deallocate(Z)
call mma_deallocate(HCF)
call mma_deallocate(BNC)
call mma_deallocate(BNS)
call mma_deallocate(Bstev)
call mma_deallocate(Akq)

return

end subroutine CRYSTALFIELD_1
