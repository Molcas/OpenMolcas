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

subroutine individual_ranks(nDIMCF,BC,BS,Hinit,LJ,iprint)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDIMCF, iprint
real(kind=wp), intent(in) :: BC(nDIMcf,0:nDIMcf), BS(nDIMcf,0:nDIMcf)
complex(kind=wp), intent(in) :: Hinit(nDIMcf,nDIMcf)
character, intent(in) :: LJ
integer(kind=iwp) :: i, ik, ikmax, info, ip, iq, ir, j, jEnd, k, m, n, nfields
real(kind=wp) :: Tnrm, TnrmKQ, wt
complex(kind=wp) :: redME
character(len=16) :: field(8)
character(len=6) :: iprog
integer(kind=iwp), allocatable :: projKQ(:), rankKQ(:)
real(kind=wp), allocatable :: ListKQ(:), Rnrm(:), RnrmKQ(:,:), Snrm(:), Winit(:), Wk(:,:), Ws(:,:), Wtmp(:)
complex(kind=wp), allocatable :: HCF(:,:,:), HCFS(:,:,:), HKQ(:,:), O(:,:), Tmp(:,:), W(:,:), Ztmp(:,:)
real(kind=wp), external :: dznrm2_

!-----------------------------------------------------------------------
call mma_allocate(Rnrm,nDIMCF,label='Rnrm')
call mma_allocate(Snrm,nDIMCF,label='Snrm')
call mma_allocate(Wk,nDIMCF,nDIMCF,label='Wk')
call mma_allocate(Ws,nDIMCF,nDIMCF,label='Ws')
call mma_allocate(Wtmp,nDIMCF,label='Wtmp')
call mma_allocate(Winit,nDIMCF,label='Winit')
call mma_allocate(RnrmKQ,[1,nDIMcf],[-nDIMcf,nDIMcf],label='RnrmKQ')
call mma_allocate(ListKQ,nDIMCF*(2*nDIMCF+1),label='ListKQ')
call mma_allocate(rankKQ,nDIMCF*(2*nDIMCF+1),label='rankKQ')
call mma_allocate(projKQ,nDIMCF*(2*nDIMCF+1),label='projKQ')

call mma_allocate(O,nDIMCF,nDIMCF,label='O')
call mma_allocate(W,nDIMCF,nDIMCF,label='W')
call mma_allocate(HKQ,nDIMCF,nDIMCF,label='HKQ')
call mma_allocate(HCF,nDIMCF,nDIMCF,nDIMCF,label='HCF')
call mma_allocate(HCFS,nDIMCF,nDIMCF,nDIMCF,label='HCFS')
call mma_allocate(Tmp,nDIMCF,nDIMCF,label='Tmp')
call mma_allocate(Ztmp,nDIMCF,nDIMCF,label='Ztmp')

Rnrm(:) = Zero
Snrm(:) = Zero
Tnrm = Zero

HCF(:,:,:) = cZero
HCFS(:,:,:) = cZero
!-----------------------------------------------------------------------
! re-construct the  initial CF matrix:
do N=2,nDIMcf-1,2
  do M=0,N
    call Liviu_ESO(nDIMcf,N,M,O,W,redME)
    if (M == 0) then
      HCF(N,:,:) = HCF(N,:,:)+BC(N,0)*O(:,:)
    else
      HCF(N,:,:) = HCF(N,:,:)+BC(N,M)*O(:,:)+BS(N,M)*W(:,:)
    end if
  end do
end do

do k=2,nDIMCF-1,2
  do ik=2,k
    ! add the ranks as follows:
    ! O2= O2
    ! O4= O2+O4
    ! O6= O2+O4+O6
    ! O8= O2+O4+O6+O8
    ! etc...
    ! compute the cumulative CF matrix
    HCFS(k,:,:) = HCFS(k,:,:)+HCF(ik,:,:)
  end do
end do

do N=2,nDIMCF-1,2
  Tmp(:,:) = HCF(N,:,:)
  Rnrm(N) = dznrm2_(nDIMcf*nDIMcf,Tmp,1)
  Tnrm = Tnrm+Rnrm(N)
  do i=2,N,2
    Snrm(N) = Snrm(N)+Rnrm(i)
  end do
end do
! compute the CF spinting of individual weight operators:
do k=2,nDIMcf-1,2
  Tmp(:,:) = HCF(k,:,:)
  call DIAG_C2(Tmp,nDIMcf,info,Wtmp,Ztmp)
  Wk(k,:) = Wtmp(:)
  Tmp(:,:) = HCFS(k,:,:)
  call DIAG_C2(Tmp,nDIMcf,info,Wtmp,Ztmp)
  Ws(k,:) = Wtmp(:)
end do
! set the initial energies as to the sum of all contributions:
call DIAG_C2(Hinit,nDIMcf,info,Winit,Ztmp)
call mma_deallocate(Wtmp)
call mma_deallocate(Tmp)

!-----------------------------------------------------------------------
! individual parameter contribution:
TnrmKQ = Zero
RnrmKQ(:,:) = Zero
ListKQ(:) = Zero
projKQ(:) = 0
rankKQ(:) = 0
ik = 0
do N=2,nDIMcf-1,2
  do M=-N,N
    ik = ik+1
    ! generate ITO operators:
    call Liviu_ITO(nDIMcf,N,abs(M),O,W,redME)
    ! generate HCF for each parameter rank and projection:
    if (M < 0) then
      HKQ(:,:) = BS(N,abs(M))*W(:,:)
    else if (M >= 0) then
      HKQ(:,:) = BC(N,abs(M))*O(:,:)
    end if
    ! find the rank value of each NM operator
    ! for further estimate the effect of the corresponding CF parameter
    RnrmKQ(N,M) = dznrm2_(nDIMcf*nDIMcf,HKQ,1)
    TnrmKQ = TnrmKQ+RnrmKQ(N,M)

    ! indexing lists
    ListKQ(ik) = RnrmKQ(N,M)
    rankKQ(ik) = N
    projKQ(ik) = M
  end do
end do
ikmax = ik
! make an ordered list of CF parameters
! following their descending effect:
call sort_KQ(ikmax,ListKQ,rankKQ,projKQ,2)
!--------- below we just print the data --------------------------------
if (LJ == 'L') then
  iprog = 'RASSCF'
else
  iprog = ' RASSI'
end if

if (iprint >= 4) then
  write(u6,'(A,F20.12)') 'Tnrm=',Tnrm
  do k=1,nDIMCF-1
    write(u6,'(2(A,i2,A,F20.12,2x))') 'Rnrm(',k,')=',Rnrm(k),'Snrm(',k,')=',Snrm(k)
  end do
end if

! print the cumulative and individual rank weight:
write(u6,'(/)')
write(u6,'(A)') 'CUMULATIVE WEIGHT OF INDIVIDUAL-RANK OPERATORS ON THE CRYSTAL FIELD SPLITTING:'
if ((nDIMCF-1) >= 2) write(u6,'(2x,A,F10.6,A)') 'O2 :------------------------------------------: ',(Snrm(2)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 4) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 :-------------------------------------: ',(Snrm(4)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 6) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 :--------------------------------: ',(Snrm(6)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 8) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 :---------------------------: ',(Snrm(8)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 10) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 :---------------------: ',(Snrm(10)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 12) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 + O12 :---------------: ',(Snrm(12)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 14) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 + O12 + O14 :---------: ',(Snrm(14)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 16) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 + O12 + O14 + O16 :---: ',(Snrm(16)/Tnrm)*100.0_wp,' %.'

!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(A)') 'ENERGY SPLITTING INDUCED BY CUMMULATIVE INDIVIDUAL-RANK OPERATORS (in cm-1).'
field(1) = '      O2       |'
field(2) = '     O2+O4     |'
field(3) = '   O2+O4+O6    |'
field(4) = '  O2+O4+O6+O8  |'
field(5) = ' O2+O4+...+O10 |'
field(6) = ' O2+O4+...+O12 |'
field(7) = ' O2+O4+...+O14 |'
field(8) = ' O2+O4+...+O16 |'
nfields = (nDIMCF-1)/2
do J=1,nfields,4
  jEnd = min(nfields,J+3)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)

  if (mod(nDIMcf,2) == 1) then
    write(u6,'(3A,I2,3x,A,5x,A,4x,A,10(5x,A,6x,A))') '  ',LJ,' =',(nDIMcf-1)/2,'|',iprog,'|',('ONLY','|',i=j,jEnd)
  else
    write(u6,'(3A,I3,A,5x,A,4x,A,10(5x,A,6x,A))') ' ',LJ,' =',nDIMcf-1,'/2 |',iprog,'|',('ONLY','|',i=j,jEnd)
  end if

  write(u6,'(10x,10A)') '|     INITIAL   |',(field(i),i=j,jEnd)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
  do i=1,nDIMcf
    write(u6,'(1x,A,1x,i2,2x,a,12(f14.8,1x,a))') 'w.f.',i,'|',(Winit(i)-Winit(1)),'|',((Ws(k,i)-Ws(k,1)),'|',k=2*j,2*jEnd,2)
  end do
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
end do ! j

write(u6,'(/)')
write(u6,'(A)') 'WEIGHT OF INDIVIDUAL-RANK OPERATORS ON THE CRYSTAL FIELD SPLITTING:'
if ((nDIMCF-1) >= 2) write(u6,'(2x,A,F10.6,A)') 'O2  :-----------------------------------------: ',(Rnrm(2)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 4) write(u6,'(2x,A,F10.6,A)') 'O4  :-----------------------------------------: ',(Rnrm(4)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 6) write(u6,'(2x,A,F10.6,A)') 'O6  :-----------------------------------------: ',(Rnrm(6)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 8) write(u6,'(2x,A,F10.6,A)') 'O8  :-----------------------------------------: ',(Rnrm(8)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 10) write(u6,'(2x,A,F10.6,A)') 'O10 :-----------------------------------------: ',(Rnrm(10)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 12) write(u6,'(2x,A,F10.6,A)') 'O12 :-----------------------------------------: ',(Rnrm(12)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 14) write(u6,'(2x,A,F10.6,A)') 'O14 :-----------------------------------------: ',(Rnrm(14)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 16) write(u6,'(2x,A,F10.6,A)') 'O16 :-----------------------------------------: ',(Rnrm(16)/Tnrm)*100.0_wp,' %.'

!-----------------------------------------------------------------------
write(u6,'(/)')
write(u6,'(A)') 'ENERGY SPLITTING INDUCED BY INDIVIDUAL-RANK OPERATORS (in cm-1).'
field(1) = '      O2       |'
field(2) = '      O4       |'
field(3) = '      O6       |'
field(4) = '      O8       |'
field(5) = '      O10      |'
field(6) = '      O12      |'
field(7) = '      O14      |'
field(8) = '      O16      |'
nfields = (nDIMCF-1)/2

do J=1,nfields,4
  jEnd = min(nfields,J+3)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)

  if (mod(nDIMcf,2) == 1) then
    write(u6,'(3A, I2,3x,A,5x,A,4x,A,10(5x,A,6x,A))') '  ',LJ,' =',(nDIMcf-1)/2,'|',iprog,'|',('ONLY','|',i=j,jEnd)
  else
    write(u6,'(3A,I3,a,5x,A,4x,A,10(5x,A,6x,A))') ' ',LJ,' =',nDIMcf-1,'/2 |',iprog,'|',('ONLY','|',i=j,jEnd)
  end if

  write(u6,'(10x,10A)') '|     INITIAL   |',(field(i),i=j,jEnd)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
  do i=1,nDIMcf
    write(u6,'(1x,A,1x,i2,2x,a,12(f14.8,1x,a))') 'w.f.',i,'|',(Winit(i)-Winit(1)),'|',((Wk(k,i)-Wk(k,1)),'|',k=2*j,2*jEnd,2)
  end do
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
end do ! j

!-----------------------------------------------------------------------
write(u6,'(/)')
write(u6,'(A)') 'WEIGHT OF INDIVIDUAL CRYSTAL FIELD PARAMETERS ON THE CRYSTAL FIELD SPLITTING: (in descending order):'
write(u6,'(A)') 'CFP are given in ITO used in J. Chem. Phys. 137, 064112 (2012).'

write(u6,'(2A)') repeat('-',55),'|'
write(u6,'(A)') '  k |  q  |         B(k,q)        |    Weight (in %)   |'
write(u6,'(A)') '----|-----|-----------------------|--------------------|'
do ik=1,ikmax
  ip = projKQ(ik)
  iq = abs(projKQ(ik))
  ir = rankKQ(ik)
  wt = 100.0_wp*ListKQ(ik)/TnrmKQ

  if (projKQ(ik) >= 0) then
    write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),(ES22.14,1x,A),F19.14,1x,A)') ir,'|',ip,'|',BC(ir,iq),'|',wt,'|'
  else
    write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),(ES22.14,1x,A),F19.14,1x,A)') ir,'|',ip,'|',BS(ir,iq),'|',wt,'|'
  end if
end do
write(u6,'(A)') '----|-----|-----------------------|--------------------|'

!-----------------------------------------------------------------------
call mma_deallocate(Rnrm)
call mma_deallocate(Snrm)
call mma_deallocate(Wk)
call mma_deallocate(Ws)
call mma_deallocate(Winit)

call mma_deallocate(O)
call mma_deallocate(W)
call mma_deallocate(HKQ)
call mma_deallocate(HCF)
call mma_deallocate(HCFS)
call mma_deallocate(Ztmp)
call mma_deallocate(RnrmKQ)
call mma_deallocate(ListKQ)
call mma_deallocate(rankKQ)
call mma_deallocate(projKQ)

return

end subroutine individual_ranks
