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

subroutine barrier(nBlock,dipIn,W,imanifold,NMULT,NDIM,doPLOT,iprint)
! the present Subroutine computes the matrix elements of the transitions from states forming the blocking barrier;
! the states of opposite magnetization are given in the input, under keyword MLTP:
! the magnetic field is applied to each group of states delared at MLTP in order to form states of a definite
! projection of M on the quantization axis.
!  N --  dimension of the barrier

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: nBlock, nMult, iprint, imanifold
integer, intent(in) :: nDim(nMult)
real(kind=8), intent(in) :: W(nBlock)
complex(kind=8), intent(in) :: dipIn(3,nBlock,nBlock)
logical, intent(in) :: doPLOT
!-----------------------------------------------------------------------
integer :: k, l, i, j, i1, i2, il, idim, Ifunct, j1, j2, iMult, ipar
integer :: nb, maxmult
real(kind=8) :: mave
character(len=90) :: string2
character(len=5) :: s1, s2
integer, allocatable :: ibas(:,:)                 ! (nmult,nBlock)
real(kind=8), allocatable :: gtens(:)             ! (3)
real(kind=8), allocatable :: maxes(:,:), E(:)     ! (3,3), E(nmult)
real(kind=8), allocatable :: wz(:,:)              ! (nmult,nBlock)
complex(kind=8), allocatable :: dipN(:,:,:)       ! (3,nBlock,nBlock)
complex(kind=8), allocatable :: dipN3(:,:,:)      ! (3,ndim(imanifold),ndim(imanifold))
complex(kind=8), allocatable :: CZ(:,:,:)         ! (nmult,nBlock,nBlock)
complex(kind=8), allocatable :: Ztr(:,:)          ! (nBlock,nBlock)
complex(kind=8), allocatable :: ML(:,:,:)         ! (3,nBlock,nBlock)
complex(kind=8), allocatable :: tmp(:,:)          ! (nBlock,nBlock)
complex(kind=8), allocatable :: MM(:)             ! (3)
complex(kind=8), allocatable :: dipso5(:,:,:,:,:) ! (3,nmult,10,nmult,10)
!-----------------------------------------------------------------------

if ((nmult > 0) .and. (nBlock > 0)) then
  call mma_allocate(ibas,nmult,nBlock,'ibas')
  call mma_allocate(wz,nmult,nBlock,'wz')
  call mma_allocate(cz,nmult,nBlock,nBlock,'cz')
  call icopy(nmult*nBlock,[0],0,ibas,1)
  call dcopy_(nmult*nBlock,[0.0_wp],0,wz,1)
  call zcopy_(nmult*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,CZ,1)
end if

if (nmult > 0) then
  call mma_allocate(E,nmult,'E')
  call mma_allocate(dipso5,3,nmult,10,nmult,10,'dipso5')
  call dcopy_(nmult,[0.0_wp],0,E,1)
  call zcopy_(3*nmult*10*nmult*10,[(0.0_wp,0.0_wp)],0,dipso5,1)
end if

if (nBlock > 0) then
  call mma_allocate(dipN,3,nBlock,nBlock,'dipN')
  call mma_allocate(ML,3,nBlock,nBlock,'ML')
  call mma_allocate(Ztr,nBlock,nBlock,'Ztr')
  call mma_allocate(tmp,nBlock,nBlock,'tmp')
  call zcopy_(3*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,dipN,1)
  call zcopy_(3*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,ML,1)
  call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,Ztr,1)
  call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,tmp,1)
end if

k = nDim(imanifold)
if (k > 0) then
  call mma_allocate(dipN3,3,k,k,'dipN3')
  call zcopy_(3*k*k,[(0.0_wp,0.0_wp)],0,dipN3,1)
end if

call mma_allocate(gtens,3,'gtens')
call mma_allocate(MM,3,'MM')
call mma_allocate(maxes,3,3,'maxes')
call dcopy_(3,[0.0_wp],0,gtens,1)
call dcopy_(3*3,[0.0_wp],0,maxes,1)
call zcopy_(3,[(0.0_wp,0.0_wp)],0,MM,1)

!-----------------------------------------------------------------------
do i=1,3
  maxes(i,i) = 1.0_wp
end do
! rotate the magnetic moment to the magnetic axes of the ground multiplet ( NDIM(1) )
if (ndim(imanifold) <= 1) then
  write(6,'(a,i2,a)') 'The manifold ',imanifold,' was chosen for determination of the quantization axis.'
  write(6,'(a     )') 'However, the size of this manifold is:'
  write(6,'(a,i1  )') 'size = ',ndim(imanifold)
  write(6,'(a     )') 'in this case, quantization axis will remain the original Z axis'
else
  do i=1,ndim(imanifold)
    do j=1,ndim(imanifold)
      do l=1,3
        dipN3(l,i,j) = dipIn(l,i,j)
      end do
    end do
  end do
  call atens(dipN3,ndim(imanifold),gtens,maxes,1)
end if
call rotmom2(dipIn(1:3,1:nBlock,1:nBlock),nBlock,maxes(1:3,1:3),dipN(1:3,1:nBlock,1:nBlock))
if (iprint > 2) then
  write(6,*)
  write(6,'(10X,A)') 'Magnetic moment (dipIN) in the original coordinate system'
  write(6,*)
  write(6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >','projection = X','projection = Y','projection = Z', &
                                         'ABS(< i | moment | j >)/3'
  write(6,*)
  do I=1,nBlock
    do J=1,nBlock
      MAVE = 0.0_wp
      MAVE = (abs(dipIN(1,I,J))+abs(dipIN(2,I,J))+abs(dipIN(3,I,J)))/3.0_wp
      write(6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)') '<',i,'| moment |',j,' >',(dipIN(L,I,J),L=1,3),MAVE
    end do
  end do
  write(6,'(///)')
  write(6,'(10X,A)') 'Magnetic moment (dipN) in the coordinate system of the magnetic axes'
  write(6,*)
  write(6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >','projection = X','projection = Y','projection = Z', &
                                         'ABS(< i | moment | j >)/3'
  write(6,*)
  do I=1,nBlock
    do J=1,nBlock
      MAVE = 0.0_wp
      MAVE = (abs(dipN(1,I,J))+abs(dipN(2,I,J))+abs(dipN(3,I,J)))/3.0_wp
      write(6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)') '<',i,'| moment |',j,' >',(dipN(L,I,J),L=1,3),MAVE
    end do
  end do
end if
!-----------------------------------------------------------------------
! determine the "CZ matrix":
Ifunct = 1
do il=1,nmult
  idim = ndim(il)
  if (idim == 0) return
  if (idim == 1) then
    CZ(il,1,1) = (1.0_wp,0.0_wp)
    WZ(il,1) = W(Ifunct)
  else ! idim > 1

    call pseudospin(dipN(1:3,Ifunct:(Ifunct+ndim(il)-1),Ifunct:(Ifunct+ndim(il)-1)),ndim(il),CZ(il,1:ndim(il),1:ndim(il)),3,1, &
                    iprint)

    if (iPrint > 2) call pa_prMat('barrier:  CZ',CZ(il,1:nDim(il),1:nDim(il)),nDim(il))
  end if
  Ifunct = Ifunct+ndim(il)
end do
!-----------------------------------------------------------------------
! compute the matrix elements between multiplets
l = 0
ibas = 0
do i=1,nmult
  do j=1,ndim(i)
    l = l+1
    ibas(i,j) = l
  end do
end do

maxmult = 0
do il=1,nmult
  if (ndim(il) > maxmult) maxmult = ndim(il)
end do
! build the transformation matrix Z(nBlock,nBlock)
call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,ZTR,1)
do iMult=1,nmult
  do j1=1,ndim(iMult)
    do j2=1,ndim(iMult)
      i = ibas(iMult,j1)
      j = ibas(iMult,j2)
      Ztr(i,j) = CZ(iMult,j1,j2)
    end do
  end do
end do

call zcopy_(3*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,ML,1)
do L=1,3
  call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,TMP,1)
  call ZGEMM_('C','N',nBlock,nBlock,nBlock,(1.0_wp,0.0_wp),Ztr,nBlock,dipN(L,:,:),nBlock,(0.0_wp,0.0_wp),TMP,nBlock)
  call ZGEMM_('N','N',nBlock,nBlock,nBlock,(1.0_wp,0.0_wp),TMP,nBlock,Ztr,nBlock,(0.0_wp,0.0_wp),ML(L,:,:),nBlock)
end do !L

if (iprint > 2) then
  write(6,*)
  write(6,'(10X,A)') 'Magnetic moment (ML) in the coordinate system of the magnetic axes'
  write(6,*)
  write(6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >','projection = X','projection = Y','projection = Z', &
                                         'ABS(< i | moment | j >)/3'
  write(6,*)
  do I=1,nBlock
    do J=1,nBlock
      MAVE = 0.0_wp
      MAVE = (abs(ML(1,I,J))+abs(ML(2,I,J))+abs(ML(3,I,J)))/3.0_wp
      write(6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)') '<',i,'| moment |',j,' >',(ML(L,I,J),L=1,3),MAVE
    end do
  end do
end if

call zcopy_(3*nmult*10*nmult*10,[(0.0_wp,0.0_wp)],0,DIPSO5,1)
do l=1,3
  do i1=1,nmult
    do j1=1,ndim(i1)
      do i2=1,nmult
        do j2=1,ndim(i2)
          i = ibas(i1,j1)
          j = ibas(i2,j2)
          dipso5(l,i1,j1,i2,j2) = ML(l,i,j)

          !write(6,'(A,5i3,A,2F20.12,3x,2I5)') 'DIPSO5(',l,i1,j1,i2,j2,')=',dipso5(  l,i1,j1,i2,j2), i,j
        end do
      end do
    end do
  end do
end do

if (doPLOT) then
  call plot_barrier(nBlock,nMult,nDIM,W,dipso5)
  write(6,'(A)') 'The following files'
  write(6,'(A)') '#-->  $WorkDir/$Project.BARRIER.plt'
  write(6,'(A)') '#-->  $WorkDir/$Project.BARRIER_ENE.dat'
  write(6,'(A)') '#-->  $WorkDir/$Project.BARRIER_TME.dat'
  write(6,'(A)') '#-->  $WorkDir/$Project.BARRIER.plt'
  write(6,'(A)') 'Have been generated successfully.'
end if

!cccccccccccccccccccccccccccccccccccc

!new print-out code:
write(6,*)
write(6,'(100A)') ('%',i=1,95)
write(6,'(15X,A)') 'AB INITIO BLOCKING BARRIER'
write(6,'(100A)') ('%',i=1,95)
write(6,'(A)') 'please, acknowledge the fact that the information printed below provides'
write(6,'(A)') 'only a qualitative relaxation path of a single-molecule magnet'
write(6,*)
! check the parity of all manifolds:
ipar = 0
i = 1
do il=1,nmult
  ipar = ipar+mod(ndim(il),2)
  if (mod(ndim(il),2) == 1) i = i+1
  ! count the number of odd manifolds
  !write(6,'(3(A,i2,2x))') 'il=',il,'ipar=',ipar,'i=',i
end do
if (i > 1) ipar = ipar/(i-1)

Ifunct = 0
do il=1,nmult
  E(il) = 0.0_wp
  do i=1,ndim(il)
    E(il) = E(il)+W(Ifunct+i)/ndim(il)
  end do
  Ifunct = Ifunct+ndim(il)
end do
write(6,'(A)') 'Zeeman eigenstates:'

! the convention to label states in the blocing barrier is the following:
!   Size of the             labelling scheme
!   manifold
!      1    ---------------                     0
!      2    ---------------                 1+,    1-
!      3    ---------------                 1+, 0, 1-
!      4    ---------------             2+, 1+,    1-, 2-
!      5    ---------------             2+, 1+, 0, 1-, 2-
!      6    ---------------         3+, 2+, 1+,    1-, 2-, 3-
!      7    ---------------         3+, 2+, 1+, 0, 1-, 2-, 3-
!      8    ---------------     4+, 3+, 2+, 1+,    1-, 2-, 3-, 4-
!      9    ---------------     4+, 3+, 2+, 1+, 0, 1-, 2-, 3-, 4-
!      9    --------------- 5+, 4+, 3+, 2+, 1+,    1-, 2-, 3-, 4-, 5-
!      9    --------------- 5+, 4+, 3+, 2+, 1+, 0, 1-, 2-, 3-, 4-, 5-
!    label(il,istate) is a Character*5 array of dimension (nmult,10)
!
!  notation:   multiplet . Lbl+
!  two Characters are assigned for multiplet

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (ipar == 0) then
  !  all multiplets have the same parity
  write(string2,'(a, i2, a)') '(A,',maxmult,'A,A)'
  write(6,string2) '-------',('--------------',i=1,maxmult),'----------------'
  if (mod(maxmult,2) == 0) then
    write(string2,'(a, i2, a,i2,a)') '(A,',maxmult/2,'(A,i2,a),',maxmult/2,'(A,i2,a),a)'
    write(6,string2) ' Mult.|',('     ',i,'+     |',i=maxmult/2,1,-1),('     ',i,'-     |',i=1,maxmult/2,1),'    E (cm-1)   |'
  else if (mod(maxmult,2) == 1) then
    write(string2,'(a,i2,a,i2,a)') '(A,',(maxmult-1)/2,'(A,i2,a),a,',(maxmult-1)/2,'(A,i2,a),a)'
    write(6,string2) ' Mult.|',('     ',i,'+     |',i=int((maxmult-1)/2),1,-1),'      0      |', &
                     ('     ',i,'-     |',i=1,int((maxmult-1)/2),1),'    E (cm-1)   |'
  end if
  write(string2,'(a, i2, a)') '(A,',maxmult,'A,A)'
  write(6,string2) '------|',('-------------|',i=1,maxmult),'---------------|'

  do il=1,nmult
    if (ndim(il) < maxmult) then
      nb = int((maxmult-ndim(il))/2)
      write(string2,'(3(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
      write(6,string2) il,'. | ',('            | ',i=1,nb),(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)), &
                       ('            | ',i=1,nb),E(il),' |'
    else
      write(string2,'(a, i2, a)') '(2x,i2,a,',maxmult,'(F11.7,a),F13.7,a)'
      write(6,string2) il,'. | ',(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),E(il),' |'
    end if
  end do !il
  write(string2,'(a, i2, a)') '(A,',maxmult,'A,A)'
  write(6,string2) '------|',('-------------|',i=1,maxmult),'---------------|'

else !ipar
  ! multiplets have different parity (even and odd)
  write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
  write(6,string2) '-------',('--------------',i=1,maxmult/2),'--------------',('--------------',i=1,maxmult/2),'----------------'
  if (mod(maxmult,2) == 0) then  !maxmult = even
    write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'(A,i2,a),a,',maxmult/2,'(A,i2,a),a)'
    write(6,string2) ' Mult.|',('     ',i,'+     |',i=maxmult/2,1,-1),'      0      |',('     ',i,'-     |',i=1,maxmult/2,1), &
                     '    E (cm-1)   |'
    write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
    write(6,string2) '------|',('-------------|',i=1,maxmult/2),'-------------|',('-------------|',i=1,maxmult/2),'---------------|'

    do il=1,nmult
      if (mod(ndim(il),2) == 0) then
        ! il = even > the same parity as maxmult
        nb = (maxmult-ndim(il))/2
        if (nb == 0) then !ndim(il) => maxmult
          write(string2,'(2(a,i2),a)') '(2x,i2,a,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),F13.7,a)'
          write(6,string2) il,'. | ',(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                           (dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),E(il),' |'
        else !nb>0, Integer
          write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),',nb,'A,F13.7,a)'
          write(6,string2) il,'. | ',('            | ',i=1,nb),(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                           (dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),('            | ',i=1,nb),E(il),' |'
        end if
      else ! il = odd  < maxmult
        nb = (maxmult+1-ndim(il))/2
        write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
        write(6,string2) il,'. | ',('            | ',i=1,nb),(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)), &
                         ('            | ',i=1,nb),E(il),' |'
      end if
    end do

  else ! maxmult = odd
    write(string2,'(a, i2, a, i2, a)') '(A,',(maxmult-1)/2,'(a,i2,a),a,',(maxmult-1)/2,'(a,i2,a),A)'
    write(6,string2) ' Mult.|',('     ',i,'+     |',i=int((maxmult-1)/2),1,-1),'      0      |', &
                     ('     ',i,'-     |',i=1,int((maxmult-1)/2),1),'    E (cm-1)   |'
    write(string2,'(a, i2, a, i2, a)')'(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
    write(6,string2) '------|',('-------------|',i=1,maxmult/2),'-------------|',('-------------|',i=1,maxmult/2),'---------------|'

    do il=1,nmult
      if (mod(ndim(il),2) == 1) then
        ! il = odd,  the parity of il = maxmult
        nb = (maxmult-ndim(il))/2
        if (nb == 0) then !ndim(il) => maxmult
          write(string2,'(2(a,i2),a)') '(2x,i2,a,',ndim(il),'(F11.7,a),F13.7,a)'
          write(6,string2) il,'. | ',(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),E(il),' |'
        else !nb>0, Integer
          write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
          write(6,string2) il,'. | ',('            | ',i=1,nb),(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)), &
                           ('            | ',i=1,nb),E(il),' |'
        end if

      else ! il = even  < maxmult
        nb = (maxmult-1-ndim(il))/2
        if (nb == 0) then !ndim(il) => maxmult
          write(string2,'(2(a,i2),a)') '(2x,i2,a,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),F13.7,a)'
          write(6,string2) il,'. | ',(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                           (dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),E(il),' |'
        else !nb>0, Integer
          write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),',nb,'A,F13.7,a)'
          write(6,string2) il,'. | ',('            | ',i=1,nb),(dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                           (dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),('            | ',i=1,nb),E(il),' |'
        end if
      end if
    end do

  end if
  write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
  write(6,string2) '------|',('-------------|',i=1,maxmult/2),'-------------|',('-------------|',i=1,maxmult/2),'---------------|'
end if ! ipar , line 271
write(6,*)
! printing of all off-diagonal matrix elements
write(6,'(A)') 'Matrix elements of the magnetic moment connecting Zeeman eigenstates'
write(6,'(A)') 'The average is done according to the formula:'
write(6,'(A)') '<i|m|j> = ( ABS(<i|m_X|j>) + ABS(<i|m_Y|j>) + ABS(<i|m_Z|j>) ) / 3'
write(6,*)
write(6,'(A)') 'MATRIX ELEMENTS BETWEEN STATES WITH OPPOSITE MAGNETIZATION'
write(6,'(A)') 'in cases with even number of electrons, these values cannot be used'
write(6,'(A)') 'check the tunnelling splitting instead'
write(6,'(4A)') '-------','-------------------------','----------------------------------------','------------------------'
write(6,'(A,5x,A,5x,A,13x,A,13x,A,8x,a,8x,a)') ' Mult.|','Matrix Element','|','Complex VALUE','|','AVERAGE','|'
write(6,'(4A)') '------|','------------------------|','---------Real------------imaginary-----|','-----------------------|'

do il=1,nmult
  if (ndim(il) == 1) then
    write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
    write(s2,'(i2,A1,A1,A1)') il,'.','0',' '
    call prbar(il,s1,s2,dipso5(1:3,il,1,il,1))
  else
    !if (mod(ndim(il),2) == 0) then  ! even multiplicity

    do i=1,int(ndim(il)/2)
      if (i > 1) write(6,'(4A)') '      |','------------------------|','---------------------------------------|', &
                                 '-----------------------|'
      write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
      do j=i,int(ndim(il)/2)
        write(s2,'(i2,A1,i1,A1)') il,'.',j,'+'
        !call prbar(il,s1,s2,dipso5(1:3,il,i,il,ndim(il)-i+1))
        call prbar(il,s1,s2,dipso5(1:3,il,i,il,j))
      end do

      do j=i,int(ndim(il)/2)
        write(s2,'(i2,A1,i1,A1)') il,'.',j,'-'
        !call prbar(il,s1,s2,dipso5(1:3,il,i,il,ndim(il)-i+1))
        call prbar(il,s1,s2,dipso5(1:3,il,i,il,ndim(il)-j+1))
      end do
    end do

    !else ! mod(ndim(il),2) == 1 , i.e. odd multiplicity
    !end if

  end if
  write(6,'(4A)') '------|','------------------------|','---------------------------------------|','-----------------------|'
end do !il
write(6,*)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(6,'(10A)') ('##########',i=1,10)
do k=1,nmult-1
  write(6,*)
  write(6,'(A,i2,a)') 'MATRIX ELEMENTS BETWEEN STATES ARISING FROM NEIGHBORING MULTIPLETS: I -> I+',k,' :'
  write(6,'(4A)') '-------','-------------------------','----------------------------------------','------------------------'
  write(6,'(A,5x,A,5x,A,13x,A,13x,A,8x,a,8x,a)') ' Mult.|','Matrix Element','|','Complex VALUE','|','AVERAGE','|'
  write(6,'(4A)') '------|','------------------------|','---------Real------------imaginary-----|','-----------------------|'
  do il=1,nmult-k
    if (mod(ndim(il),2) == 0) then
      do i=ndim(il)/2,1,-1
        if (mod(ndim(il+k),2) == 0) then
          do j=ndim(il+k)/2,1,-1

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
            call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,j))

          end do

          do j=1,ndim(il+k)/2

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
            call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,ndim(il+k)-j+1))

          end do

        else !ndim(il+k) is odd
          if (ndim(il+k) == 1) then

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,1))

          else

            do j=(ndim(il+k)-1)/2,1,-1
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,j))
            end do

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(1:3) = dipso5(1:3,il,i,il+k,(ndim(il+k)+1)/2)
            call prbar(il,s1,s2,MM(1:3))

            do j=1,(ndim(il+k)-1)/2
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
              MM(1:3) = dipso5(1:3,il,i,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM(1:3))
            end do

          end if ! ndim(il+k) = 1
        end if ! ndim(il+k)
        if (i > 1) write(6,'(4A)') '      |','------------------------|','---------------------------------------|', &
                                   '-----------------------|'
      end do !i

    else ! ndim(il) = odd

      if (ndim(il) == 1) then
        if (mod(ndim(il+k),2) == 0) then
          do j=ndim(il+k)/2,1,-1
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
            MM(1:3) = dipso5(1:3,il,1,il+k,j)
            call prbar(il,s1,s2,MM(1:3))
          end do

          do j=1,ndim(il+k)/2
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
            MM(1:3) = dipso5(1:3,il,1,il+k,ndim(il+k)-j+1)
            call prbar(il,s1,s2,MM(1:3))
          end do

        else !ndim(il+k) is odd
          if (ndim(il+k) == 1) then
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(1:3) = dipso5(1:3,il,1,il+k,1)
            call prbar(il,s1,s2,MM(1:3))
          else
            do j=(ndim(il+k)-1)/2,1,-1
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              MM(1:3) = dipso5(1:3,il,1,il+k,j)
              call prbar(il,s1,s2,MM(1:3))
            end do
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(1:3) = dipso5(1:3,il,1,il+k,(ndim(il+k)+1)/2)
            call prbar(il,s1,s2,MM(1:3))

            do j=1,(ndim(il+k)-1)/2
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
              MM(1:3) = dipso5(1:3,il,1,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM(1:3))
            end do
          end if ! ndim(il+k) = 1
        end if ! ndim(il+k), parity

      else !ndim(il) > 1, odd

        do i=(ndim(il)-1)/2,1,-1
          if (mod(ndim(il+k),2) == 0) then
            do j=ndim(il+k)/2,1,-1
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              MM(1:3) = dipso5(1:3,il,i,il+k,j)
              call prbar(il,s1,s2,MM(1:3))
            end do
            do j=1,ndim(il+k)/2
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
              MM(1:3) = dipso5(1:3,il,i,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM(1:3))
            end do
          else !ndim(il+k) is odd
            if (ndim(il+k) == 1) then
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
              MM(1:3) = dipso5(1:3,il,i,il+k,1)
              call prbar(il,s1,s2,MM(1:3))
            else
              do j=(ndim(il+k)-1)/2,1,-1
                write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
                write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
                MM(1:3) = dipso5(1:3,il,i,il+k,j)
                call prbar(il,s1,s2,MM(1:3))
              end do
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
              MM(1:3) = dipso5(1:3,il,i,il+k,(ndim(il+k)+1)/2)
              call prbar(il,s1,s2,MM(1:3))
              do j=1,(ndim(il+k)-1)/2
                write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
                write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
                MM(1:3) = dipso5(1:3,il,i,il+k,ndim(il+k)-j+1)
                call prbar(il,s1,s2,MM(1:3))
              end do
            end if ! ndim(il+k) = 1
          end if ! ndim(il+k), parity
          write(6,'(4A)') '      |','------------------------|','---------------------------------------|', &
                          '-----------------------|'
        end do ! i

        if (mod(ndim(il+k),2) == 0) then
          do j=ndim(il+k)/2,1,-1
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
            MM(1:3) = dipso5(1:3,il,(ndim(il)+1)/2,il+k,j)
            call prbar(il,s1,s2,MM(1:3))
          end do
          do j=1,ndim(il+k)/2
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
            MM(1:3) = dipso5(1:3,il,(ndim(il)+1)/2,il+k,ndim(il+k)-j+1)
            call prbar(il,s1,s2,MM(1:3))
          end do
        else !ndim(il+k) is odd
          if (ndim(il+k) == 1) then
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(1:3) = dipso5(1:3,il,(ndim(il)+1)/2,il+k,1)
            call prbar(il,s1,s2,MM(1:3))
          else
            do j=(ndim(il+k)-1)/2,1,-1
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              MM(1:3) = dipso5(1:3,il,(ndim(il)+1)/2,il+k,j)
              call prbar(il,s1,s2,MM(1:3))
            end do
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(1:3) = dipso5(1:3,il,(ndim(il)+1)/2,il+k,(ndim(il+k)+1)/2)
            call prbar(il,s1,s2,MM(1:3))
            do j=1,(ndim(il+k)-1)/2
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
              MM(1:3) = dipso5(1:3,il,(ndim(il)+1)/2,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM(1:3))
            end do
          end if ! ndim(il+k) = 1
        end if ! ndim(il+k), parity

      end if ! ndim(il), size

    end if !ndim(il), parity
    write(6,'(4A)') '------|','------------------------|','---------------------------------------|','-----------------------|'
  end do ! il
end do ! k

!-----------------------------------------------------------------------
if ((nmult > 0) .and. (nBlock > 0)) then
  call mma_deallocate(ibas)
  call mma_deallocate(wz)
  call mma_deallocate(cz)
end if

if (nmult > 0) then
  call mma_deallocate(E)
  call mma_deallocate(dipso5)
end if

if (nBlock > 0) then
  call mma_deallocate(dipN)
  call mma_deallocate(ML)
  call mma_deallocate(Ztr)
  call mma_deallocate(tmp)
end if

k = nDim(imanifold)
if (k > 0) call mma_deallocate(dipN3)

call mma_deallocate(gtens)
call mma_deallocate(MM)
call mma_deallocate(maxes)

return

end subroutine barrier
