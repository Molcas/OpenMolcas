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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBlock, imanifold, nMult, nDim(nMult), iprint
complex(kind=wp), intent(in) :: dipIn(3,nBlock,nBlock)
real(kind=wp), intent(in) :: W(nBlock)
logical(kind=iwp), intent(in) :: doPLOT
integer(kind=iwp) :: dimi, i, i1, i2, Ifunct, il, iMult, ipar, j, j1, j2, k, l, maxmult, nb
real(kind=wp) :: gtens(3), mave, maxes(3,3)
complex(kind=wp) :: MM(3)
character(len=90) :: string2
character(len=5) :: s1, s2
integer(kind=iwp), allocatable :: ibas(:,:)
real(kind=wp), allocatable :: E(:), wz(:,:)
complex(kind=wp), allocatable :: CZ(:,:,:), dipN(:,:,:), dipN3(:,:,:), dipso5(:,:,:,:,:), ML(:,:,:), tmp(:,:), tmp2(:,:), &
                                 tmp3(:,:,:), Ztr(:,:)

if ((nmult > 0) .and. (nBlock > 0)) then
  call mma_allocate(ibas,nmult,nBlock,'ibas')
  call mma_allocate(wz,nmult,nBlock,'wz')
  call mma_allocate(CZ,nmult,nBlock,nBlock,'CZ')
  ibas(:,:) = 0
  wz(:,:) = Zero
  cz(:,:,:) = cZero
end if

if (nmult > 0) then
  call mma_allocate(E,nmult,'E')
  call mma_allocate(dipso5,3,nmult,10,nmult,10,'dipso5')
  E(:) = Zero
  dipso5(:,:,:,:,:) = cZero
end if

if (nBlock > 0) then
  call mma_allocate(dipN,3,nBlock,nBlock,'dipN')
  call mma_allocate(ML,3,nBlock,nBlock,'ML')
  call mma_allocate(Ztr,nBlock,nBlock,'Ztr')
  call mma_allocate(tmp,nBlock,nBlock,'tmp')
  dipN(:,:,:) = cZero
  ML(:,:,:) = cZero
  Ztr(:,:) = cZero
  tmp(:,:) = cZero
end if

k = nDim(imanifold)
if (k > 0) then
  call mma_allocate(dipN3,3,k,k,'dipN3')
  dipN3(:,:,:) = cZero
end if

gtens(:) = Zero
MM(:) = Zero

!-----------------------------------------------------------------------
call unitmat(maxes,3)
! rotate the magnetic moment to the magnetic axes of the ground multiplet ( NDIM(1) )
if (ndim(imanifold) <= 1) then
  write(u6,'(a,i2,a)') 'The manifold ',imanifold,' was chosen for determination of the quantization axis.'
  write(u6,'(a     )') 'However, the size of this manifold is:'
  write(u6,'(a,i1  )') 'size = ',ndim(imanifold)
  write(u6,'(a     )') 'in this case, quantization axis will remain the original Z axis'
else
  dipN3(:,:,:) = dipIn(:,:,:)
  call atens(dipN3,ndim(imanifold),gtens,maxes,1)
end if
call rotmom2(dipIn,nBlock,maxes,dipN)
if (iprint > 2) then
  write(u6,*)
  write(u6,'(10X,A)') 'Magnetic moment (dipIN) in the original coordinate system'
  write(u6,*)
  write(u6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >','projection = X','projection = Y','projection = Z', &
                                          'ABS(< i | moment | j >)/3'
  write(u6,*)
  do I=1,nBlock
    do J=1,nBlock
      MAVE = (abs(dipIN(1,I,J))+abs(dipIN(2,I,J))+abs(dipIN(3,I,J)))/Three
      write(u6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)') '<',i,'| moment |',j,' >',(dipIN(L,I,J),L=1,3),MAVE
    end do
  end do
  write(u6,'(///)')
  write(u6,'(10X,A)') 'Magnetic moment (dipN) in the coordinate system of the magnetic axes'
  write(u6,*)
  write(u6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >','projection = X','projection = Y','projection = Z', &
                                          'ABS(< i | moment | j >)/3'
  write(u6,*)
  do I=1,nBlock
    do J=1,nBlock
      MAVE = (abs(dipN(1,I,J))+abs(dipN(2,I,J))+abs(dipN(3,I,J)))/Three
      write(u6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)') '<',i,'| moment |',j,' >',(dipN(L,I,J),L=1,3),MAVE
    end do
  end do
end if
!-----------------------------------------------------------------------
! determine the "CZ matrix":
Ifunct = 1
do il=1,nmult
  dimi = ndim(il)
  if (dimi == 0) return
  if (dimi == 1) then
    CZ(il,1,1) = cOne
    WZ(il,1) = W(Ifunct)
  else ! dimi > 1

    call mma_allocate(tmp2,ndim(il),ndim(il),label='tmp2')
    call mma_allocate(tmp3,3,ndim(il),ndim(il),label='tmp3')
    tmp3(:,:,:) = dipN(:,Ifunct:Ifunct+ndim(il)-1,Ifunct:Ifunct+ndim(il)-1)
    call pseudospin(tmp3,ndim(il),tmp2,3,1,iprint)
    CZ(il,1:ndim(il),1:ndim(il)) = tmp2(:,:)
    call mma_deallocate(tmp2)
    call mma_deallocate(tmp3)

    if (iPrint > 2) call pa_prMat('barrier:  CZ',CZ(il,1:nDim(il),1:nDim(il)),nDim(il))
  end if
  Ifunct = Ifunct+ndim(il)
end do
!-----------------------------------------------------------------------
! compute the matrix elements between multiplets
l = 0
ibas(:,:) = 0
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
do iMult=1,nmult
  do j1=1,ndim(iMult)
    do j2=1,ndim(iMult)
      i = ibas(iMult,j1)
      j = ibas(iMult,j2)
      Ztr(i,j) = CZ(iMult,j1,j2)
    end do
  end do
end do

call mma_allocate(tmp2,nBlock,nblock,label='tmp2')
do L=1,3
  tmp2(:,:) = dipN(L,:,:)
  call ZGEMM_('C','N',nBlock,nBlock,nBlock,cOne,Ztr,nBlock,tmp2,nBlock,cZero,TMP,nBlock)
  call ZGEMM_('N','N',nBlock,nBlock,nBlock,cOne,TMP,nBlock,Ztr,nBlock,cZero,tmp2,nBlock)
  ML(L,:,:) = tmp2(:,:)
end do !L
call mma_deallocate(tmp2)

if (iprint > 2) then
  write(u6,*)
  write(u6,'(10X,A)') 'Magnetic moment (ML) in the coordinate system of the magnetic axes'
  write(u6,*)
  write(u6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >','projection = X','projection = Y','projection = Z', &
                                          'ABS(< i | moment | j >)/3'
  write(u6,*)
  do I=1,nBlock
    do J=1,nBlock
      MAVE = (abs(ML(1,I,J))+abs(ML(2,I,J))+abs(ML(3,I,J)))/Three
      write(u6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)') '<',i,'| moment |',j,' >',(ML(L,I,J),L=1,3),MAVE
    end do
  end do
end if

do l=1,3
  do i1=1,nmult
    do j1=1,ndim(i1)
      do i2=1,nmult
        do j2=1,ndim(i2)
          i = ibas(i1,j1)
          j = ibas(i2,j2)
          dipso5(l,i1,j1,i2,j2) = ML(l,i,j)

          !write(u6,'(A,5i3,A,2F20.12,3x,2I5)') 'DIPSO5(',l,i1,j1,i2,j2,')=',dipso5(  l,i1,j1,i2,j2), i,j
        end do
      end do
    end do
  end do
end do

if (doPLOT) then
  call plot_barrier(nBlock,nMult,nDIM,W,dipso5)
  write(u6,'(A)') 'The following files'
  write(u6,'(A)') '#-->  $WorkDir/$Project.BARRIER.plt'
  write(u6,'(A)') '#-->  $WorkDir/$Project.BARRIER_ENE.dat'
  write(u6,'(A)') '#-->  $WorkDir/$Project.BARRIER_TME.dat'
  write(u6,'(A)') '#-->  $WorkDir/$Project.BARRIER.plt'
  write(u6,'(A)') 'Have been generated successfully.'
end if

!cccccccccccccccccccccccccccccccccccc

!new print-out code:
write(u6,*)
write(u6,'(A)') repeat('%',95)
write(u6,'(15X,A)') 'AB INITIO BLOCKING BARRIER'
write(u6,'(A)') repeat('%',95)
write(u6,'(A)') 'please, acknowledge the fact that the information printed below provides'
write(u6,'(A)') 'only a qualitative relaxation path of a single-molecule magnet'
write(u6,*)
! check the parity of all manifolds:
ipar = 0
i = 1
do il=1,nmult
  ipar = ipar+mod(ndim(il),2)
  if (mod(ndim(il),2) == 1) i = i+1
  ! count the number of odd manifolds
  !write(u6,'(3(A,i2,2x))') 'il=',il,'ipar=',ipar,'i=',i
end do
if (i > 1) ipar = ipar/(i-1)

Ifunct = 0
do il=1,nmult
  do i=1,ndim(il)
    E(il) = E(il)+W(Ifunct+i)/ndim(il)
  end do
  Ifunct = Ifunct+ndim(il)
end do
write(u6,'(A)') 'Zeeman eigenstates:'

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
!    label(il,istate) is a Character(len=5) array of dimension (nmult,10)
!
!  notation:   multiplet . Lbl+
!  two Characters are assigned for multiplet

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (ipar == 0) then
  !  all multiplets have the same parity
  write(string2,'(a, i2, a)') '(A,',maxmult,'A,A)'
  write(u6,string2) '-------',('--------------',i=1,maxmult),'----------------'
  if (mod(maxmult,2) == 0) then
    write(string2,'(a, i2, a,i2,a)') '(A,',maxmult/2,'(A,i2,a),',maxmult/2,'(A,i2,a),a)'
    write(u6,string2) ' Mult.|',('     ',i,'+     |',i=maxmult/2,1,-1),('     ',i,'-     |',i=1,maxmult/2,1),'    E (cm-1)   |'
  else if (mod(maxmult,2) == 1) then
    write(string2,'(a,i2,a,i2,a)') '(A,',(maxmult-1)/2,'(A,i2,a),a,',(maxmult-1)/2,'(A,i2,a),a)'
    write(u6,string2) ' Mult.|',('     ',i,'+     |',i=int((maxmult-1)/2),1,-1),'      0      |', &
                      ('     ',i,'-     |',i=1,int((maxmult-1)/2),1),'    E (cm-1)   |'
  end if
  write(string2,'(a, i2, a)') '(A,',maxmult,'A,A)'
  write(u6,string2) '------|',('-------------|',i=1,maxmult),'---------------|'

  do il=1,nmult
    if (ndim(il) < maxmult) then
      nb = int((maxmult-ndim(il))/2)
      write(string2,'(3(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
      write(u6,string2) il,'. | ',('            | ',i=1,nb),(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)), &
                        ('            | ',i=1,nb),E(il),' |'
    else
      write(string2,'(a, i2, a)') '(2x,i2,a,',maxmult,'(F11.7,a),F13.7,a)'
      write(u6,string2) il,'. | ',(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),E(il),' |'
    end if
  end do !il
  write(string2,'(a, i2, a)') '(A,',maxmult,'A,A)'
  write(u6,string2) '------|',('-------------|',i=1,maxmult),'---------------|'

else !ipar
  ! multiplets have different parity (even and odd)
  write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
  write(u6,string2) '-------',('--------------',i=1,maxmult/2),'--------------',('--------------',i=1,maxmult/2),'----------------'
  if (mod(maxmult,2) == 0) then  !maxmult = even
    write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'(A,i2,a),a,',maxmult/2,'(A,i2,a),a)'
    write(u6,string2) ' Mult.|',('     ',i,'+     |',i=maxmult/2,1,-1),'      0      |',('     ',i,'-     |',i=1,maxmult/2,1), &
                      '    E (cm-1)   |'
    write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
    write(u6,string2) '------|',('-------------|',i=1,maxmult/2),'-------------|',('-------------|',i=1,maxmult/2), &
                      '---------------|'

    do il=1,nmult
      if (mod(ndim(il),2) == 0) then
        ! il = even > the same parity as maxmult
        nb = (maxmult-ndim(il))/2
        if (nb == 0) then !ndim(il) => maxmult
          write(string2,'(2(a,i2),a)') '(2x,i2,a,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),F13.7,a)'
          write(u6,string2) il,'. | ',(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                            (real(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),E(il),' |'
        else !nb>0, Integer
          write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),',nb,'A,F13.7,a)'
          write(u6,string2) il,'. | ',('            | ',i=1,nb),(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                            (real(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),('            | ',i=1,nb),E(il),' |'
        end if
      else ! il = odd  < maxmult
        nb = (maxmult+1-ndim(il))/2
        write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
        write(u6,string2) il,'. | ',('            | ',i=1,nb),(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)), &
                          ('            | ',i=1,nb),E(il),' |'
      end if
    end do

  else ! maxmult = odd
    write(string2,'(a, i2, a, i2, a)') '(A,',(maxmult-1)/2,'(a,i2,a),a,',(maxmult-1)/2,'(a,i2,a),A)'
    write(u6,string2) ' Mult.|',('     ',i,'+     |',i=(maxmult-1)/2,1,-1),'      0      |', &
                      ('     ',i,'-     |',i=1,(maxmult-1)/2,1),'    E (cm-1)   |'
    write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
    write(u6,string2) '------|',('-------------|',i=1,maxmult/2),'-------------|',('-------------|',i=1,maxmult/2), &
                      '---------------|'

    do il=1,nmult
      if (mod(ndim(il),2) == 1) then
        ! il = odd,  the parity of il = maxmult
        nb = (maxmult-ndim(il))/2
        if (nb == 0) then !ndim(il) => maxmult
          write(string2,'(2(a,i2),a)') '(2x,i2,a,',ndim(il),'(F11.7,a),F13.7,a)'
          write(u6,string2) il,'. | ',(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),E(il),' |'
        else !nb>0, Integer
          write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
          write(u6,string2) il,'. | ',('            | ',i=1,nb),(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)), &
                            ('            | ',i=1,nb),E(il),' |'
        end if

      else ! il = even  < maxmult
        nb = (maxmult-1-ndim(il))/2
        if (nb == 0) then !ndim(il) => maxmult
          write(string2,'(2(a,i2),a)') '(2x,i2,a,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),F13.7,a)'
          write(u6,string2) il,'. | ',(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                            (real(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),E(il),' |'
        else !nb>0, Integer
          write(string2,'(4(a,i2),a)') '(2x,i2,a,',nb,'A,',ndim(il)/2,'(F11.7,a),    A,',ndim(il)/2,'(F11.7,a),',nb,'A,F13.7,a)'
          write(u6,string2) il,'. | ',('            | ',i=1,nb),(real(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),'            | ', &
                            (real(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,ndim(il)),('            | ',i=1,nb),E(il),' |'
        end if
      end if
    end do

  end if
  write(string2,'(a, i2, a, i2, a)') '(A,',maxmult/2,'A,A,',maxmult/2,'A,A)'
  write(u6,string2) '------|',('-------------|',i=1,maxmult/2),'-------------|',('-------------|',i=1,maxmult/2),'---------------|'
end if ! ipar , line 271
write(u6,*)
! printing of all off-diagonal matrix elements
write(u6,'(A)') 'Matrix elements of the magnetic moment connecting Zeeman eigenstates'
write(u6,'(A)') 'The average is done according to the formula:'
write(u6,'(A)') '<i|m|j> = ( ABS(<i|m_X|j>) + ABS(<i|m_Y|j>) + ABS(<i|m_Z|j>) ) / 3'
write(u6,*)
write(u6,'(A)') 'MATRIX ELEMENTS BETWEEN STATES WITH OPPOSITE MAGNETIZATION'
write(u6,'(A)') 'in cases with even number of electrons, these values cannot be used'
write(u6,'(A)') 'check the tunnelling splitting instead'
write(u6,'(4A)') '-------','-------------------------','----------------------------------------','------------------------'
write(u6,'(A,5x,A,5x,A,13x,A,13x,A,8x,a,8x,a)') ' Mult.|','Matrix Element','|','Complex VALUE','|','AVERAGE','|'
write(u6,'(4A)') '------|','------------------------|','---------Real------------imaginary-----|','-----------------------|'

do il=1,nmult
  if (ndim(il) == 1) then
    write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
    write(s2,'(i2,A1,A1,A1)') il,'.','0',' '
    call prbar(il,s1,s2,dipso5(:,il,1,il,1))
  else
    !if (mod(ndim(il),2) == 0) then  ! even multiplicity

    do i=1,ndim(il)/2
      if (i > 1) write(u6,'(4A)') '      |','------------------------|','---------------------------------------|', &
                                  '-----------------------|'
      write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
      do j=i,ndim(il)/2
        write(s2,'(i2,A1,i1,A1)') il,'.',j,'+'
        !call prbar(il,s1,s2,dipso5(:,il,i,il,ndim(il)-i+1))
        call prbar(il,s1,s2,dipso5(:,il,i,il,j))
      end do

      do j=i,ndim(il)/2
        write(s2,'(i2,A1,i1,A1)') il,'.',j,'-'
        !call prbar(il,s1,s2,dipso5(:,il,i,il,ndim(il)-i+1))
        call prbar(il,s1,s2,dipso5(:,il,i,il,ndim(il)-j+1))
      end do
    end do

    !else ! mod(ndim(il),2) == 1 , i.e. odd multiplicity
    !end if

  end if
  write(u6,'(4A)') '------|','------------------------|','---------------------------------------|','-----------------------|'
end do !il
write(u6,*)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(u6,'(10A)') ('##########',i=1,10)
do k=1,nmult-1
  write(u6,*)
  write(u6,'(A,i2,a)') 'MATRIX ELEMENTS BETWEEN STATES ARISING FROM NEIGHBORING MULTIPLETS: I -> I+',k,' :'
  write(u6,'(4A)') '-------','-------------------------','----------------------------------------','------------------------'
  write(u6,'(A,5x,A,5x,A,13x,A,13x,A,8x,a,8x,a)') ' Mult.|','Matrix Element','|','Complex VALUE','|','AVERAGE','|'
  write(u6,'(4A)') '------|','------------------------|','---------Real------------imaginary-----|','-----------------------|'
  do il=1,nmult-k
    if (mod(ndim(il),2) == 0) then
      do i=ndim(il)/2,1,-1
        if (mod(ndim(il+k),2) == 0) then
          do j=ndim(il+k)/2,1,-1

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
            call prbar(il,s1,s2,dipso5(:,il,i,il+k,j))

          end do

          do j=1,ndim(il+k)/2

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
            call prbar(il,s1,s2,dipso5(:,il,i,il+k,ndim(il+k)-j+1))

          end do

        else !ndim(il+k) is odd
          if (ndim(il+k) == 1) then

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            call prbar(il,s1,s2,dipso5(:,il,i,il+k,1))

          else

            do j=(ndim(il+k)-1)/2,1,-1
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              call prbar(il,s1,s2,dipso5(:,il,i,il+k,j))
            end do

            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(:) = dipso5(:,il,i,il+k,(ndim(il+k)+1)/2)
            call prbar(il,s1,s2,MM)

            do j=1,(ndim(il+k)-1)/2
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
              MM(:) = dipso5(:,il,i,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM)
            end do

          end if ! ndim(il+k) = 1
        end if ! ndim(il+k)
        if (i > 1) write(u6,'(4A)') '      |','------------------------|','---------------------------------------|', &
                                    '-----------------------|'
      end do !i

    else ! ndim(il) = odd

      if (ndim(il) == 1) then
        if (mod(ndim(il+k),2) == 0) then
          do j=ndim(il+k)/2,1,-1
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
            MM(:) = dipso5(:,il,1,il+k,j)
            call prbar(il,s1,s2,MM)
          end do

          do j=1,ndim(il+k)/2
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
            MM(:) = dipso5(:,il,1,il+k,ndim(il+k)-j+1)
            call prbar(il,s1,s2,MM)
          end do

        else !ndim(il+k) is odd
          if (ndim(il+k) == 1) then
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(:) = dipso5(:,il,1,il+k,1)
            call prbar(il,s1,s2,MM)
          else
            do j=(ndim(il+k)-1)/2,1,-1
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              MM(:) = dipso5(:,il,1,il+k,j)
              call prbar(il,s1,s2,MM)
            end do
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(:) = dipso5(:,il,1,il+k,(ndim(il+k)+1)/2)
            call prbar(il,s1,s2,MM)

            do j=1,(ndim(il+k)-1)/2
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
              MM(:) = dipso5(:,il,1,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM)
            end do
          end if ! ndim(il+k) = 1
        end if ! ndim(il+k), parity

      else !ndim(il) > 1, odd

        do i=(ndim(il)-1)/2,1,-1
          if (mod(ndim(il+k),2) == 0) then
            do j=ndim(il+k)/2,1,-1
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              MM(:) = dipso5(:,il,i,il+k,j)
              call prbar(il,s1,s2,MM)
            end do
            do j=1,ndim(il+k)/2
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
              MM(:) = dipso5(:,il,i,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM)
            end do
          else !ndim(il+k) is odd
            if (ndim(il+k) == 1) then
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
              MM(:) = dipso5(:,il,i,il+k,1)
              call prbar(il,s1,s2,MM)
            else
              do j=(ndim(il+k)-1)/2,1,-1
                write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
                write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
                MM(:) = dipso5(:,il,i,il+k,j)
                call prbar(il,s1,s2,MM)
              end do
              write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
              write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
              MM(:) = dipso5(:,il,i,il+k,(ndim(il+k)+1)/2)
              call prbar(il,s1,s2,MM)
              do j=1,(ndim(il+k)-1)/2
                write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
                write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
                MM(:) = dipso5(:,il,i,il+k,ndim(il+k)-j+1)
                call prbar(il,s1,s2,MM)
              end do
            end if ! ndim(il+k) = 1
          end if ! ndim(il+k), parity
          write(u6,'(4A)') '      |','------------------------|','---------------------------------------|', &
                           '-----------------------|'
        end do ! i

        if (mod(ndim(il+k),2) == 0) then
          do j=ndim(il+k)/2,1,-1
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
            MM(:) = dipso5(:,il,(ndim(il)+1)/2,il+k,j)
            call prbar(il,s1,s2,MM)
          end do
          do j=1,ndim(il+k)/2
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
            MM(:) = dipso5(:,il,(ndim(il)+1)/2,il+k,ndim(il+k)-j+1)
            call prbar(il,s1,s2,MM)
          end do
        else !ndim(il+k) is odd
          if (ndim(il+k) == 1) then
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(:) = dipso5(:,il,(ndim(il)+1)/2,il+k,1)
            call prbar(il,s1,s2,MM)
          else
            do j=(ndim(il+k)-1)/2,1,-1
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
              MM(:) = dipso5(:,il,(ndim(il)+1)/2,il+k,j)
              call prbar(il,s1,s2,MM)
            end do
            write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
            write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
            MM(:) = dipso5(:,il,(ndim(il)+1)/2,il+k,(ndim(il+k)+1)/2)
            call prbar(il,s1,s2,MM)
            do j=1,(ndim(il+k)-1)/2
              write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
              write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
              MM(:) = dipso5(:,il,(ndim(il)+1)/2,il+k,ndim(il+k)-j+1)
              call prbar(il,s1,s2,MM)
            end do
          end if ! ndim(il+k) = 1
        end if ! ndim(il+k), parity

      end if ! ndim(il), size

    end if !ndim(il), parity
    write(u6,'(4A)') '------|','------------------------|','---------------------------------------|','-----------------------|'
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

return

end subroutine barrier
