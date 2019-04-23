************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine natureorb(nsbas,nsall,nsdel,ngsm,den1,lden,cmo,lcmo,
     *                     bsbl,lenb,cno,occ,nmo,pror)
      implicit none
      integer :: nsbas(8),nsall(8),nsdel(8)
      integer :: ngsm,lden,lcmo,lenb,nmo
      real*8 :: cmo(lcmo),cno(lcmo),occ(nmo)
      real*8 :: den1(lden)
      character bsbl(lenb)*1
      real*8 :: pror
!-----------------------------------------------------
      integer :: nsfrz(8),nsort(nmo)
      integer :: nc,nc0,nc1,nc2,nc3,nc4,nc5
      integer :: im,i,j
      real*8 :: buff(nmo**2)
      real*8 :: val
      character*128 :: header

      header="MRCISD Natrual orbital"
      nc0=1
      do im=1,ngsm
        if(nsall(im).eq.0) cycle
        nc1=nsall(im)*(nsall(im)+1)/2
        nc0=nc0+nc1
      enddo
!      open(100,file="dat2")
!      open(200,file="dat1")
!      do i=1,1293
!        read(100,*) den1(i)
!        write(200,"(f14.8)") den1(i)
!      enddo
!      close(100)
!      close(200)


      val=0
      cno=cmo
      nc0=1
      nc1=1
      nc2=1
      do im=1,ngsm
        if(nsbas(im).eq.0) cycle
        nsfrz(im)=nsbas(im)-nsall(im)-nsdel(im)
! For density matrix, we only have correlated orbital density matrix
        buff=0.d0
        nc=nsall(im)*(nsall(im)+1)/2
! Dagnolize density matrix in symmetry block im and transform MO
! Copy density matrix
        buff(1:nc)=den1(nc2:nc2+nc-1)
        nc3=nc1+nsfrz(im)*nsbas(im)
        call jacob(buff,cno(nc3),nsall(im),nsbas(im))
! OCC num from diagonal element
        nc3=nc0
        occ(nc3:nc3+nsfrz(im)-1)=2.d0
        nc3=nc3+nsfrz(im)
        do i=1,nsall(im)
          nc=i*(i+1)/2
          occ(nc3)=buff(nc)
          nc3=nc3+1
        enddo

! Sort nature orbital in OCC num decreasing order
        buff=0.d0
! Sort and copy correlated nature orb
        do i=1,nsall(im)
          nsort(i)=i
        enddo
        nc=nc0+nsfrz(im)-1
        do i=1,nsall(im)
          nc3=i
          do j=i+1,nsall(im)
            if(occ(nc+j).gt.occ(nc+i)) then
              val=occ(nc+j)
              occ(nc+j)=occ(nc+i)
              occ(nc+i)=val
              nc3=j
            endif
          enddo
          nsort(i)=nc3
        enddo

        buff=0.d0
        nc3=nsbas(im)**2
        buff(1:nc3)=cno(nc1:nc1+nc3-1)
! Copy sorted orbital
        nc=nc1+nsfrz(im)*nsbas(im)
        nc5=1+nsfrz(im)*nsbas(im)
        do i=1,nsall(im)
          j=nsort(i)
          nc3=nc+(i-1)*nsbas(im)
          nc4=nc5+(j-1)*nsbas(im)
          cno(nc3:nc3+nsbas(im)-1)=buff(nc4:nc4+nsbas(im)-1)
        enddo

! do nothing for deleted orbital
!        nc3=1+(nsfrz(im)+nsall(im))*nsbas(im)
!        nc=nsdel(im)*nsbas(im)-1
!        cno(nc3:nc3+nc)=buff(nc3:nc3+nc)

!100     continue
        nc0=nc0+nsbas(im)
        nc1=nc1+nsbas(im)**2
        nc2=nc2+nsall(im)*(nsall(im)+1)/2
      enddo

!      nc0=1
!      do im=1,ngsm
!        if(nsbas(im).eq.0) cycle
!        print*, " "
!        print*, "SYM",im
!        print"(10(1xf8.4))", occ(nc0:nc0+nsbas(im)-1)
!        nc0=nc0+nsbas(im)
!      enddo

      call primo(Header,.true.,.false.,pror,0.d0,ngsm,nsbas,nsbas,
     *           bsbl,[val],occ,cno,-1)

      return
      end
