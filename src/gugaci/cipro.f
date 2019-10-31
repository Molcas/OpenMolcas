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
      subroutine cipro()
#include "drt_h.fh"
#include "grad_h.fh"
#include "files_gugaci.fh"
      parameter (maxmolcasorb=5000,maxpro=50)
      dimension idx_idisk0(64),idx_idisk1(max_root+1)
      dimension cmo(max_orb**2),cno(max_orb**2),occ(max_orb)
!     *          denao(max_orb,max_orb)
      dimension ipcom(maxpro)
      REAL*8, pointer :: omat(:),denao(:),vprop(:,:,:)
      character bsbl(2*4*maxmolcasorb)*1
      character*8 label
      character*8 pname(maxpro),ptyp(maxpro)
      dimension pnuc(maxpro)
      dimension pgauge(3,maxpro)
      dimension idummy(1)

! mrci nature orbital
      idisk=0
      call idafile(lucimo,2,idx_idisk0,64,idisk)
!      print*, "read head"
      idisk=idx_idisk0(2)
      nc=2*4*maxmolcasorb
      call cdafile(lucimo,2,bsbl,nc,idisk)

!      open(100,file="dat1")
      nc0=0
      nc1=0
      nc2=0
      nmo=0
      do im=1,ng_sm
        nc0=nc0+nlsm_bas(im)**2
        nc1=nc1+nlsm_all(im)*(nlsm_all(im)+1)/2
        nc2=nc2+nlsm_bas(im)*(nlsm_bas(im)+1)/2
        nmo=nmo+nlsm_bas(im)
      enddo
! read property labels
      iopt=8
      npro=0
      do i=1,100
        label="undef"
        call irdone(irec,iopt+1,label,ipc,idummy,isymlb)
        if(irec.ne.0) exit
        iopt=16
        if(mod(isymlb,2).eq.0) cycle
        npro=npro+1
        pname(npro)=label
        ipcom(npro)=ipc
        ptyp(npro)="HERM"
        if(label.eq."VELOCITY") ptyp(npro)="ANTI"
        if(label.eq."ANGMOM  ") ptyp(npro)="ANTI"
        if(label(1:5).eq."MLTPV") ptyp(npro)="ANTI"
        if(npro.ge.maxpro) exit
      enddo

!      print"(10(1X,a8))", pname(1:npro)
!      print*, ipcom(1:npro)
!      print*, ptyp(1:npro)

      allocate(omat(nc2))
      allocate(denao(nc0))
      allocate(vprop(mroot,mroot,npro))
!      allocate(denao(nmo,nmo))
      idisk=idx_idisk0(3)
      call ddafile(lucimo,2,cmo,nc0,idisk)
! read overlap matirx
      call rdone(irec,6,"MLTPL  0",1,omat,idummy(1))
      idisk=0
      call idafile(luciden,2,idx_idisk1,max_root+1,idisk)

      icall=0
      do iroot=1,mroot
! read density matrix
        idisk=idx_idisk1(iroot)
        call ddafile(luciden,2,denm1,nc1,idisk)
        call natureorb(nlsm_bas,nlsm_all,nlsm_del,ng_sm,denm1,
     *                 nc1,cmo,nc0,bsbl,nc,cno,occ,nmo,pror)

! mulliken population
        call xflush(6)
        write(6,'(a,i2)')' mulliken charges for state nr ',iroot
!        call charge(nsym,nbas,name,cno,occ,smat,2,.true.,.true.)
        call charge(ng_sm,nlsm_bas,bsbl,cno,occ,omat,2,.true.,.false.)
        write(6,*)' ',('*',i=1,70)
        call xflush(6)
! transform mo density matrix to ao basis
!        print"(10i8)", nc0,nmo,nlsm_bas
        call transden(ng_sm,nlsm_bas,denao,cno,nc0,occ,nmo)

! calculated properties
        call calprop(ng_sm,nlsm_bas,mroot,iroot,iroot,nc2,npro,
     *               pname,ipcom,ptyp,
     *               denao,nc0,vprop,pgauge,pnuc,icall)
! print property
      enddo
!      close(100)


! this code copy from molcas
      if(npro.gt.0) then
c write expectation values:
        write(6,*)
        write(6,*)' expectation values of various operators:'
        write(6,*)
     *  '(note: electronic multipoles include a negative sign.)'
        do 110 iprop=1,npro
          if(ptyp(iprop).eq.'ANTI') goto 110
          do 105 ista=1,mroot,4
            iend=min(ista+3,mroot)
            write(6,*)
            write(6,'(1x,a,a8,a,i4)')
     *    '   property: ',pname(iprop),
     *    '   component:',ipcom(iprop)
            write(6,'(1x,a,3f16.8)')
     *    '    gauge origin:',(pgauge(i,iprop),i=1,3)
            write(6,'(1x,a,i8,3i16)')
     *    '            root:',(i,i=ista,iend)
            write(6,'(1x,a,4f16.8)')
     *    '      electronic:',(vprop(i,i,iprop),i=ista,iend)
            write(6,'(1x,a,4f16.8)')
     *    '         nuclear:',(pnuc(iprop),i=ista,iend)
            write(6,'(1x,a,4f16.8)')
     *    '           total:',(pnuc(iprop)+vprop(i,i,iprop),i=ista,iend)
105       continue
110     continue
        write(6,*)
      end if


      iopt=8
      npro=0
      do i=1,100
        label="undef"
        call irdone(irec,iopt+1,label,ipc,idummy,isymlb)
        if(irec.ne.0) exit
        iopt=16
!       if(mod(isymlb,2).eq.0) cycle
        npro=npro+1
        pname(npro)=label
        ipcom(npro)=ipc
        ptyp(npro)="HERM"
        if(label.eq."VELOCITY") ptyp(npro)="ANTI"
        if(label.eq."ANGMOM  ") ptyp(npro)="ANTI"
        if(label(1:5).eq."MLTPV") ptyp(npro)="ANTI"
        if(npro.ge.maxpro) exit
      enddo

!      print"(10(1x,a8))", pname(1:npro)
!      print"(10i4)", ipcom(1:npro)
!      print"(10(1x,a8))", ptyp(1:npro)

      nsiz=0
      Call Molcas_BinaryOpen_Vanilla(110,"soint.dat")
      do i=1,npro
        if(pname(i)(1:4).ne."AMFI") cycle
        omat=0.d0
        call irdone(irtc,1,pname(i),ipcom(i),idummy,isymlb)
        if (irtc.eq.0) nsiz=idummy(1)
        call rdone(irtc,0,pname(i),ipcom(i),omat,isymlb)
        if(nsiz.gt.nc2) then
          write(6,*) "in subroutine cipro, read so int error"
          call abend
c          stop 1999
        endif
        write(6,"(a8,1x,2i4)") "nsiz=",i,nsiz
        write(6,"(5(1x,f12.8))") omat(1:nsiz)
        write(110) nsiz
        write(110) omat(1:nsiz)
      enddo
      close(110)
      deallocate(omat)
      deallocate(denao)
      deallocate(vprop)
      return
      end

      subroutine calprop(ngsm,nsbas,mroot,istate,jstate,nsi,npro,pname,
     *                   ipcom,ptyp,aden,lmo,
     *                   vprop,pgauge,pnuc,icall)
      implicit none
!-----------------------------------------------------
      integer ngsm,mroot,istate,jstate,nsi,npro,lmo,icall
      character*8 pname(npro),ptyp(npro)
      integer :: ipcom(npro),nsbas(ngsm)
      real*8 :: smat(nsi),amat(nsi),aden(lmo)
      real*8 :: vprop(mroot,mroot,npro)
      real*8 :: pgauge(3,npro),pnuc(npro)
!-----------------------------------------------------
      integer nsiz
      integer i,j,im,nc,nc0,nc1,irtc,isymlb
      real*8 pint(nsi+4)
      real*8 val,sgn,ddot_
      integer idummy(1)

! we have two kind of densitry matrix, symm or anti-symm
! compress density matrix
      nc0=0
      nc1=1
      nc=0
      do im=1,ngsm
        if(nsbas(im).eq.0) cycle
        do i=1,nsbas(im)
          do j=1,i-1
            smat(nc1)=aden(nc0+(j-1)*nsbas(im)+i)
     *               +aden(nc0+(i-1)*nsbas(im)+j)
            amat(nc1)=aden(nc0+(j-1)*nsbas(im)+i)
     *               -aden(nc0+(i-1)*nsbas(im)+j)
!            smat(nc1)=aden(nc+i,nc+j)+aden(nc+j,nc+i)
!            amat(nc1)=aden(nc+i,nc+j)-aden(nc+j,nc+i)
            nc1=nc1+1
          enddo
          smat(nc1)=aden(nc0+(i-1)*nsbas(im)+i)
!          smat(nc1)=aden(nc+i,nc+i)
          amat(nc1)=0.d0
          nc1=nc1+1
        enddo
        nc=nc+nsbas(im)
        nc0=nc0+nsbas(im)**2
      enddo

!      open(100,file="dat1")
!      do i=1,1293
!        write(100,"(i8,1x,f14.8)") i,smat(i)
!      enddo
!      close(100)

      nsiz=nsi
! read property int and calculated property
      do i=1,npro
        call irdone(irtc,1,pname(i),ipcom(i),idummy,isymlb)
        if (irtc.eq.0) nsiz=idummy(1)
        call rdone(irtc,0,pname(i),ipcom(i),pint,isymlb)
C        print*, "nsiz",nsiz
        if(icall.eq.0) then
          pgauge(1,i)=pint(nsiz+1)
          pgauge(2,i)=pint(nsiz+2)
          pgauge(3,i)=pint(nsiz+3)
          pnuc(i)=pint(nsiz+4)
        endif
        if(isymlb.ne.1) then
          write(6,*) "error calcualte property,need debug"
          call abend
!          stop 888
        endif
!        print*, "pop int"
!        print"(10(1x,f8.5))",pint(1:100)
!        print*, "sden "
!        print"(10(1x,f8.5))",smat(1:100)

! calculate and save property
        sgn=1.d0
        if(pname(i)(1:5).eq."MLTPL") sgn=-1.d0
        if(ptyp(i)(1:4).eq."HERM") then
          val=sgn*ddot_(nsiz,smat,1,pint,1)
          vprop(istate,jstate,i)=val
          vprop(jstate,istate,i)=val
        else
          val=sgn*ddot_(nsiz,amat,1,pint,1)
          vprop(istate,jstate,i)=val
          vprop(jstate,istate,i)=-val
        endif
!        write(6,*) nsiz,"val=",val
      enddo

      icall=1

      return
      end

      subroutine transden(ngsm,nsbas,denao,cno,lmo,occ,loc)
! calculate density matrix in ao basis
      implicit none
      integer ngsm,lmo,loc
      integer :: nsbas(ngsm)
      real*8  :: denao(lmo),cno(lmo),occ(loc)
!      real*8  :: denao(loc,loc),cno(lmo),occ(loc)
!-----------------------------------------------
      integer i,im,ni,nc,nc0,nc1
      real*8 val

      denao=0.d0
!      print*, "occ",nsbas(1:ngsm)
!      print"(10(1x,f8.5))",occ
      nc=1
      nc0=0
      nc1=1
      do im=1,ngsm
        ni=nsbas(im)
        if(ni.eq.0) cycle
!        print*, nc,ni
        do i=1,ni
          val=occ(i+nc0)
          call dger(ni,ni,val,cno(nc1),1,cno(nc1),1,
     *              denao(nc),ni)
          nc1=nc1+ni
        enddo
        nc0=nc0+ni
        nc=nc+ni**2
      enddo

      return
      end
