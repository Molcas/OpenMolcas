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
      subroutine frenkelexc(Frenkeltri, dim, nst1, nst2)

      use Constants, only: Zero, One, Two, Three, auToEV
      use Definitions, only: wp, iwp, u6
      use frenkel_global_vars, only: iTyp, jTyp, nestla, nestlb,
     &                               valst, excl
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL(kind=wp) (A-H,O-Z)
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "Morsel.fh"
#include "SysDef.fh"
#include "rassi.fh"
#include "jobin.fh"
#include "symmul.fh"

      integer(kind=iwp) :: d1lines, d2lines, dim, a, b,
     &                     nst1, nst2, ntrans, tri, iPL, run
      real(kind=wp), intent(inout) :: Frenkeltri(dim*(dim+1)/2)
      real(kind=wp), allocatable :: Frenkelsq(:), Frenkelquad(:,:),
     &                              Freninterm(:,:), Freninterm2(:,:)
      real(kind=wp), allocatable :: D1X(:), D1Y(:), D1Z(:), D2X(:),
     &               D2Z(:), Dip1(:,:,:), Dip2(:,:,:),
     &               EigEn(:),EDIFF(:), OSCSTR(:),Rassi1(:),
     &               Rassi2(:),Frenkeldia(:,:),D2Y(:),
     &               DipFrx(:,:), DipFry(:,:), DipFrz(:,:),
     &               DipFrintermx(:,:), DipFrintermy(:,:),
     &               DipFrintermz(:,:),DipFrintermx2(:,:),
     &               DipFrintermy2(:,:),DipFrintermz2(:,:),
     &               EigEnHa(:), E_FRENKEL(:)
      integer(kind=iwp), allocatable :: I1(:), F1(:), I2(:), F2(:)
      integer(kind=iwp), external :: isFreeUnit
      character(len=13) :: filnam, filnam1, filnam2, filnam3
      character(len=8), allocatable :: Hamelembra(:), Hamelemket(:)
#ifdef _DEBUGPRINT_RASSI_
      logical :: debug_rassi_code = .true.
#else
      logical :: debug_rassi_code = .false.
#endif


222   FORMAT (5X,2(1X,I4),5X,3(1X,ES18.8))
2     FORMAT(5X,(1x,I4,3x,I4,3x,F18.8,3x,F18.8,3x,
     &       F18.8,3x,F18.8,3x,F18.8))
3     FORMAT(6X,A,5x,A,6x,A,6x,A,19x,A,19x,A,19x,A)
443   FORMAT(3x,I4,3x,I4,3x,F18.8,3x,F18.8,3x,F18.8)
551   FORMAT(5X,A,4x,A,10x,A,20x,A,20x,A)

      call GETPRINTLEVEL
      iPL = iPrintLevel(-1)

      call mma_allocate(Frenkelsq,dim*dim)
      call mma_allocate(Frenkeldia,dim,dim)
      call mma_allocate(Freninterm,dim,dim)
      call mma_allocate(Freninterm2,dim,dim)
      call mma_allocate(Frenkelquad,dim,dim)
      call mma_allocate(Rassi1,nst1)
      call mma_allocate(Rassi2,nst2)
      call mma_allocate(E_FRENKEL,dim)

      ! read in rassi energies of both monomers to add them
      ! to the diagonal elements
      if (debug_rassi_code) then
        write(u6,*) 'iTyp=', iTyp, 'jTyp=', jTyp
      end if

      if (iTyp < jTyp) then
        write(u6,*) 'Monomer B was calculated first.'
        write(u6,*) 'Index 1 refers to A.'
      else
        write(u6,*) 'Monomer A was calculated first.'
        write(u6,*) 'Index 1 refers to B.'
        iTyp=2
        jTyp=1
      end if

      write(filnam2,'(A,I1)') 'stE', iTyp
      LuT3 = 11
      LuT1 = isFreeUnit(LuT3)

      call molcas_open(LuT1,filnam2)
      write(u6,'(A,1X,I1,A,I3.3)') 'Number of states of system'
     &      ,iTyp,'=',nst1
      write(u6,'(A,1X,I1,A,I3.3)') 'Number of states of system'
     &      ,jTyp,'=',nst2
      do i=1,nst1
        read(LuT1,*)Rassi1(i)
      end do
      close(LuT1)

      write(filnam3,'(A,I1)') 'stE', jTyp
      LuT3 = 11
      LuT1 = isFreeUnit(LuT3)

      call molcas_open(LuT1,filnam3)
      do i=1,nst2
        read(LuT1,*)Rassi2(i)
      end do
      close(LuT1)

      ! add rassi energies to diag of exciton coupl matrix
      do i=1,nst1
        write(u6,*)'State energy of system 1', Rassi1(i)
      end do
      do i=1,nst2
        write(u6,*)'State energy of system 2', Rassi2(i)
      end do

      tri=0
      do a=1,nst1
        do b=1,nst2
          if ((a/=1) .or. (b/=1)) then
            if (EXCL) then
              if ((ALL(nestla /= a)) .or. (ALL(nestlb /= b))) then
                cycle
              end if
            else
              if((a <= valst) .and. (b <= valst)) cycle
              if((a > valst) .and. (b > valst)) cycle
            end if
          end if
          tri = tri + 1
          Frenkeltri((tri*(tri+1))/2) = Frenkeltri((tri*(tri+1))/2)
     &          + Rassi1(a) + Rassi2(b)
          if (debug_rassi_code) then
            write(u6,*)'tri counter', tri
            write(u6,*) 'add energies 1 and 2', Rassi1(a), Rassi2(b)
          end if
        end do
      end do

      ! subtract gs energy
      gs = Frenkeltri(1)
      do i=1,dim
        Frenkeltri((i*(i+1))/2) = Frenkeltri((i*(i+1))/2) - gs
      end do
      ! convert to eV
      Frenkeltri(:) = Frenkeltri(:)*auToEV

      ! make a full Hamilton matrix from the lower triangular
      call square(Frenkeltri,Frenkelsq,1,dim,dim)
      if (debug_rassi_code) then
       write(u6,*) 'Frenkelsq:'
       do i=1,dim
         write(u6,'(1000ES18.8)') (Frenkelsq(j+(i-1)*dim),j=1,dim)
       enddo
      endif
      ! form a real quadratic matrix from Frenkelsq -> Frenkelquad
      run = 0
      do i=1,dim
        do j=1,dim
          run= run + 1
          Frenkelquad(i,j) = Frenkelsq(run)
        end do
      end do


      call mma_allocate(Hamelembra,dim)
      call mma_allocate(Hamelemket,dim)
      tri=0
      do a=1,nst1
        do b=1,nst2
          if ((a /= 1) .or. (b /= 1)) then
            if (EXCL) then
              if ((ALL(nestla /= a)) .or. (All(nestlb /= b))) then
                cycle
              end if
            else
              if ((a <= valst) .and. (b <= valst)) cycle
              if ((a > valst) .and. (b > valst)) cycle
            end if
          end if
          tri = tri + 1
          write(Hamelemket(tri),'(A,I3.3,I3.3,A)') '|',a,b,'>'
          write(Hamelembra(tri),'(A,I3.3,I3.3,A)') '<',a,b,'|'
        end do
      end do

      write(u6,*)'Frenkel Hamiltonian [eV]:'
      write(u6,'(11X)',advance='NO')
      write(u6,'(100(A,8X))',advance='YES') (Hamelemket(i),i=1,dim)

      do i=1,dim
        write(u6,'(A,1X,100(ES16.8))') Hamelembra(i),
     &        (Frenkelquad(i,j),j=1,dim)
      end do

      write(u6,*) 'Frenkel eigenvectors (column wise):'
      write(u6,'(11X)',advance="NO")
      write(u6,'(100(I3.3,13X))',advance="YES") (i,i=1,dim)

      call NIdiag_New(Frenkeltri, Frenkeldia, dim, dim)

      do i=1,dim
        write(u6,'(A,1X,100(ES16.8))') Hamelembra(i),
     &  (Frenkeldia(i,j),j=1,dim)
      end do

      write(u6,*) 'Hamiltonian eigenvalues [eV]:'
      do i=1,dim
        E_FRENKEL(i) = Frenkeltri((i*(i+1))/2)
        write(u6,*) Frenkeltri((i*(i+1))/2)
      end do

      call Add_Info('E_FRENKEL',E_FRENKEL,dim,6)

      call mma_deallocate(E_FRENKEL)
      call mma_deallocate(Hamelembra)
      call mma_deallocate(Hamelemket)


      if (debug_rassi_code) then
        call dgemm_('N','N',dim,dim,dim,One,Frenkelquad,dim,
     &              Frenkeldia,dim,Zero,Freninterm,dim)

        write(u6,*) 'first trafo H*U'
        do i=1,dim
          write(u6,'(1000ES18.8)') (Freninterm(i,j),j=1,dim)
        end do
        call dgemm_('T','N',dim,dim,dim,One,Frenkeldia,dim,
     &              Freninterm,dim,Zero,Freninterm2,dim)
        write(u6,*) 'second trafo U^(t)HU (should be diagonal)'
        do i=1,dim
          write(u6,'(1000ES18.8)') (Freninterm2(i,j),j=1,dim)
        end do
      end if

      ntrans = 0
      do i=1,dim-1
        ntrans = ntrans + i
      end do

      d1lines = 0
      d2lines = 0
      ! search for number of lines in the dipvec files
      write(filnam, '(A,I1)') 'dip_vec', iTyp
      LuT3 = 11
      LuT1 = isFreeUnit(LuT3)
      call molcas_open(LuT1, filnam)
      do
        read(LuT1, *, iostat=io)
        if (io/=0) exit
        d1lines = d1lines + 1
      end do
      close(LuT1)

      write(filnam1, '(A,I1)') 'dip_vec', jTyp
      LuT3 = 11
      LuT1 = isFreeUnit(LuT3)
      call molcas_open(LuT1, filnam1)
      do
        read(LuT1,*,iostat=io)
        if (io /= 0) exit
        d2lines = d2lines + 1
      end do
      close(LuT1)

      if (debug_rassi_code) then
        write(u6,*) 'd1lines ', d1lines
        write(u6,*) 'd2lines ', d2lines
        write(u6,*) 'ntrans' , ntrans
      end if

      call mma_allocate(I1,d1lines)
      call mma_allocate(F1,d1lines)
      call mma_allocate(D1X,d1lines)
      call mma_allocate(D1Y,d1lines)
      call mma_allocate(D1Z,d1lines)
      call mma_allocate(I2,d2lines)
      call mma_allocate(F2,d2lines)
      call mma_allocate(D2X,d2lines)
      call mma_allocate(D2Y,d2lines)
      call mma_allocate(D2Z,d2lines)
      call mma_allocate(Dip1,nst1,nst1,3)
      call mma_allocate(Dip2,nst2,nst2,3)
      call mma_allocate(EigEn,dim)
      call mma_allocate(EigEnHa,dim)
      call mma_allocate(EDIFF,ntrans)
      call mma_allocate(OSCSTR,ntrans)
      call mma_allocate(DipFrx,dim,dim)
      call mma_allocate(DipFry,dim,dim)
      call mma_allocate(DipFrz,dim,dim)
      call mma_allocate(DipFrintermx,dim,dim)
      call mma_allocate(DipFrintermx2,dim,dim)
      call mma_allocate(DipFrintermy,dim,dim)
      call mma_allocate(DipFrintermy2,dim,dim)
      call mma_allocate(DipFrintermz,dim,dim)
      call mma_allocate(DipFrintermz2,dim,dim)


      write(filnam, '(A,I1)') 'dip_vec', iTyp
      LuT3 = 11
      LuT1 = isFreeUnit(LuT3)
      call molcas_open(LuT1, filnam)
      do i=1,d1lines
        read(LuT1,222) I1(i), F1(i), D1X(i), D1Y(i), D1Z(i)
      end do
      close(LuT1)

      write(filnam1, '(A,I1)') 'dip_vec', jTyp
      LuT3 = 11
      LuT1 = isFreeUnit(LuT3)
      call molcas_open(LuT1, filnam1)
      do i=1,d2lines
        read(LuT1,222) I2(i), F2(i), D2X(i), D2Y(i), D2Z(i)
      end do
      close(LuT1)

      ! fill dipole matrices and use symmetry
      Dip1(:,:,:)=Zero
      do i=1,d1lines
        Dip1(I1(i),F1(i),1) = D1X(i)
        Dip1(I1(i),F1(i),2) = D1Y(i)
        Dip1(I1(i),F1(i),3) = D1Z(i)
        Dip1(F1(i),I1(i),1) = D1X(i)
        Dip1(F1(i),I1(i),2) = D1Y(i)
        Dip1(F1(i),I1(i),3) = D1Z(i)
      end do

      write(u6, *) 'Transition dipole vectors of system 1:'
      write(u6, 551) 'from','to','x','y','z'
      do i=1,d1lines
        write(u6,443)  I1(i), F1(i),
     &        (Dip1(I1(i),F1(i),j),j=1,3)
      end do

      Dip2(:,:,:) = Zero
      do i=1,d2lines
        Dip2(I2(i),F2(i),1) = D2X(i)
        Dip2(I2(i),F2(i),2) = D2Y(i)
        Dip2(I2(i),F2(i),3) = D2Z(i)
        Dip2(F2(i),I2(i),1) = D2X(i)
        Dip2(F2(i),I2(i),2) = D2Y(i)
        Dip2(F2(i),I2(i),3) = D2Z(i)
      end do

      write(u6, *) 'Transition dipole vectors of system 2:'
      write(u6, 551) 'from','to','x','y','z'
      do i=1,d2lines
        write(u6,443)  I2(i), F2(i),
     &        (Dip2(I2(i),F2(i),j),j=1,3)
      end do

      write(u6, *) 'relative excitonic state eigenenergies [eV]:'
      do i=1,dim
        EigEn(i) = Frenkeltri((i*(i+1))/2) - Frenkeltri(1)
        write(u6, *) 'EigEn:',i , EigEn(i)
      end do

      if (iPL >= 3) then
        write(u6,*) 'relative excitonic state '//
     &              'eigenenergies in Hartree:'
        do i=1,dim
          EigEnHa(i) = (Frenkeltri((i*(i+1))/2)-Frenkeltri(1))/auToEV
          write(u6,*) 'EigEn:',i , EigEnHa(i)
        end do
      end if

      ! create x
      a = 0
      b = 0
      DipFrx(:,:) = Zero
      do i=1,nst1
        do k=1,nst2
          if ((i /= 1) .or. (k /= 1)) then
            if (excl) then
              if ((All(nestla /= i)) .or. (ALL(nestlb /= k))) then
                cycle
              end if
            else
              if ((i <= valst) .and. (k <= valst)) cycle
              if ((i > valst) .and. (k > valst)) cycle
            end if
          end if
          a = a + 1
          b = 0
          do j=1,nst1
            do l=1,nst2
              if ((j /= 1) .or. (l /= 1)) then
                if (EXCL) then
                  if ((ALL(nestla /= j)) .or. (ALL(nestlb /= l))) then
                    cycle
                  end if
                else
                  if ((j <= valst) .and. (l <= valst)) cycle
                  if ((j > valst) .and. (l > valst)) cycle
                end if
              end if
              b = b + 1
              if (k == l) then
                DipFrx(a,b) = Dip1(i,j,1)
              end if
              if (i == j) then
                DipFrx(a,b) = Dip2(k,l,1)
              end if
            end do
          end do
        end do
      end do

      if (iPL >= 3) then
        write(u6,*) 'DipFrx:'
        do i=1,dim
          write(u6,'(100ES18.8)') (DipFrx(i,j),j=1,dim)
        end do
      end if
      ! transform dipole matrix in dipole basis
      call dgemm_('N','N',dim,dim,dim,One,DipFrx,dim,
     &            Frenkeldia,dim,Zero,DipFrintermx,dim)

      call dgemm_('T','N',dim,dim,dim,One,Frenkeldia,dim,
     &            DipFrintermx,dim,Zero,DipFrintermx2,dim)

      if (iPL >= 3) then
        write(u6,*) 'trafo U^(t)DU, dipole mtx in exc. basis, X'
        do i=1,dim
          write(u6,'(100ES18.8)') (DipFrintermx2(i,j),j=1,dim)
        end do
      end if
      ! create y
      a = 0
      b = 0
      DipFry(:,:) = Zero
      do i=1,nst1
        do k=1,nst2
          if ((i /= 1) .or. (k /= 1)) then
            if (EXCL) then
              if ((ALL(nestla /= i)) .or. (ALL(nestlb /= k))) then
                cycle
              end if
            else
              if ((i <= valst) .and. (k <= valst)) cycle
              if ((i > valst) .and. (k > valst)) cycle
            end if
          end if
          a = a + 1
          b = 0
          do j=1,nst1
            do l=1,nst2
              if ((j /= 1) .or. (l /= 1)) then
                if (EXCL) then
                  if ((ALL(nestla /= j)) .or. (ALL(nestlb /= l))) then
                    cycle
                  end if
                else
                  if ((j <= valst) .and. (l <= valst)) cycle
                  if ((j > valst) .and. (l > valst)) cycle
                end if
              end if
              b = b + 1
              if (k == l) then
                DipFry(a,b)=Dip1(I,J,2)
              end if
              if (i == j) then
                DipFry(a,b)=Dip2(K,L,2)
              end if
            end do
          end do
        end do
      end do

      if (iPL >= 3) then
        write(u6,*) 'DipFry:'
        do i=1,dim
          write(u6,'(100ES18.8)') (DipFry(i,j),j=1,dim)
        end do
      end if

      ! transform dipole matrix in dipole basis
      call dgemm_('N','N',dim,dim,dim,One,DipFry,dim,
     &            Frenkeldia,dim,Zero,DipFrintermy,dim)

      call dgemm_('T','N',dim,dim,dim,One,Frenkeldia,dim,
     &            DipFrintermy,dim,Zero,DipFrintermy2,dim)
      if (iPL >= 3) then
        write(u6,*) 'trafo U^(t)DU, dipole mtx in exc. basis, Y'
        do i=1,dim
          write(u6,'(100ES18.8)') (DipFrintermy2(i,j),j=1,dim)
        end do
      end if

      ! create z
      a = 0
      b = 0
      DipFrz(:,:) = Zero
      do i=1,nst1
        do k=1,nst2
          if ((i /= 1) .or. (k /= 1)) then
            if (excl) then
              if ((ALL(nestla /= i)) .or. (ALL(nestlb /= k))) then
                cycle
              end if
            else
              if ((i <= valst) .and. (k <= valst)) cycle
              if ((i > valst) .and. (k > valst)) cycle
            end if
          end if
          a = a + 1
          b = 0
          do j=1,nst1
            do l=1,nst2
              if ((j /= 1) .or. (l /= 1)) then
                if (excl) then
                  if ((ALL(nestla /= j)) .or. (ALL(nestlb /= l))) then
                    cycle
                  end if
                else
                  if ((j <= valst) .and. (l <= valst)) cycle
                  if ((j > valst) .and. (l > valst)) cycle
                end if
              end if
              b = b + 1
              if (k == l) then
                if ((k > 1) .and. (l > 1)) cycle
                DipFrz(a,b)=Dip1(I,J,3)
              end if
              if (i == j) then
                DipFrz(a,b) = Dip2(K,L,3)
              end if
            end do
          end do
        end do
      end do
      if (iPL >= 3) then
        write(u6,*) 'DipFrz:'
        do i=1,dim
          write(u6,'(100ES18.8)') (DipFrz(i,j),j=1,dim)
        end do
      end if

      ! transform dipole matrix in dipole basis
      call dgemm_('N','N',dim,dim,dim,One,DipFrz,dim,
     &            Frenkeldia,dim,Zero,DipFrintermz,dim)

      call dgemm_('T','N',dim,dim,dim,One,Frenkeldia,dim,
     &            DipFrintermz,dim,Zero,DipFrintermz2,dim)
      if (iPL >= 3) then
        write(u6,*) 'trafo U^(t)DU, dipole mtx in exc. basis, Z'
        do i=1,dim
          write(u6,'(100ES18.8)') (DipFrintermz2(i,j),j=1,dim)
        end do
      end if

      ! calc and print spectrum
      run = 0
      write(u6,*) 'Excitonic absorption spectrum'
      write(u6,3) 'from','to','excitation energy [eV]',
     &       'Dx','Dy','Dz','osc.str.'
      do i=1,dim
        do j=i+1,dim
          run = run + 1
          EDIFF(run) = EigEn(j) - EigEn(i)
          dipnorm = DipFrintermx2(i,j)*DipFrintermx2(i,j) +
     &          DipFrintermy2(i,j)*DipFrintermy2(i,j) +
     &          DipFrintermz2(i,j)*DipFrintermz2(i,j)
          OSCSTR(run) = (Two/Three)*EDIFF(run)*dipnorm
          OSCSTR(run) = OSCSTR(run)/auToEV
          write(u6,2) i,j, EDIFF(run), DipFrintermx2(i,j),
     &      DipFrintermy2(i,j), DipFrintermz2(i,j), OSCSTR(run)
        end do
      end do

      call Add_Info('FRENKEL_OSCSTR',OSCSTR,dim*(dim-1)/2,6)

      if (iPL >= 3) then
        run = 0
        write(u6,*) 'Excitonic absorption spectrum'
        write(u6,3) 'from','to','excitation energy [Ha]',
     &       'Dx','Dy','Dz','osc.str.'
        do i=1,dim
          do j=i+1,dim
            run = run + 1
            EDIFF(run) = (EigEn(j) - EigEn(i))/auToEV
            dipnorm = DipFrintermx2(i,j)*DipFrintermx2(i,j) +
     &          DipFrintermy2(i,j)*DipFrintermy2(i,j) +
     &          DipFrintermz2(i,j)*DipFrintermz2(i,j)
            OSCSTR(run) = (Two/Three)*EDIFF(run)*dipnorm
            write(u6,2) i,j, EDIFF(run), DipFrintermx2(i,j),
     &        DipFrintermy2(i,j), DipFrintermz2(i,j), OSCSTR(run)
          end do
        end do
      end if

      call mma_deallocate(I1)
      call mma_deallocate(F1)
      call mma_deallocate(D1X)
      call mma_deallocate(D1Y)
      call mma_deallocate(D1Z)
      call mma_deallocate(I2)
      call mma_deallocate(F2)
      call mma_deallocate(D2X)
      call mma_deallocate(D2Y)
      call mma_deallocate(D2Z)
      call mma_deallocate(Dip1)
      call mma_deallocate(Dip2)
      call mma_deallocate(EigEn)
      call mma_deallocate(EigEnHa)
      call mma_deallocate(EDIFF)
      call mma_deallocate(OSCSTR)
      call mma_deallocate(Frenkelsq)
      call mma_deallocate(Frenkeldia)
      call mma_deallocate(Freninterm)
      call mma_deallocate(Freninterm2)
      call mma_deallocate(Frenkelquad)
      call mma_deallocate(Rassi1)
      call mma_deallocate(Rassi2)
      call mma_deallocate(DipFrx)
      call mma_deallocate(DipFry)
      call mma_deallocate(DipFrz)
      call mma_deallocate(DipFrintermx)
      call mma_deallocate(DipFrintermx2)
      call mma_deallocate(DipFrintermy)
      call mma_deallocate(DipFrintermy2)
      call mma_deallocate(DipFrintermz)
      call mma_deallocate(DipFrintermz2)

      end subroutine frenkelexc
