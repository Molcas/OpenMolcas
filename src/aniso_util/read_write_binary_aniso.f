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
      Subroutine read_binary_aniso( nss, nstate, multiplicity, eso,
     &                              esfs, U, MM, MS, ML, DM, angmom,
     &                              edmom, amfi, HSO )
      Implicit None
      Integer, parameter :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer            :: nss, nstate
      Integer            :: multiplicity(nstate)
      Real(kind=8)      :: eso(nss), esfs(nstate)
      Real(kind=8)      :: angmom(3,nstate,nstate)
      Real(kind=8)      ::  edmom(3,nstate,nstate)
      Real(kind=8)      ::   amfi(3,nstate,nstate)
      Complex(kind=8)   :: MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss)
      Complex(kind=8)   :: DM(3,nss,nss)
      Complex(kind=8)   ::   U(nss,nss)
      Complex(kind=8)   :: HSO(nss,nss)
      ! local variables:
#include "stdalloc.fh"
      Integer            :: i, j, l
      Integer            :: luaniso, idisk, idum(1)
      Real(kind=8)      :: g_e
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)

      Call qEnter('read_binary_aniso')
      g_e=2.0023193043718_wp
      ! initialize:
      multiplicity=0
      eso=0.0_wp
      esfs=0.0_wp
      angmom=0.0_wp
       edmom=0.0_wp
        amfi=0.0_wp
      U =(0.0_wp,0.0_wp)
      MS=(0.0_wp,0.0_wp)
      ML=(0.0_wp,0.0_wp)
      MM=(0.0_wp,0.0_wp)
      DM=(0.0_wp,0.0_wp)
      HSO=(0.0_wp,0.0_wp)
      luaniso=8
      idisk=0
      ! get the information from binary "$project.aniso" file:
      Call daname(luaniso,'POLYFILE')
      Call idafile(luaniso,2,idum,1,idisk)
      nstate=idum(1)
      Call idafile(luaniso,2,idum,1,idisk)
      nss=idum(1)

      Call mma_allocate(tmpR,nss,nss,'tmpR')
      Call mma_allocate(tmpI,nss,nss,'tmpI')

      Call idafile(luaniso,2,multiplicity,nstate,idisk)
      Call ddafile(luaniso,2,eso,nss,idisk)
      Call ddafile(luaniso,2,esfs,nstate,idisk)

      ! spin-orbit mixing coefficients:
      tmpR=0.0_wp
      tmpI=0.0_wp
      Call ddafile(luaniso,2,tmpR,nss**2,idisk)
      Call ddafile(luaniso,2,tmpI,nss**2,idisk)
      Do i=1,nss
        Do j=1,nss
          U(i,j)=cmplx(tmpR(i,j),tmpI(i,j),wp)
        End Do
      End Do

      ! spin-orbit Hamiltonian
      tmpR=0.0_wp
      tmpI=0.0_wp
      Call ddafile(luaniso,2,tmpR,nss**2,idisk)
      Call ddafile(luaniso,2,tmpI,nss**2,idisk)
      Do i=1,nss
        Do j=1,nss
          HSO(i,j)=cmplx(tmpR(i,j),tmpI(i,j),wp)
        End Do
      End Do

      ! angmom
      Call ddafile(luaniso,3,angmom,3*nstate*nstate,idisk)
      ! edmom
      Call ddafile(luaniso,3, edmom,3*nstate*nstate,idisk)
      ! amfi
      Call ddafile(luaniso,3,  amfi,3*nstate*nstate,idisk)

      ! magnetic moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Call ddafile(luaniso,2,tmpR,nss**2,idisk)
        Call ddafile(luaniso,2,tmpI,nss**2,idisk)
        Do i=1,nss
          Do j=1,nss
            MM(l,i,j)=cmplx(tmpR(i,j),tmpI(i,j),wp)
          End Do
        End Do
      End Do

      ! spin moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Call ddafile(luaniso,2,tmpR,nss**2,idisk)
        Call ddafile(luaniso,2,tmpI,nss**2,idisk)
        Do i=1,nss
          Do j=1,nss
            MS(l,i,j)=cmplx(tmpR(i,j),tmpI(i,j),wp)
          End Do
        End Do
      End Do
      ! generate magnetic moment on the fly:
      Do l=1,3
        Do i=1,nss
          Do j=1,nss
            ML(l,i,j)=-MM(l,i,j)-g_e*MS(l,i,j)
          End Do
        End Do
      End Do

      ! electric dipole moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Call ddafile(luaniso,2,tmpR,nss**2,idisk)
        Call ddafile(luaniso,2,tmpI,nss**2,idisk)
        Do i=1,nss
          Do j=1,nss
            DM(l,i,j)=cmplx(tmpR(i,j),tmpI(i,j),wp)
          End Do
        End Do
      End Do

      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)
      Call daclos(luaniso)

      Call qExit('read_binary_aniso')
      Return
      End Subroutine read_binary_aniso



      Subroutine write_binary_aniso( nss, nstate, multiplicity, eso,
     &                               esfs, U, MM, MS, DM, ANGMOM,
     &                               EDMOM, AMFI, HSO )
      Implicit None
      Integer, parameter :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer            :: nss, nstate
      Integer            :: multiplicity(nstate)
      Real(kind=8)      :: eso(nss), esfs(nstate)
      Real(kind=8)      :: angmom(3,nstate,nstate)
      Real(kind=8)      ::  edmom(3,nstate,nstate)
      Real(kind=8)      ::   amfi(3,nstate,nstate)
      Complex(kind=8)   :: U(nss,nss), HSO(nss,nss)
      Complex(kind=8)   :: MM(3,nss,nss), MS(3,nss,nss), DM(3,nss,nss)
      ! local variables:
#include "stdalloc.fh"
      Integer            :: i, j, l
      Integer            :: luaniso, idisk
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)

      Call qEnter('write_binary_aniso')
      LUANISO=8
      ! open the binary $Project.aniso file
      Call DANAME(LUANISO,'POLYFILE')
      ! writing data to it:
      idisk=0
      Call idafile(luaniso,1,[nstate],1,idisk)
      Call idafile(luaniso,1,[nss],1,idisk)
      ! spin multiplicity ofthe spin free states:
      Call idafile(luaniso,1,multiplicity,nstate,idisk)
      ! spin-orbit energies:
      Call ddafile(luaniso,1,eso,nss,idisk)
      ! spin-free energies:
      Call ddafile(luaniso,1,esfs,nstate,idisk)

      Call mma_allocate(tmpR,nss,nss,'tmpR')
      Call mma_allocate(tmpI,nss,nss,'tmpI')
      ! spin-orbit mixing coefficients:
      tmpR=0.0_wp
      tmpI=0.0_wp
      Do i=1,nss
        Do j=1,nss
          tmpR(i,j)= dble( U(i,j) )
          tmpI(i,j)=aimag( U(i,j) )
        End Do
      End Do
      Call ddafile(luaniso,1,tmpR,nss*nss,idisk)
      Call ddafile(luaniso,1,tmpI,nss*nss,idisk)

      ! spin-orbit Hamiltonian
      tmpR=0.0_wp
      tmpI=0.0_wp
      Do i=1,nss
        Do j=1,nss
          tmpR(i,j)= dble( HSO(i,j) )
          tmpI(i,j)=aimag( HSO(i,j) )
        End Do
      End Do
      Call ddafile(luaniso,1,tmpR,nss*nss,idisk)
      Call ddafile(luaniso,1,tmpI,nss*nss,idisk)

      ! angmom
      Call ddafile(luaniso,1,angmom,3*nstate*nstate,idisk)
      ! edmom
      Call ddafile(luaniso,1, edmom,3*nstate*nstate,idisk)
      ! amfi
      Call ddafile(luaniso,1,  amfi,3*nstate*nstate,idisk)

      ! magnetic moment:
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do i=1,nss
          Do j=1,nss
            tmpR(i,j)= dble( MM(l,i,j) )
            tmpI(i,j)=aimag( MM(l,i,j) )
          End Do
        End Do
        Call ddafile(luaniso,1,tmpR,nss*nss,idisk)
        Call ddafile(luaniso,1,tmpI,nss*nss,idisk)
      End Do

      ! spin moment:
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do i=1,nss
          Do j=1,nss
            tmpR(i,j)= dble( MS(l,i,j) )
            tmpI(i,j)=aimag( MS(l,i,j) )
          End Do
        End Do
        Call ddafile(luaniso,1,tmpR,nss*nss,idisk)
        Call ddafile(luaniso,1,tmpI,nss*nss,idisk)
      End Do

      ! electric dipole moment:
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do i=1,nss
          Do j=1,nss
            tmpR(i,j)= dble( DM(l,i,j) )
            tmpI(i,j)=aimag( DM(l,i,j) )
          End Do
        End Do
        Call ddafile(luaniso,1,tmpR,nss*nss,idisk)
        Call ddafile(luaniso,1,tmpI,nss*nss,idisk)
      End Do

      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)
      ! close the binary $Project.aniso file
      Call daclos(luaniso)

      Call qExit('write_binary_aniso')
      Return
      End



