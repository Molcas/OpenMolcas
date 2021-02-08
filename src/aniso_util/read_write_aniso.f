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
      Subroutine read_formatted_aniso( input_file_name, nss, nstate,
     &                                multiplicity, eso, esfs,
     &                                U, MM, MS, ML, DM, ANGMOM, EDMOM,
     &                                AMFI, HSO )
      Implicit None
#include "stdalloc.fh"
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(inout)        :: nss, nstate
      Integer, intent(out)          :: multiplicity(nstate)
      Real(kind=8), intent(out)    :: eso(nss), esfs(nstate)
      Real(kind=8), intent(out)    ::  edmom(3,nstate,nstate)
      Real(kind=8), intent(out)    :: angmom(3,nstate,nstate)
      Real(kind=8), intent(out)    ::   amfi(3,nstate,nstate)
      Complex(kind=8), intent(out) :: MM(3,nss,nss)
      Complex(kind=8), intent(out) :: MS(3,nss,nss)
      Complex(kind=8), intent(out) :: ML(3,nss,nss)
!     electric dipole moment
      Complex(kind=8), intent(out) :: DM(3,nss,nss)
      Complex(kind=8), intent(out) ::   U(nss,nss)
      Complex(kind=8), intent(out) :: HSO(nss,nss)
      Character(Len=180)            :: input_file_name
      ! local variables:
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Real(kind=8) :: g_e
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)
      External      :: IsFreeUnit
      logical :: DBG

      dbg=.false.

      If(dbg) write(6,'(A)') 'Enterring read_formatted_aniso'
      call xFlush(6)
      g_e=2.00231930437180_wp
c   set to zero all arrays:
      multiplicity=0
      eso=0.0_wp
      esfs=0.0_wp
      edmom=0.0_wp
      angmom=0.0_wp
      MM=(0.0_wp,0.0_wp)
      MS=(0.0_wp,0.0_wp)
      ML=(0.0_wp,0.0_wp)
      DM=(0.0_wp,0.0_wp)
       U=(0.0_wp,0.0_wp)
c  read the file "aniso.input":
      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,trim(input_file_name))
c compatibility with the present version: of aniso_i.input file
      read(LuAniso,*) nstate, nss
      If(dbg) write(6,'(A,2I6)') 'nstate, nss:',nstate, nss
      call xFlush(6)
      read(LuAniso,*) (eso(j),j=1,nss)
      If(dbg) then
         write(6,'(A)') 'ESO:'
         write(6,'(5ES24.14)') (eso(j),j=1,nss)
      End If
      read(LuAniso,*) (multiplicity(j),j=1,nstate)
      If(dbg) then
         write(6,'(A)') '(multiplicity(j),j=1,nstate)'
         write(6,'(50I3)') (multiplicity(j),j=1,nstate)
      End If
      call xFlush(6)

      Call mma_allocate(tmpR,nss,nss,'tmpR')
      Call mma_allocate(tmpI,nss,nss,'tmpI')
      ! magnetic moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MM(l,j1,j2) = cmplx( tmpR(j1,j2), tmpI(j1,j2), wp )
          End Do
        End Do
      End Do
      call xFlush(6)

      ! spin moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do

      ! spin-free energies
      read(LuAniso,*) (esfs(j),j=1,nstate)
      If(dbg) then
         write(6,'(A)') 'ESFS:'
         write(6,'(5ES24.14)') (esfs(j),j=1,nstate)
      End If

      ! U matrix
      tmpR=0.0_wp
      tmpI=0.0_wp
      Do j1=1,nss
        read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
      End Do

      Do j1=1,nss
        Do j2=1,nss
          U(j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
        End Do
      End Do

      ! angmom
      Do l=1,3
        Do j1=1,nstate
          Read(LuAniso,*) ( angmom(l,j1,j2), j2=1,nstate )
        End Do
      End Do

      ! compute the orbital moment
      Do l=1,3
        Do j1=1,nss
          Do j2=1,nss
            ML(l,j1,j2) = -MM(l,j1,j2) - MS(l,j1,j2)*g_e
          End Do
        End Do
      End Do

      ! DMmom
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            DM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do

      ! edmom
      Do l=1,3
        Do j1=1,nstate
          Read(LuAniso,*) (edmom(l,j1,j2),j2=1,nstate)
        End Do
      End Do

      ! amfi
      Do l=1,3
        Do j1=1,nstate
          Read(LuAniso,*) (amfi(l,j1,j2),j2=1,nstate)
        End Do
      End Do

      ! HSO matrix
      tmpR=0.0_wp
      tmpI=0.0_wp
      Do j1=1,nss
        read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
      End Do

      Do j1=1,nss
        Do j2=1,nss
          HSO(j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
        End Do
      End Do

      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)

      Close(LuAniso)
      Return
      End



      Subroutine read_formatted_aniso_old( input_file_name, nss, nstate,
     &                                multiplicity, eso, MM, MS, ML )
      Implicit None
#include "stdalloc.fh"
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(inout)        :: nss, nstate
      Integer, intent(out)          :: multiplicity(nstate)
      Real(kind=8), intent(out)    :: eso(nss)
      Complex(kind=8), intent(out) :: MM(3,nss,nss)
      Complex(kind=8), intent(out) :: MS(3,nss,nss)
      Complex(kind=8), intent(out) :: ML(3,nss,nss)
      Character(Len=180)            :: input_file_name
      ! local variables:
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Real(kind=8) :: g_e
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)
      External      :: IsFreeUnit

      g_e=2.00231930437180_wp
c   set to zero all arrays:
      multiplicity=0
      eso=0.0_wp
        MM=(0.0_wp,0.0_wp)
        MS=(0.0_wp,0.0_wp)
        ML=(0.0_wp,0.0_wp)
c  read the file "aniso.input":
      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,trim(input_file_name))
c compatibility with the present version: of aniso_i.input file
      read(LuAniso,*) nstate, nss
      read(LuAniso,*) (eso(j),j=1,nss)
      read(LuAniso,*) (multiplicity(j),j=1,nstate)

      Call mma_allocate(tmpR,nss,nss,'tmpR')
      Call mma_allocate(tmpI,nss,nss,'tmpI')
      ! magnetic moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do

      ! spin moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do
      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)

      ! compute the orbital moment
      Do l=1,3
        Do j1=1,nss
          Do j2=1,nss
            ML(l,j1,j2) = -MM(l,j1,j2) - MS(l,j1,j2)*g_e
          End Do
        End Do
      End Do

c      If (iprint.gt.4) Then
c        Write(6,'(10a12)') (('------------'),j=1,10)
c        Write(6,'(15x,a,i2,a)') 'aniso_.input'
c        Write(6,'(10a12)') (('------------'),j=1,10)
c        Write(6,'(5x,a,i6)') 'nstate = ',nstate
c        Write(6,'(5x,a,i6)') '   nss = ',nss
c        Write(6,'(5x,a)')    ' eso(j): '
c        Write(6,'(10(f12.5,1x))') (eso(j),j=1,nss)
c        Write(6,'(5x,a,i2,a)') 'multiplicity(j):'
c        Write(6,'(40i3)') (multiplicity(j),j=1,nstate)
c        Write(6,'(5x,a,i2,a)') 'dipso(l,j1,j2):'
c        Do l=1,3
c          Write(6,'(5x,a,i2)') 'axis= ',l
c          Do j1=1,nss
c            Write(6,'(20(2f20.14,1x))') (MM(l,j1,j2), j2=1,nss)
c          End Do
c        End Do
c        Write(6,*)
c        Write(6,'(5x,a,i2,a)') 's_so(l,j1,j2):'
c        Do l=1,3
c          Write(6,'(5x,a,i2)') 'axis= ',l
c          Do j1=1,nss
c            Write(6,'(20(2f20.14,1x))') (MS(l,j1,j2), j2=1,nss)
c          End Do
c        End Do
c      End If
      close(LuAniso)
      Return
      End


      Subroutine read_formatted_aniso_poly( input_file_name, nss,
     &                                      nstate, nLoc, eso, MM, MS,
     &                                      iReturn )
      Implicit None
#include "stdalloc.fh"
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(inout)        :: nss, nstate, nLoc, iReturn
      ! nLoc is the maximal value of the array nss(1:nneq)
      Real(kind=8), intent(out)    :: eso(nLoc)
      Complex(kind=8), intent(out) :: MM(3,nLoc,nLoc)
      Complex(kind=8), intent(out) :: MS(3,nLoc,nLoc)
      Character(len=180)            :: input_file_name
      ! local variables:
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Integer       :: multiplicity(nLoc)
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)
      External      :: IsFreeUnit

c   set to zero all arrays:
      iReturn=0
      multiplicity=0
      eso(:)=0.0_wp
      MM(:,:,:)=(0.0_wp,0.0_wp)
      MS(:,:,:)=(0.0_wp,0.0_wp)
c  read the file "aniso.input":
      LuAniso=IsFreeUnit(40)
      Call molcas_open( LuAniso, input_file_name )
c compatibility with the present version: of aniso_i.input file
      read(LuAniso,*) nstate, nss
      read(LuAniso,*) (eso(j),j=1,nss)
      read(LuAniso,*) (multiplicity(j),j=1,nstate)

      Call mma_allocate(tmpR,nss,nss,'tmpR')
      Call mma_allocate(tmpI,nss,nss,'tmpI')
      ! magnetic moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do

      ! spin moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do
      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)
c
c      If (iprint.gt.4) Then
c        Write(6,'(10a12)') (('------------'),j=1,10)
c        Write(6,'(15x,a,i2,a)') 'aniso_.input'
c        Write(6,'(10a12)') (('------------'),j=1,10)
c        Write(6,'(5x,a,i6)') 'nstate = ',nstate
c        Write(6,'(5x,a,i6)') '   nss = ',nss
c        Write(6,'(5x,a)')    ' eso(j): '
c        Write(6,'(10(f12.5,1x))') (eso(j),j=1,nss)
c        Write(6,'(5x,a,i2,a)') 'multiplicity(j):'
c        Write(6,'(40i3)') (multiplicity(j),j=1,nstate)
c        Write(6,'(5x,a,i2,a)') 'dipso(l,j1,j2):'
c        Do l=1,3
c          Write(6,'(5x,a,i2)') 'axis= ',l
c          Do j1=1,nss
c            Write(6,'(20(2f20.14,1x))') (MM(l,j1,j2), j2=1,nss)
c          End Do
c        End Do
c        Write(6,*)
c        Write(6,'(5x,a,i2,a)') 's_so(l,j1,j2):'
c        Do l=1,3
c          Write(6,'(5x,a,i2)') 'axis= ',l
c          Do j1=1,nss
c            Write(6,'(20(2f20.14,1x))') (MS(l,j1,j2), j2=1,nss)
c          End Do
c        End Do
c      End If
      Close(LuAniso)
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer_array(multiplicity)
#endif
      End









      Subroutine write_formatted_aniso( nss, nstate, multiplicity, eso,
     &                                  esfs, U, MM, MS, DM, angmom,
     &                                  edmom, amfi, HSO )

      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: nss, nstate, multiplicity(nstate)
      Real(kind=8), intent(in)    :: eso(nss), esfs(nstate)
      Real(kind=8), intent(in)    :: angmom(3,nstate,nstate)
      Real(kind=8), intent(in)    ::  edmom(3,nstate,nstate)
      Real(kind=8), intent(in)    ::   amfi(3,nstate,nstate)
      Complex(kind=8), intent(in) :: MM(3,nss,nss)
      Complex(kind=8), intent(in) :: MS(3,nss,nss)
      Complex(kind=8), intent(in) :: DM(3,nss,nss)
      Complex(kind=8), intent(in) ::    U(nss,nss)
      Complex(kind=8), intent(in) ::  HSO(nss,nss)
      ! local stuff
      Integer                      :: l,i,j,LuAniso,IsFreeUnit
      External                     :: IsFreeUnit

      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,'ANISOINPUT')
      Write(LuAniso,'(2i10)') nstate, nss
      Write(LuAniso,'(5ES24.14)') (eso(i),i=1,nss)
      Write(LuAniso,'(30i4)') (multiplicity(i),i=1,nstate)
      Do l=1,3
         Do i=1,nss
            Write(LuAniso,'(5ES24.14)') (MM(l,i,j),j=1,nss)
         End Do
      End Do
      Do l=1,3
         Do i=1,nss
            Write(LuAniso,'(5ES24.14)') (MS(l,i,j),j=1,nss)
         End Do
      End Do
      ! add data at the end, so that we do not break the functionality
      ! with the present format:
      Write(LuAniso,'(5ES24.14)') (esfs(i),i=1,nstate)
      Do i=1,nss
         Write(LuAniso,'(5ES24.14)') (U(i,j) ,j=1,nss)
      End Do
      ! angmom
      Do l=1,3
        Do i=1,nstate
          Write(LuAniso,'(5ES24.14)') (angmom(l,i,j),j=1,nstate)
        End Do
      End Do

      ! DMmom
      Do l=1,3
        Do i=1,nss
          Write(LuAniso,'(5ES24.14)') (DM(l,i,j),j=1,nss)
        End Do
      End Do

      ! edmom
      Do l=1,3
        Do i=1,nstate
          Write(LuAniso,'(5ES24.14)') (edmom(l,i,j),j=1,nstate)
        End Do
      End Do

      ! amfi
      Do l=1,3
        Do i=1,nstate
          Write(LuAniso,'(5ES24.14)') (amfi(l,i,j),j=1,nstate)
        End Do
      End Do

      ! HSO matrix
      Do i=1,nss
         Write(LuAniso,'(5ES24.14)') (HSO(i,j) ,j=1,nss)
      End Do

      Close(LuAniso)
      Return
      End




      Subroutine read_aniso_old_exch( input_file_name, nss, eso,
     &                                MM, MS, ML )
      Implicit None
#include "stdalloc.fh"
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)           :: nss
      Real(kind=8), intent(out)    :: eso(nss)
      Complex(kind=8), intent(out) :: MM(3,nss,nss)
      Complex(kind=8), intent(out) :: MS(3,nss,nss)
      Complex(kind=8), intent(out) :: ML(3,nss,nss)
      Character(Len=180)            :: input_file_name
      ! local variables:
      Integer       :: nss_local, nstate_local
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Real(kind=8) :: g_e
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:), tmp(:)
      External      :: IsFreeUnit
!     in this subroutine nss is input data




      g_e=2.00231930437180_wp
c   set to zero all arrays:
      eso=0.0_wp
        MM=(0.0_wp,0.0_wp)
        MS=(0.0_wp,0.0_wp)
        ML=(0.0_wp,0.0_wp)
      nss_local=0
      nstate_local=0

c  read the file "aniso.input":
      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,trim(input_file_name))
c compatibility with the present version: of aniso_i.input file
      read(LuAniso,*) nstate_local, nss_local
!---------------------------------------------------------
      Call mma_allocate(tmp,nss_local,'tmp')
      Call dcopy_(nss_local,[0.0_wp],0,tmp,1)
      ! local spin-orbit energy
      read(LuAniso,*) (tmp(j),j=1,nss_local)
      ! copy the lowest nss states to eso:
      Do j=1,nss
        eso(j) = tmp(j)
      End Do
      Call mma_deallocate(tmp)
!---------------------------------------------------------
      read(LuAniso,*) (l,j=1,nstate_local)
!---------------------------------------------------------

      Call mma_allocate(tmpR,nss_local,nss_local,'tmpR')
      Call mma_allocate(tmpI,nss_local,nss_local,'tmpI')
      ! magnetic moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss_local
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss_local )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do

      ! spin moment
      Do l=1,3
        tmpR=0.0_wp
        tmpI=0.0_wp
        Do j1=1,nss_local
          read(LuAniso,*) ( tmpR(j1,j2), tmpI(j1,j2), j2=1,nss_local )
        End Do
        Do j1=1,nss
          Do j2=1,nss
            MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
          End Do
        End Do
      End Do
      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)

      ! compute the orbital moment
      Do l=1,3
        Do j1=1,nss
          Do j2=1,nss
            ML(l,j1,j2) = -MM(l,j1,j2) - MS(l,j1,j2)*g_e
          End Do
        End Do
      End Do

      close(LuAniso)
      Return
      End




      Subroutine write_formatted_aniso_poly( filename, nss, eso, MM, MS)

      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: nss
      Real(kind=8), intent(in)    :: eso(nss)
      Complex(kind=8), intent(in) :: MM(3,nss,nss)
      Complex(kind=8), intent(in) :: MS(3,nss,nss)
      Character(len=180), intent(in):: filename
      ! local stuff
      Integer                      :: l,i,j,LuAniso,IsFreeUnit
      Integer                      :: nstate, multiplicity
      External                     :: IsFreeUnit

      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,filename)
      nstate=1
      multiplicity=1
      Write(LuAniso,'(2i10)') nstate, nss
      Write(LuAniso,'(5ES24.14)') (eso(i),i=1,nss)
      Write(LuAniso,'(30i4)') multiplicity
      Do l=1,3
         Do i=1,nss
            Write(LuAniso,'(5ES24.14)') (MM(l,i,j),j=1,nss)
         End Do
      End Do
      Do l=1,3
         Do i=1,nss
            Write(LuAniso,'(5ES24.14)') (MS(l,i,j),j=1,nss)
         End Do
      End Do
      Close(LuAniso)
      Return
      End









      Subroutine write_new_formatted_aniso(
     &                               nss, nstate, multiplicity, eso_au,
     &                               esfs_au, U, MM, MS, DM, angmom,
     &                               edmom, amfi, HSO )

      Implicit None
      Integer, Parameter          :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)         :: nss, nstate, multiplicity(nstate)
      Real(kind=8), intent(in)    :: eso_au(nss), esfs_au(nstate)
      Real(kind=8), intent(in)    :: angmom(3,nstate,nstate)
      Real(kind=8), intent(in)    ::  edmom(3,nstate,nstate)
      Real(kind=8), intent(in)    ::   amfi(3,nstate,nstate)
      Complex(kind=8), intent(in) :: MM(3,nss,nss)
      Complex(kind=8), intent(in) :: MS(3,nss,nss)
      Complex(kind=8), intent(in) :: DM(3,nss,nss)
      Complex(kind=8), intent(in) ::    U(nss,nss)
      Complex(kind=8), intent(in) ::  HSO(nss,nss)
      ! local stuff
#include "stdalloc.fh"
      Integer                     :: njob, mxjob, mult, iss, ipar, ist
      Integer                     :: data_file_format
      Integer                     :: i,IsFreeUnit,Lu,Lutmp
      Character(LEN=30)           :: fmt_int, fmt_real, fmt_key
      External                    :: IsFreeUnit
      Integer, allocatable        :: szproj(:), jbnum(:), mltplt(:)
      Integer, allocatable        :: nroot(:)
      Character(len=128)          :: Filename
      Character(len=1024)         :: molcas,fname,molcasversion


      !-------------------------------------------------------------
      ! some preparations
      njob=0
      mxjob=0
      Call get_iScalar('NJOB_SINGLE',njob)
      Call get_iScalar('MXJOB_SINGLE',mxjob)
      ! allocate temporary memory:
      Call mma_allocate(jbnum,nstate,'jbnum')
      Call mma_allocate(mltplt,mxjob,'mltplt')
      ! get the information from RUNFILE:
      mltplt=0
      jbnum=0
      Call get_iArray('MLTP_SINGLE',mltplt,mxjob)
      Call get_iArray('JBNUM_SINGLE',jbnum,nstate)
      !-------------------------------------------------------------


      !-------------------------------------------------------------
      ! prepare the szproj index table
      Call mma_allocate(szproj,nss,'szproj')
      szproj(1:nss)=0
      iss=0
      ipar=mod(multiplicity(1),2)
      Do Ist=1,nstate
         Mult=Multiplicity(Ist)
         Do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            If( (Ipar==0) .AND. (I==0)) Go To 310
               Iss=Iss+1
               szproj(iss)=I
  310       Continue
         End Do ! i
      End Do ! ist

      Call get_iScalar('MXJOB_SINGLE',mxjob)
      Call mma_allocate(nroot,mxjob,'nroot')
      nroot=0
      Call get_iArray('NSTAT_SINGLE',nroot,mxjob)


      !-------------------------------------------------------------
      ! Get the MOLCAS version: index table
      CALL getenvf('MOLCAS ',molcas)
      WRITE(fname,'(A)') trim(molcas)//'/.molcasversion'
      CALL molcas_open(Lutmp,fname)
      READ(Lutmp,'(A180)') molcasversion
      CLOSE(Lutmp)


      !-------------------------------------------------------------
      ! write the data to the new aniso file

      FileName='ANISOFILE'
      Lu=IsFreeUnit(81)

      Call molcas_open(Lu,FileName)

      data_file_format=2021

      fmt_key='(A)'
      fmt_real='(5ES22.14,1x)'
      fmt_int='(40(I0,1x))'

      WRITE(Lu,fmt_key) '# OPENMOLCAS interface to ANISO'
      !-------------------------------------------------------------
      ! ORIGIN of DATA file
      WRITE(Lu,fmt_key) '$source'
      WRITE(Lu,'(2A)')  'MOLCAS  ',trim(molcas)
      WRITE(Lu,'(2A)')  'VERSION ',trim(molcasversion)
      WRITE(Lu,'(A)')
      !-------------------------------------------------------------
      ! DATA FILE FORMAT VERSION:
      WRITE(Lu,fmt_key) '$format'
      WRITE(Lu,fmt_int)  data_file_format
      WRITE(Lu,'(A)')
      !-------------------------------------------------------------
      ! Number of spin orbit states
      WRITE(Lu,fmt_key) '$nss'
      WRITE(Lu,fmt_int)  nss
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Number of spin free states
      WRITE(Lu,fmt_key) '$nstate'
      WRITE(Lu,fmt_int)  nstate
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Number of spin multiplicities
      WRITE(Lu,fmt_key)  '$nmult'
      WRITE(Lu,fmt_int)  njob
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Values of the spin multiplicity
      WRITE(Lu,fmt_key)  '$imult'
      WRITE(Lu,fmt_int)  njob
      WRITE(Lu,fmt_int)  mltplt(1:njob)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Number of roots in each spin multiplicity
      WRITE(Lu,fmt_key)  '$nroot'
      WRITE(Lu,fmt_int)  njob
      WRITE(Lu,fmt_int)  nroot(1:njob)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! SZ projection of the states defining the ording
      WRITE(Lu,fmt_key)  '$szproj'
      WRITE(Lu,fmt_int)  nss
      WRITE(Lu,fmt_int)  szproj(1:nss)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! spin multiplicity for all states
      WRITE(Lu,fmt_key)  '$multiplicity'
      WRITE(Lu,fmt_int)  nstate
      WRITE(Lu,fmt_int)  multiplicity(1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Eigenvalues of the HBO matrix (spin-free energies from RASSCF)
      ! atomic units
      WRITE(Lu,fmt_key)  '$eso'
      WRITE(Lu,fmt_int)  nss
      WRITE(Lu,fmt_real) (eso_au(i),i=1,nss)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Eigenvalues of the HBO matrix (spin-free energies from RASSCF)
      ! atomic units
      WRITE(Lu,fmt_key)  '$esfs'
      WRITE(Lu,fmt_int)  nstate
      WRITE(Lu,fmt_real) (esfs_au(i),i=1,nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Angular momentum operator in the basis of spin free states,
      ! component X,Y,Z
      WRITE(Lu,fmt_key)  '$angmom_xi'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) angmom(1,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$angmom_yi'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) angmom(2,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$angmom_zi'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) angmom(3,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! AMFI/SOMF operator in the basis of spin free states
      ! (see eq. 35 in [1]).
      ! Note that the meaning of these matrices is different in
      ! MOLCAS and ORCA
      WRITE(Lu,fmt_key)  '$amfi_x'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) amfi(1,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$amfi_y'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) amfi(2,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$amfi_z'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) amfi(3,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Electric transition-dipole moments in the basis of spin
      ! free states
      ! The dipole moments, where bra and ket have the same state,
      ! are not computed. ??
      WRITE(Lu,fmt_key)  '$edmom_x'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) edmom(1,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edmom_y'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) edmom(2,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edmom_z'
      WRITE(Lu,fmt_int)  nstate, nstate
      WRITE(Lu,fmt_real) edmom(3,1:nstate,1:nstate)
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Magnetic dipole moment in the basis of SO states
      WRITE(Lu,fmt_key)  '$magn_xr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real)  DBLE(MM(1,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$magn_xi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(MM(1,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$magn_yr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real)  DBLE(MM(2,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$magn_yi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(MM(2,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$magn_zr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real)  DBLE(MM(3,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$magn_zi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(MM(3,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Spin moment in the basis of SO states
      WRITE(Lu,fmt_key)  '$spin_xr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(MS(1,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$spin_xi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(MS(1,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$spin_yr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(MS(2,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$spin_yi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(MS(2,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$spin_zr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(MS(3,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$spin_zi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(MS(3,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Electric transition-dipole moments in the basis of SO states
      WRITE(Lu,fmt_key)  '$edipm_xr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(DM(1,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edipm_xi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(DM(1,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edipm_yr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(DM(2,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edipm_yi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(DM(2,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edipm_zr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(DM(3,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      WRITE(Lu,fmt_key)  '$edipm_zi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(DM(3,1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! Eigenvectors of the HBO+SOC matrix
      WRITE(Lu,fmt_key)  '$eigenr'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(U(1:nss,1:nss))
      WRITE(Lu,'(A)')
      WRITE(Lu,fmt_key)  '$eigeni'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(U(1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      ! The HBO+SOC matrix, atomic units
      WRITE(Lu,fmt_key)  '$hsor'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DBLE(HSO(1:nss,1:nss))
      WRITE(Lu,'(A)')
      WRITE(Lu,fmt_key)  '$hsoi'
      WRITE(Lu,fmt_int)  nss, nss
      WRITE(Lu,fmt_real) DIMAG(HSO(1:nss,1:nss))
      WRITE(Lu,'(A)')
      FLUSH(Lu)
      !-------------------------------------------------------------
      FLUSH(Lu)
      CLOSE(Lu)

      Call mma_deallocate(szproj)
      Call mma_deallocate(jbnum)
      Call mma_deallocate(mltplt)
      Call mma_deallocate(nroot)

      Return
      End






      Subroutine read_formatted_new_aniso( input_file_name, nss, nstate,
     &                                multiplicity, eso, esfs,
     &                                U, MM, MS, ML, DM, ANGMOM, EDMOM,
     &                                AMFI, HSO, eso_au, esfs_au )
      USE io_data

      Implicit None
#include "stdalloc.fh"
      Integer, Parameter           :: wp=kind(0.d0)
      Integer, intent(inout)       :: nss, nstate
      Integer, intent(out)         :: multiplicity(nstate)
      Real(kind=8), intent(out)    :: eso(nss), esfs(nstate)
      Real(kind=8), intent(out)    :: eso_au(nss), esfs_au(nstate)
      Real(kind=8), intent(out)    ::  edmom(3,nstate,nstate)
      Real(kind=8), intent(out)    :: angmom(3,nstate,nstate)
      Real(kind=8), intent(out)    ::   amfi(3,nstate,nstate)
      Complex(kind=8), intent(out) :: MM(3,nss,nss)
      Complex(kind=8), intent(out) :: MS(3,nss,nss)
      Complex(kind=8), intent(out) :: ML(3,nss,nss)
!     electric dipole moment
      Complex(kind=8), intent(out) :: DM(3,nss,nss)
      Complex(kind=8), intent(out) ::   U(nss,nss)
      Complex(kind=8), intent(out) :: HSO(nss,nss)
      Character(Len=180)           :: input_file_name
      ! local variables:
      Integer       :: l,i,j,LuAniso,IsFreeUnit
      Real(kind=8)  :: g_e,conv_au_to_cm1, gtens(3), maxes(3,3)
      External      :: IsFreeUnit
      logical       :: DBG

      dbg=.false.
      If(dbg) write(6,'(A)') 'Enterring read_formatted_aniso_new'
      If(dbg) call xFlush(6)
      g_e=2.00231930437180_wp
      conv_au_to_cm1=2.194746313702E5_wp
!     set to zero all arrays:
      multiplicity=0
      eso=0.0_wp
      esfs=0.0_wp
      edmom=0.0_wp
      angmom=0.0_wp
      MM=(0.0_wp,0.0_wp)
      MS=(0.0_wp,0.0_wp)
      ML=(0.0_wp,0.0_wp)
      DM=(0.0_wp,0.0_wp)
       U=(0.0_wp,0.0_wp)
!     read the data file:
      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,input_file_name)
      Call read_magnetic_moment(LuAniso,nss,MM)
      Call read_electric_moment(LuAniso,nss,DM)
      Call read_spin_moment(LuAniso,nss,MS)
      Call read_angmom(LuAniso,nstate,angmom)
      Call read_amfi(LuAniso,nstate,amfi)
      Call read_edipmom(LuAniso,nstate,EDMOM)
      Call read_multiplicity(LuAniso,nstate,multiplicity)
      Call read_eso(LuAniso,nss,eso_au)
      Call read_esfs(LuAniso,nstate,esfs_au)
      Call read_hso(LuAniso,nss,HSO)
      Call read_eigen(LuAniso,nss,U)

      ! compute the relative spin-orbit energies in cm-1
      do i=1,nss
        eso(i)=(eso_au(i)-eso_au(1))*conv_au_to_cm1
      end do

      ! compute the relative spin-free energies in cm-1
      do i=1,nstate
        esfs(i)=(esfs_au(i)-esfs_au(1))*conv_au_to_cm1
      end do
      write(6,*) esfs

      ! compute the orbital moment
      do l=1,3
        do i=1,nss
          do j=1,nss
            ML(l,i,j) = -MM(l,i,j) - MS(l,i,j)*g_e
          end do
        end do
      end do

      !read_nmult,
      !read_imult,
      !read_format,
      !read_nroot,
      !read_szproj,

      Close(LuAniso)
      Return
      End


