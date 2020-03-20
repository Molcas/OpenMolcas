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
      Real(kind=wp), intent(out)    :: eso(nss), esfs(nstate)
      Real(kind=wp), intent(out)    ::  edmom(3,nstate,nstate)
      Real(kind=wp), intent(out)    :: angmom(3,nstate,nstate)
      Real(kind=wp), intent(out)    ::   amfi(3,nstate,nstate)
      Complex(kind=wp), intent(out) :: MM(3,nss,nss)
      Complex(kind=wp), intent(out) :: MS(3,nss,nss)
      Complex(kind=wp), intent(out) :: ML(3,nss,nss)
!     electric dipole moment
      Complex(kind=wp), intent(out) :: DM(3,nss,nss)
      Complex(kind=wp), intent(out) ::   U(nss,nss)
      Complex(kind=wp), intent(out) :: HSO(nss,nss)
      Character(180)                :: input_file_name
      ! local variables:
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Real(kind=wp) :: g_e
      Real(kind=wp), allocatable :: tmpR(:,:), tmpI(:,:)
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
      Real(kind=wp), intent(out)    :: eso(nss)
      Complex(kind=wp), intent(out) :: MM(3,nss,nss)
      Complex(kind=wp), intent(out) :: MS(3,nss,nss)
      Complex(kind=wp), intent(out) :: ML(3,nss,nss)
      Character(180)                :: input_file_name
      ! local variables:
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Real(kind=wp) :: g_e
      Real(kind=wp), allocatable :: tmpR(:,:), tmpI(:,:)
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
      Real(kind=wp), intent(out)    :: eso(nLoc)
      Complex(kind=wp), intent(out) :: MM(3,nLoc,nLoc)
      Complex(kind=wp), intent(out) :: MS(3,nLoc,nLoc)
      Character(len=180)            :: input_file_name
      ! local variables:
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Integer       :: multiplicity(nLoc)
      Real(kind=wp), allocatable :: tmpR(:,:), tmpI(:,:)
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
      End









      Subroutine write_formatted_aniso( nss, nstate, multiplicity, eso,
     &                                  esfs, U, MM, MS, DM, angmom,
     &                                  edmom, amfi, HSO )

      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: nss, nstate, multiplicity(nstate)
      Real(kind=wp), intent(in)    :: eso(nss), esfs(nstate)
      Real(kind=wp), intent(in)    :: angmom(3,nstate,nstate)
      Real(kind=wp), intent(in)    ::  edmom(3,nstate,nstate)
      Real(kind=wp), intent(in)    ::   amfi(3,nstate,nstate)
      Complex(kind=wp), intent(in) :: MM(3,nss,nss)
      Complex(kind=wp), intent(in) :: MS(3,nss,nss)
      Complex(kind=wp), intent(in) :: DM(3,nss,nss)
      Complex(kind=wp), intent(in) ::    U(nss,nss)
      Complex(kind=wp), intent(in) ::  HSO(nss,nss)
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
      Real(kind=wp), intent(out)    :: eso(nss)
      Complex(kind=wp), intent(out) :: MM(3,nss,nss)
      Complex(kind=wp), intent(out) :: MS(3,nss,nss)
      Complex(kind=wp), intent(out) :: ML(3,nss,nss)
      Character(180)                :: input_file_name
      ! local variables:
      Integer       :: itmp, nss_local, nstate_local
      Integer       :: l,j,j1,j2,LuAniso,IsFreeUnit
      Real(kind=wp) :: g_e
      Real(kind=wp), allocatable :: tmpR(:,:), tmpI(:,:), tmp(:)
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
      read(LuAniso,*) (itmp,j=1,nstate_local)
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
      Real(kind=wp), intent(in)    :: eso(nss)
      Complex(kind=wp), intent(in) :: MM(3,nss,nss)
      Complex(kind=wp), intent(in) :: MS(3,nss,nss)
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

