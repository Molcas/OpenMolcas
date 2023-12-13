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
      Integer, parameter        :: wp=kind(0.d0)
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
      Integer, parameter        :: wp=kind(0.d0)
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
      Integer, parameter        :: wp=kind(0.d0)
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
      Integer, parameter        :: wp=kind(0.d0)
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
      Integer, parameter        :: wp=kind(0.d0)
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
      Integer, parameter        :: wp=kind(0.d0)
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

      Use Constants, only: Angstrom
      Implicit None
      Integer, parameter        :: wp=kind(0.d0)
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
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
      Integer                     :: njob, mxjob, mult, iss, ipar, ist
      Integer                     :: data_file_format
      Integer                     :: i,IsFreeUnit,Lu,Lutmp
      !Character(LEN=30)           :: fmt_int, fmt_real, fmt_key
      External                    :: IsFreeUnit
      Integer, allocatable        :: szproj(:), jbnum(:), mltplt(:)
      Integer, allocatable        :: nroot(:)
      Character(len=128)          :: Filename
      Character(len=1024)         :: molcas,fname,molcasversion
      LOGICAL                     :: dbg
      !
      Integer                     :: nAtoms, iAt, l
      Character(LEN=LENIN)        :: AtomLbl(MxAtom)
      Real*8, Allocatable         :: xyz(:,:)
      dbg=.false.

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
      Lutmp=IsFreeUnit(89)
      CALL molcas_open(Lutmp,fname)
      READ(Lutmp,'(A180)') molcasversion
      CLOSE(Lutmp)


      !-------------------------------------------------------------
      ! add coordinates
      nAtoms=0
      Call Get_iScalar('Unique atoms',nAtoms)
      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
      Call mma_allocate(xyz,3,8*nAtoms)
      Call Get_dArray('Unique Coordinates',xyz,3*nAtoms)

      !-------------------------------------------------------------
      ! write the data to the new aniso file
      FileName='ANISOFILE'
      Lu=IsFreeUnit(81)

      Call molcas_open(Lu,FileName)
      data_file_format=2021

      !fmt_key='(A)'
      !fmt_real='(5ES22.14,1x)'
      !fmt_int='(40(I0,1x))'

      WRITE(Lu,'(        A)') '# OPENMOLCAS interface to ANISO'
      !-------------------------------------------------------------
      ! ORIGIN of DATA file
      WRITE(Lu,'(        A)') '$source '
      WRITE(Lu,'(       2A)') 'MOLCAS  ',trim(molcas)
      WRITE(Lu,'(       2A)') 'VERSION ',trim(molcasversion)
      WRITE(Lu,'(        A)')
      !-------------------------------------------------------------
      ! DATA FILE FORMAT VERSION:
      WRITE(Lu,'(        A)') '$format '
      WRITE(Lu,'(40(I0,1x))')  data_file_format
      WRITE(Lu,'(        A)')
      !-------------------------------------------------------------
      ! Atom labels and coordinates
      WRITE(Lu,'(        A)') '$natoms '
      WRITE(Lu,'(40(I0,1x))')  nAtoms
      WRITE(Lu,'(        A)')
      WRITE(Lu,'(        A)') '$atomlbl'
      WRITE(Lu,'(40(I0,1x))')  nAtoms
      WRITE(Lu,'(40(A8,1x))')  (AtomLbl(iAt),iAt=1,nAtoms)
      WRITE(Lu,'(        A)')
      WRITE(Lu,'(        A)') '$coords (in Angstrom)'
      WRITE(Lu,'(40(I0,1x))')  nAtoms
      DO iAt=1,nAtoms
        WRITE(Lu,'(i3,1x,A8,1x,3(ES22.14,1x))')
     &                    iAt, AtomLbl(iAt), (Angstrom*XYZ(l,iAt),l=1,3)
      END DO
      WRITE(Lu,'(        A)')
      !-------------------------------------------------------------
      ! Number of spin orbit states
      CALL write_nss ( LU, nss, dbg )
      !-------------------------------------------------------------
      ! Number of spin free states
      CALL write_nstate ( LU, nstate, dbg )
      !-------------------------------------------------------------
      ! Number of spin multiplicities
      CALL write_nmult ( LU, njob, dbg )
      !-------------------------------------------------------------
      ! Values of the spin multiplicity
      CALL write_imult ( LU, njob, mltplt, dbg )
      !-------------------------------------------------------------
      ! Number of roots in each spin multiplicity
      CALL write_nroot ( LU, njob, nroot, dbg )
      !-------------------------------------------------------------
      ! SZ projection of the states defining the ording
      CALL write_szproj ( LU, nss, szproj, dbg )
      !-------------------------------------------------------------
      ! spin multiplicity for all states
      CALL write_multiplicity ( LU, nstate, multiplicity, dbg )
      !-------------------------------------------------------------
      ! Eigenvalues of the HBO matrix (spin-free energies from RASSCF)
      ! atomic units
      CALL write_eso ( LU, nss, eso_au, dbg )
      !-------------------------------------------------------------
      ! Eigenvalues of the HBO matrix (spin-free energies from RASSCF)
      ! atomic units
      CALL write_esfs ( LU, nstate, esfs_au, dbg )
      !-------------------------------------------------------------
      ! Angular momentum operator in the basis of spin free states,
      ! component X,Y,Z
      CALL write_angmom ( LU, nstate, angmom, dbg )
      !-------------------------------------------------------------
      ! AMFI/SOMF operator in the basis of spin free states
      ! (see eq. 35 in:
      ! D. Ganyushin and F. Neese, J. Chem. Phys. 138, 104113, 2013.
      ! Note that the meaning of AMFI integrals is different in
      ! MOLCAS and ORCA
      CALL write_amfi ( LU, nstate, amfi, dbg )
      !-------------------------------------------------------------
      ! Electric transition-dipole moments in the basis of spin
      ! free states
      ! The dipole moments, where bra and ket have the same state,
      ! are not computed. ??
      CALL write_edipmom ( LU, nstate, edmom, dbg )
      !-------------------------------------------------------------
      ! Magnetic dipole moment in the basis of SO states
      CALL write_magnetic_moment ( LU, nss, MM, dbg )
      !-------------------------------------------------------------
      ! Spin moment in the basis of SO states
      CALL write_spin_moment ( Lu, nss, MS, dbg )
      !-------------------------------------------------------------
      ! Electric transition-dipole moments in the basis of SO states
      CALL write_electric_moment ( LU, nss, DM, dbg )
      !-------------------------------------------------------------
      ! Eigenvectors of the HBO+SOC matrix
      CALL write_eigen ( LU, nss, U, dbg )
      !-------------------------------------------------------------
      ! The HBO+SOC matrix, atomic units
      CALL write_hso ( LU, nss, HSO, dbg )
      !-------------------------------------------------------------
      FLUSH(Lu)
      CLOSE(Lu)

      Call mma_deallocate(szproj)
      Call mma_deallocate(jbnum)
      Call mma_deallocate(mltplt)
      Call mma_deallocate(nroot)
      Call mma_deallocate(xyz)

      Return
      End






      Subroutine read_formatted_new_aniso( input_file_name, nss, nstate,
     &                                multiplicity, eso, esfs,
     &                                U, MM, MS, ML, DM, ANGMOM, EDMOM,
     &                                AMFI, HSO, eso_au, esfs_au )

      Implicit None
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
      Real(kind=8)  :: g_e,conv_au_to_cm1
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
      Call read_magnetic_moment(LuAniso,nss,MM,dbg)
      Call read_electric_moment(LuAniso,nss,DM,dbg)
      Call read_spin_moment(LuAniso,nss,MS,dbg)
      Call read_angmom(LuAniso,nstate,angmom,dbg)
      Call read_amfi(LuAniso,nstate,amfi,dbg)
      Call read_edipmom(LuAniso,nstate,EDMOM,dbg)
      Call read_multiplicity(LuAniso,nstate,multiplicity,dbg)
      Call read_eso(LuAniso,nss,eso_au,dbg)
      Call read_esfs(LuAniso,nstate,esfs_au,dbg)
      Call read_hso(LuAniso,nss,HSO,dbg)
      Call read_eigen(LuAniso,nss,U,dbg)

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


      Subroutine read_formatted_aniso_poly_NEW (
     &                                      input_file_name, nss,
     &                                      nstate, eso, MM, MS,
     &                                      iReturn )
      Implicit None
#include "stdalloc.fh"
      Integer, Parameter           :: wp=kind(0.d0)
      Integer, intent(inout)       :: nss, nstate, iReturn
      ! nLoc is the maximal value of the array nss(1:nneq)
      Real(kind=8), intent(out)    :: eso(nss)
      Complex(kind=8), intent(out) :: MM(3,nss,nss)
      Complex(kind=8), intent(out) :: MS(3,nss,nss)
      Character(len=180)           :: input_file_name
      ! local variables:
      Real(wp), allocatable        :: eso_au(:)
      Real(wp) :: conv_au_to_cm1
      Integer  :: i,j,l,LuAniso,IsFreeUnit
      External :: IsFreeUnit
      Logical  :: dbg

      dbg=.false.
      iReturn=0
      conv_au_to_cm1=2.194746313702E5_wp
      eso(1:nss)=0.0_wp
      MM(1:3,1:nss,1:nss)=(0.0_wp,0.0_wp)
      MS(1:3,1:nss,1:nss)=(0.0_wp,0.0_wp)
      Call mma_allocate(eso_au,nss,'eso_au')
      eso_au=0.0_wp

      LuAniso=IsFreeUnit(81)
      Call molcas_open(LuAniso,input_file_name)

      Call read_nss( LuAniso, nss, dbg )
      If(dbg) Write(6,*) 'read_formatted_aniso_poly_NEW: nss=',nss
      Call read_nstate( LuAniso, nstate, dbg )
      If(dbg) Write(6,*) 'read_formatted_aniso_poly_NEW: nstate=',nstate
      Call read_eso( LuAniso, nss, eso_au, dbg )
      If(dbg) Write(6,*) 'read_formatted_aniso_poly_NEW: eso_au=',
     &                   (eso_au(i),i=1,nss)
      Call read_magnetic_moment ( LuAniso, nss,
     &                            MM(1:3,1:nss,1:nss), dbg )
      IF (dbg) WRITE (6,*) 'Call read_spin_moment'
      FLUSH(6)
      Call read_spin_moment ( LuAniso, nss,
     &                        MS, dbg )

      ! compute the relative spin-orbit energies in cm-1
      do i=1,nss
        eso(i)=(eso_au(i)-eso_au(1))*conv_au_to_cm1
      end do
      Call mma_deallocate(eso_au)
      Close(LuAniso)

      IF (dbg) THEN
        WRITE(6,*) 'read_formatted_aniso_poly_NEW:  nss: ',nss
        WRITE(6,*) 'read_formatted_aniso_poly_NEW:   MM: '
        DO l=1,3
          WRITE(6,'(A,I0)') 'projection: L=',l
          DO i=1,nss
             WRITE(6,'(10(2F8.4,2x))')  ( MM(l,i,j), j=1,nss )
          END DO
        END DO

        WRITE(6,*) 'read_formatted_aniso_poly_NEW:   MS'
        DO l=1,3
          WRITE(6,'(A,I0)') 'projection: L=',l
          DO i=1,nss
             WRITE(6,'(10(2F8.4,2x))')  ( MS(l,i,j),j=1,nss )
          END DO
        END DO
      END IF


      Return
      End

