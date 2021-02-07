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
!
! Library of routines for Dynamix module
!
!***********************************************************************
!      SUBROUTINE DxRdNAtomStnd(natom)
!      SUBROUTINE DxRdNAtomHbrd(natom)
!      SUBROUTINE DxRdStnd(natom,atom,xyz,force)
!      SUBROUTINE DxRdHbrd(natom,atom,xyz,force)
!      SUBROUTINE DxPtTableCo(caption,time,natom,atom,xyz,lastline,M,fo)
!      SUBROUTINE DxPtTableWithoutMassForce(caption,time,natom,atom,xyz)
!      SUBROUTINE DxRdOut(pcoo,POUT,natom)
!      SUBROUTINE DxRdVel(vel,natom)
!      SUBROUTINE DxWtVel(vel,natom3)
!      SUBROUTINE DxCoord(natom,atom,xyz,hybrid)
!      SUBROUTINE DxEnergies(time,Epot,Ekin,Etot)
!      SUBROUTINE Put_Velocity(vel,natom3)
!      SUBROUTINE Get_Velocity(vel,natom3)
!      SUBROUTINE Get_NHC(NHC,nh)
!      SUBROUTINE Put_NHC(NHC,nh)
!***********************************************************************
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
!
!     Read the number of atoms. This needs to be done separately, because
!     dynamic-size matrices depend on the result.
!
      SUBROUTINE DxRdNAtomStnd(natom)
      IMPLICIT NONE
      INTEGER natom
      CALL Get_nAtoms_Full(natom)
      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
!
!     Read the number of atoms. This needs to be done separately, because
!     dynamic-size matrices depend on the result.
!
      SUBROUTINE DxRdNAtomHbrd(natom)
      IMPLICIT NONE
      EXTERNAL IsFreeUnit
      INTEGER natom,file,IsFreeUnit
      CHARACTER filname*80
      file=81
      file=IsFreeUnit(file)
      filname='fixforce.dmx'
      Call Molcas_Open(file,filname)
      READ(file,100) natom
      CLOSE(file)
 100  FORMAT(I6)
      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxRdStnd(natom,atom,xyz,force)
!
!     Read in the atom names, coordinates and forces from RUNFILE for
!     the standard QM molecular dynamics.
!
!     IS 14/06-2007
!
      IMPLICIT NONE
#include "Molcas.fh"
#include "WrkSpc.fh"
      INTEGER       natom
      REAL*8        xyz(natom*3),force(natom*3),conv
      CHARACTER     atom(natom)*2
      PARAMETER     (conv=-1.0d0)
!
!     The parameter conv converts the gradients (Hartree/Bohr)
!                                  to forces (Hartree/Bohr)
!
      CALL Get_Name_Full(atom)
      CALL Get_Coord_Full(xyz,natom)
!
!     Read the gradients from RUNFILE and convert them to forces.
!
      CALL Get_Grad_Full(force,natom)
      call dscal_(3*natom,conv,force,1)
!
      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxRdHbrd(natom,atom,xyz,force)
!
!     Read in the atom names, coordinates and forces from files
!     fixforce.dmx and prmcrd2.
!
      IMPLICIT NONE
#include "Molcas.fh"
#include "constants2.fh"
      External IsFreeUnit
      INTEGER       i,j,natom,natom2,file,IsFreeUnit
      REAL*8        xyz(natom*3),force(natom*3),conv,a2bohr
      CHARACTER     filname*80,atom(natom)*2
      PARAMETER     (a2bohr=1.0d0/Angstrom)
      PARAMETER     (conv=Angstrom/CONV_AU_TO_KJ_PER_MOLE_)
!
!     The parameter conv converts the forces from Hartree/Bohr
!                        to kJ/mole/Agstrom   => 1/4961.475514610d0
!                 a2bohr converts the coordinates from Angstrom
!                        to Bohr              => 1/0.52917720859
!
      file=81
      file=IsFreeUnit(file)
      filname='fixforce.dmx'
      Call Molcas_Open(file,filname)
      READ(file,*)
      DO i=1, natom
         READ(file,101)(force((i-1)*3+j),j=1,3),atom(i)
      END DO
      CLOSE(file)
!
!     Convert forces from kJ/mole/Angstrom to Hartree/Bohr
!
      call dscal_(3*natom,conv,force,1)
!
!     Read coordinates from fiel 'prmcrd2' and check for consistency
!     with forces.
!
      file=IsFreeUnit(file)
      filname='prmcrd2'
      Call Molcas_Open(file,filname)
      READ(file,'(/,I6)') natom2
      IF (natom2.ne.natom) STOP "Inconsistency between coordinates"
      READ(file,102)(xyz(i),i=1,natom*3)
      CLOSE(file)
!
!     Convert coordinates from Angstrom to Bohr
!
      call dscal_(3*natom,a2bohr,xyz,1)

 101  FORMAT(3E21.14,A)
 102  FORMAT(6F12.7)

      RETURN

      CALL Abend()

      END

!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxPtTableCo(caption,time,natom,atom,xyz,lastline,M,fo)
!
!     Prints a nice table in the output file
!
      IMPLICIT NONE
#include "Molcas.fh"
      INTEGER     i,j,natom,k
      CHARACTER   caption*15, lastline*80, atom(natom)*2
      REAL*8      time,xyz(natom*3),M(natom),fo(natom*3)

      DO i=1, 3
        WRITE(6,*)
      END DO

      WRITE(6,100) caption,' (time = ',time,' a.u.):'
      WRITE(6,102)'-------------------------------------------------'// &
     &            '---------------------------------------------'
      WRITE(6,102)'      No. Atom    X          Y          Z        '// &
     &            'Mass       F(x)         F(y)         F(z)'
      WRITE(6,102)'-------------------------------------------------'// &
     &            '---------------------------------------------'

      DO i=1, natom
        WRITE(6,101) i,atom(i),(xyz(3*(i-1)+j),j=1,3)                   &
     &              ,M(i),(fo(3*(i-1)+k),k=1,3)
      END DO

      WRITE(6,102)'-------------------------------------------------'// &
     &            '---------------------------------------------'
      WRITE(6,102) trim(lastline)

      DO i=1, 3
       WRITE(6,*)
      END DO

 100  FORMAT(A22,A7,F8.1,A)
 101  FORMAT(6X,I4,A3,3(1X,F10.6),1X,ES9.2,3(1X,ES12.5))
 102  FORMAT(1X,A)

      RETURN

      END

!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxPtTableWithoutMassForce(caption,time,natom,atom,xyz)
!
!     Prints a nice table in the output file for the initial velcoities
!
      IMPLICIT NONE
#include "Molcas.fh"
      INTEGER     i,j,natom
      CHARACTER   caption*15, atom(natom)*2
      REAL*8      time,xyz(natom*3)

      DO i=1, 3
        WRITE(6,*)' '
      END DO

      WRITE(6,102) caption,' (time = ',time,' a.u.):'
      WRITE(6,*)'----------------------------------------------'
      WRITE(6,*)'     No.  Atom    X          Y          Z     '
      WRITE(6,*)'----------------------------------------------'

      DO i=1, natom
        WRITE(6,103) '      ',i,atom(i),(xyz(3*(i-1)+j),j=1,3)
      END DO

      WRITE(6,*)'----------------------------------------------'

      DO i=1, 3
       WRITE(6,*)' '
      END DO

 102  FORMAT(A22,A7,F8.1,A)
 103  FORMAT(A6,I4,A3,3F11.6)

      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxRdOut(pcoo,POUT,natom)
!
!     This Subroutine reads in the coordinates to project out from the file 'out.00N.xyz'
!
      IMPLICIT NONE
#include "stdalloc.fh"
#include "real.fh"
      EXTERNAL      IsFreeUnit
      INTEGER       i,j,p,file,natom,POUT,IsFreeUnit,mn
      REAL*8        pcoo(POUT,natom*3)
      CHARACTER*80  filname
      CHARACTER*180 OutLine, Get_Ln
      REAL*8, ALLOCATABLE :: UMat(:,:),VMat(:,:),S(:),SqrtM(:)
!
      file=81
      file=IsFreeUnit(file)
      DO P=1,POUT
        WRITE(filname,'(A,I3.3,A)') 'out.',p,'.xyz'
        CALL Molcas_Open(file,filname)
        DO i=1,natom
          OutLine = Get_Ln(file)
          CALL UpCase(OutLine)
          DO j=1,3
            pcoo(p,3*(i-1)+j)=0.0D0
            CALL Get_F(j,pcoo(p,3*(i-1)+j),1)
          END DO
        END DO
      END DO
      CLOSE(file)
!
! Orthonormalize the vectors in mass-weighted coordinates
!
      CALL mma_Allocate(SqrtM,natom)
      CALL GetMassDx(SqrtM,natom)
      DO i=1, natom
        SqrtM(i)=Sqrt(SqrtM(i))
        j=3*(i-1)+1
        pcoo(:,j:j+2)=pcoo(:,j:j+2)*SqrtM(i)
      END DO
      mn=MIN(POUT,3*natom)
      CALL mma_Allocate(UMat,POUT,mn)
      CALL mma_Allocate(VMat,mn,3*natom)
      CALL mma_Allocate(S,mn)
      CALL large_svd(POUT,3*natom,pcoo,UMat,VMat,S)
      CALL dgemm_('N','N',POUT,3*natom,mn,One,UMat,POUT,VMat,mn,        &
     &                                    Zero,pcoo,POUT)
      DO i=1, natom
        j=3*(i-1)+1
        pcoo(:,j:j+2)=pcoo(:,j:j+2)/SqrtM(i)
      END DO
      CALL mma_deAllocate(UMat)
      CALL mma_deAllocate(VMat)
      CALL mma_deAllocate(S)
      CALL mma_deAllocate(SqrtM)
!
      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxRdIn(pcoo,PIN,natom)
!
!     This Subroutine reads in the coordinates to keep in from the file 'in.00N.xyz'
!
      IMPLICIT NONE
#include "stdalloc.fh"
#include "real.fh"
      EXTERNAL      IsFreeUnit
      INTEGER       i,j,p,file,natom,PIN,IsFreeUnit,mn
      REAL*8        pcoo(PIN,natom*3)
      CHARACTER*80  filname
      CHARACTER*180 OutLine, Get_Ln
      REAL*8, ALLOCATABLE :: UMat(:,:),VMat(:,:),S(:),SqrtM(:)
!
      file=81
      file=IsFreeUnit(file)
      DO P=1,PIN
        WRITE(filname,'(A,I3.3,A)') 'in.',p,'.xyz'
        CALL Molcas_Open(file,filname)
        DO i=1,natom
          OutLine = Get_Ln(file)
          CALL UpCase(OutLine)
          DO j=1,3
            pcoo(p,3*(i-1)+j)=0.0D0
            CALL Get_F(j,pcoo(p,3*(i-1)+j),1)
          END DO
        END DO
      END DO
      CLOSE(file)
!
! Orthonormalize the vectors in mass-weighted coordinates
!
      CALL mma_Allocate(SqrtM,natom)
      CALL GetMassDx(SqrtM,natom)
      DO i=1, natom
        SqrtM(i)=Sqrt(SqrtM(i))
        j=3*(i-1)+1
        pcoo(:,j:j+2)=pcoo(:,j:j+2)*SqrtM(i)
      END DO
      mn=MIN(PIN,3*natom)
      CALL mma_Allocate(UMat,PIN,mn)
      CALL mma_Allocate(VMat,mn,3*natom)
      CALL mma_Allocate(S,mn)
      CALL large_svd(PIN,3*natom,pcoo,UMat,VMat,S)
      CALL dgemm_('N','N',PIN,3*natom,mn,One,UMat,PIN,VMat,mn,          &
     &                                    Zero,pcoo,PIN)
      DO i=1, natom
        j=3*(i-1)+1
        pcoo(:,j:j+2)=pcoo(:,j:j+2)/SqrtM(i)
      END DO
      CALL mma_deAllocate(UMat)
      CALL mma_deAllocate(VMat)
      CALL mma_deAllocate(S)
      CALL mma_deAllocate(SqrtM)
!
      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxRdVel(vel,natom)
!
!     This Subroutine reads in the velocities from the file 'velocity.xyz'
!
      IMPLICIT NONE
#include "Molcas.fh"
      EXTERNAL      IsFreeUnit
      INTEGER       i,j,file,natom,IsFreeUnit
      REAL*8        vel(natom*3)
      CHARACTER*80  filname
      CHARACTER*180 VelLine, Get_Ln
!
      file=81
      file=IsFreeUnit(file)
      filname='velocity.xyz'
      CALL Molcas_Open(file,filname)
      DO i=1,natom
        VelLine = Get_Ln(file)
        CALL UpCase(VelLine)
        DO j=1,3
          vel(3*(i-1)+j)=0.0D0
          CALL Get_F(j,vel(3*(i-1)+j),1)
        END DO
      END DO
!      READ(file,100)(vel(i),i=1,3*natom3)
      CLOSE(file)
!
!100  FORMAT(3D18.10)

      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxWtVel(vel,natom3)
!
!     This Subroutine writes the velocities to the file 'velocity.xyz'
!
      IMPLICIT NONE
#include "Molcas.fh"
      EXTERNAL    IsFreeUnit
      INTEGER     i,file,natom3,IsFreeUnit
      REAL*8      vel(natom3)
      CHARACTER*80  filname

      file=81
      file=IsFreeUnit(file)
      filname='velocity.xyz'
      CALL Molcas_Open(file,filname)
      WRITE(file,100)(vel(i),i=1,natom3)
      CLOSE(file)

 100  FORMAT(3D18.10)

      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxCoord(natom,atom,xyz,hybrid)
!
!     Write the coordinates to a file in xyz format
!
      IMPLICIT NONE
#include "Molcas.fh"
#include "MD.fh"
#include "constants2.fh"
      External IsFreeUnit
      INTEGER     i,j,file,natom,IsFreeUnit
      REAL*8      xyz(natom*3)
      CHARACTER   atom(natom)*2,filename*9
      LOGICAL     hybrid,Exist
!
      IF (.NOT.hybrid) THEN
         file=82
         file=IsFreeUnit(file)
         filename='md.xyz'
         CALL OpnFl(filename,file,Exist)
         CALL Append_file(file)
         WRITE(file,'(I5,/)') natom
         DO i=1, natom
            WRITE(file,100) atom(i),(Angstrom*xyz(3*(i-1)+j),j=1,3)
         END DO
         CLOSE(file)
      ELSE
         file=IsFreeUnit(file)
         filename='md.prmcrd'
         CALL OpnFl(filename,file,Exist)
         CALL Append_file(file)
         WRITE(file,'(/,I6)') natom
         WRITE(file,'(6F12.7)') (Angstrom*xyz(i),i=1,3*natom)
         CLOSE(file)
         file=IsFreeUnit(file)
         filename='vmd.mdcrd'
         CALL OpnFl(filename,file,Exist)
         CALL Append_file(file)
         WRITE(file,'(10F8.3)') (Angstrom*xyz(i),i=1,3*natom)
         CLOSE(file)
      ENDIF
 100  FORMAT(1X,A2,3F15.8)
      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE DxEnergies(time,Epot,Ekin,Etot)
!
!     Write out a summary of energies in CSV format
!
      IMPLICIT NONE
      External IsFreeUnit
#include "WrkSpc.fh"
      INTEGER      file,nEnergies,i,n,ipEnergies
      INTEGER      IsFreeUnit
      REAL*8       time,Epot,Ekin,Etot
      CHARACTER    filename*12,frmt*24
      LOGICAL      Exist,RootCheck
!
      filename='md.energies'
!
      RootCheck=.False.
!
      CALL Qpg_iScalar('Relax CASSCF root',RootCheck)
!
      file=82
      file=IsFreeUnit(file)
      CALL OpnFl(filename,file,Exist)
      CALL Append_file(file)
      IF (.NOT.Exist) THEN
         WRITE(file,'(3X,A4,10X,A4,2(17X,A4))') 'time','Epot',          &
     &                                          'Ekin','Etot'
      END IF
      frmt='(F8.2, (2X,D19.12))'
!
      IF (RootCheck) THEN
         CALL Get_iScalar('Number of roots',nEnergies)
         CALL GetMem('MS energies','ALLO','REAL',ipEnergies,nEnergies)
         CALL Get_dArray('Last energies',Work(ipEnergies),nEnergies)
         n = nEnergies + 3
         WRITE(frmt(7:7),'(I1)') n
         WRITE(file,frmt) time,Epot,Ekin,Etot,                          &
     &                   (Work(ipEnergies-1+i),i=1,nEnergies)
         CALL GetMem('MS energies','FREE','REAL',ipEnergies,nEnergies)
      ELSE
         n = 3
         WRITE(frmt(7:7),'(I1)') n
         WRITE(file,frmt) time,Epot,Ekin,Etot
      END IF
      CLOSE(file)
      RETURN

      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE Put_Velocity(vel,natom3)
      IMPLICIT NONE
#include "WrkSpc.fh"
#include "Molcas.fh"
      INTEGER      natom3
      REAL*8       vel(natom3)
!
!     Writes the velocities on RUNFILE
!
      CALL Put_dArray('Velocities',vel,natom3)
!      CALL GetMem('Velocities','FREE','REAL',vel,natom3)
      RETURN
!
      END
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE Get_Velocity(vel,natom3)
      IMPLICIT NONE
#include "WrkSpc.fh"
#include "Molcas.fh"
      INTEGER      natom3
      REAL*8       vel(natom3)
!
!     Reads the velocities from RUNFILE
!
!      CALL GetMem('Velocities','ALLO','REAL',vel,natom3)
      CALL Get_dArray('Velocities',vel,natom3)
!
      RETURN
!
      END
!

!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE Get_NHC(NHC,nh)
      IMPLICIT NONE
#include "WrkSpc.fh"
#include "Molcas.fh"
      INTEGER      nh
      REAL*8       NHC(nh)
!
!     Reads the extra degrees of freedom from RUNFILE
!
!      CALL GetMem('NOSEHOOVER','ALLO','REAL',NHC,nh)
      CALL Get_dArray('NOSEHOOVER',NHC,nh)
!
      RETURN
!
      END
!
!
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
!
      SUBROUTINE Put_NHC(NHC,nh)
      IMPLICIT NONE
#include "WrkSpc.fh"
#include "Molcas.fh"
      INTEGER      nh
      REAL*8       NHC(nh)
!
!     Writes the extra degrees of Freedom on RUNFILE
!
      CALL Put_dArray('NOSEHOOVER',NHC,nh)
!      CALL GetMem('NOSEHOOVER','FREE','REAL',NHC,nh)
      RETURN
!
      END
!
