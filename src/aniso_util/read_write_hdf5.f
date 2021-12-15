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
#ifdef _HDF5_
      SUBROUTINE read_hdf5_init(file_h5,nstate,nss)
      Implicit None
#include "stdalloc.fh"
#include "mh5.fh"
      Character(Len=180),intent(in) ::    file_h5
      Integer, intent(out)      ::    nstate,nss
      ! local variables:
      Integer                   ::    i,fileid
      Character                 ::    tmp*256, sFile*128
      Character(Len=180)        ::    tmp2
      Integer, allocatable      ::    spin_mult(:)
      Character(Len=5)          ::    molcas_module_kind
      Logical                   ::    Exist
      Logical                   ::    DBG

      DBG=.false.

      WRITE (6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)
      NSS=0
      NSTATE=0
      Exist=.false.
      ! Check if it is a rassi-h5 file:

      ! Check the file exists in $WorkDir
      Call f_inquire(trim(file_h5),Exist)
      IF (Exist) WRITE (6,*) 'file ',trim(file_h5),' exists!!!'
      ! if not present, look for it in the submit directory:
      IF (.not.Exist) THEN
        Call getenvf('MOLCAS_SUBMIT_DIR',tmp)
        IF (tmp.ne.' ') THEN
          i=index(tmp,' ')
          IF (i.gt.0) THEN
            sFile=trim(tmp(1:i-1)//'/'//file_h5)
            Call f_inquire(sFile,Exist)
          END IF
        END IF
        ! if still not present, warn the user and abort:
        IF (.not.Exist) THEN
          Call WarningMessage(2,'File '//
     &                        trim(file_h5)//' is not found')
          Call Quit_OnUserError()
        END IF
      END IF
!----------------------------------------------------------------------|
      ! open the file
      fileid = mh5_open_file_r(trim(file_h5))
      ! check if it is an HDF5 file produced by RASSI
      Call mh5_fetch_attr(fileid,'MOLCAS_MODULE',tmp2)
      molcas_module_kind=trim(tmp2)
      IF (DBG) WRITE (6,'(A,A)') 'read_hdf5::  molcas_module=',
     &                           molcas_module_kind
      IF ( molcas_module_kind(1:5).ne.'RASSI') THEN
          Call WarningMessage(2,'Input HDF5 file '//trim(file_h5)//
     &                        ' is not produced by RASSI')
          Call Quit_OnUserError()
      END IF

!----------------------------------------------------------------------|
      ! read number of spin free states
      Call mh5_fetch_attr(fileid,'NSTATE',nstate)
      IF (DBG) WRITE (6,'(A,I6)') 'read_hdf5::  NSTATE=', NSTATE
      Call Put_iScalar('NSTATE_SINGLE   ',NSTATE)


      ! read spin multiplicity of each state:
      Call mma_allocate(spin_mult,nstate,'nstate')
      Call mh5_fetch_attr(fileid,'STATE_SPINMULT',spin_mult(1:nstate))
      IF (DBG) WRITE (6,'(A)') 'spin_mult'
      IF (DBG) WRITE (6,'(20I4)') (spin_mult(i),i=1,nstate)
      ! compute the number of spin-orbit states:
      nss=0
      DO i=1,nstate
        nss=nss+spin_mult(i)
      END DO
      IF (DBG) WRITE (6,'(A,I6)') 'read_hdf5::     NSS=',NSS
      Call Put_iScalar('NSS_SINGLE      ',NSS)
      Call mma_deallocate(spin_mult)

!----------------------------------------------------------------------|
      ! close the file
      Call mh5_close_file(fileid)

      RETURN
      END SUBROUTINE read_hdf5_init






      SUBROUTINE read_hdf5_all(file_h5, nss, nstate,
     &                         multiplicity, eso, esfs,
     &                         U, MM, MS, ML, DM, ANGMOM, EDMOM,
     &                         amfi, HSO )
      Implicit None
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
#include "mh5.fh"
      Integer, intent(in)           :: nstate,nss
      Integer, intent(out)          :: multiplicity(nstate)

      Real(kind=8), intent(out)    :: esfs(nstate)
      Real(kind=8), intent(out)    :: eso(nss)
      Real(kind=8), intent(out)    ::  edmom(3,nstate,nstate)
      Real(kind=8), intent(out)    ::   amfi(3,nstate,nstate)
      Real(kind=8), intent(out)    :: angmom(3,nstate,nstate)
      Complex(kind=8), intent(out) :: MM(3,nss,nss)
      Complex(kind=8), intent(out) :: MS(3,nss,nss)
      Complex(kind=8), intent(out) :: ML(3,nss,nss)
!     electric dipole moment
      Complex(kind=8), intent(out) :: DM(3,nss,nss)
      Complex(kind=8), intent(out) :: U(nss,nss)
      Complex(kind=8), intent(out) :: HSO(nss,nss)

      Real(kind=8)                 :: AU2CM
      Real(kind=8), allocatable    :: etmp(:)
      Real(kind=8), allocatable    :: RR(:,:), RI(:,:)
      Real(kind=8), allocatable    :: AL(:,:,:)
      Integer                       :: fileid,jend,INRM
      Character(Len=180)            :: file_h5
      Real(kind=8)                 :: RNRM
      Real(kind=8), external       :: dnrm2_, dznrm2_
      Complex(kind=8), external    :: spin
      ! local variables:
      Integer       :: iss, ibas(nstate,-50:50)
      Integer       :: i, j, i1, j1, ist, jst, mult, multI, multJ
      Integer       :: l, ipar
      Real(kind=8) :: g_e
      Complex(kind=8), allocatable    :: tmp(:,:)

!      Logical :: Exist
      Logical :: found_edmom, found_angmom, found_hso, found_amfi,
     &           found_sos_coeff, found_eso, found_esfs, found_mult
      Logical :: DBG

      DBG  =.false.
      AU2CM=219474.6313702_wp
      g_e=2.00231930437180_wp
      found_edmom=.false.
      found_angmom=.false.
      found_hso=.false.
      found_amfi=.false.
      found_sos_coeff=.false.
      found_eso=.false.
      found_esfs=.false.
      found_mult=.false.

      WRITE (6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)

      ! the presence of the RASSI-HDF5 has already been made in
      ! read_hdf5_init subroutine, executed earlier
      ! nss, nstate are already known
!----------------------------------------------------------------------|
      ! open the file
      fileid = mh5_open_file_r(trim(file_h5))
      IF (DBG) WRITE (6,'(A,I24)') 'read_hdf5_all:: fileid=',fileid
!----------------------------------------------------------------------|
      ! read spin multiplicity of each state:
c      IF (mh5_exists_dset(fileid,'STATE_SPINMULT')) THEN
         found_mult=.true.
         Call mh5_fetch_attr(fileid,'STATE_SPINMULT',
     &                       multiplicity(1:nstate))
         IF (DBG) WRITE (6,'(A)') 'read_hdf5_all:: multiplicity'
         IF (DBG) WRITE (6,'(20I4)') (multiplicity(i),i=1,nstate)
         INRM=0
         INRM=SUM(multiplicity(:))
         IF ( INRM == 0 ) THEN
            Call WarningMessage(2,'STATE_SPINMULT array read from '//
     &                          'HDF5 file has norm = zero')
            WRITE (6,*) 'Norm=',INRM
         END IF
c      ELSE
c         Call WarningMessage(2,'State multiplicity array was not '//
c     &                         'found on HDF5 file')
c      END IF

!----------------------------------------------------------------------|
      ! read the spin free energies (cm-1)
      Call mma_allocate(etmp,nstate,'etmp)')
      Call dcopy_(nstate,[0.0_wp],0,etmp,1)
      IF (mh5_exists_dset(fileid,'SFS_ENERGIES')) THEN
         found_esfs=.true.
         Call mh5_fetch_dset_array_real(fileid,'SFS_ENERGIES',etmp)
         RNRM=0.0_wp
         RNRM=dnrm2_(nstate,etmp,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'ESFS read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         ! compute the energeis in cm-1:
         DO i=1,nstate
            esfs(i)=(etmp(i)-etmp(1))*AU2CM
         END DO
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: esfs'
            DO i=1,nstate,4
               jEND=MIN(nstate,i+3)
               WRITE (6,'(4ES24.14)') (esfs(j),j=i,jEnd)
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'Spin-free energies were not '//
     &                         'found on HDF5 file')
      END IF
      Call mma_deallocate(etmp)


!----------------------------------------------------------------------|
      ! read the spin-orbit energies (cm-1)
      Call mma_allocate(etmp,nss,'tmp)')
      Call dcopy_(nss,[0.0_wp],0,etmp,1)
      IF (mh5_exists_dset(fileid,'SOS_ENERGIES')) THEN
         found_eso=.true.
         Call mh5_fetch_dset_array_real(fileid,'SOS_ENERGIES',etmp)
         RNRM=0.0_wp
         RNRM=dnrm2_(nss,etmp,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'ESO read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         ! compute the energeis in cm-1:
         DO i=1,nss
            eso(i)=(etmp(i)-etmp(1))*AU2CM
         END DO
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: eso'
            DO i=1,nss,4
               jEnd=MIN(nss,i+3)
               WRITE (6,'(4ES24.14)') (eso(j),j=i,jEnd)
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'Spin-orbit energies were not '//
     &                         'found on HDF5 file')
      END IF
      Call mma_deallocate(etmp)


!----------------------------------------------------------------------|
      ! read the spin-orbit mixing coefficient matrix:
      Call mma_allocate(RR,nss,nss,'RR')
      Call mma_allocate(RI,nss,nss,'RI')
      Call dcopy_(nss*nss,[0.0_wp],0,RR,1)
      Call dcopy_(nss*nss,[0.0_wp],0,RI,1)

      IF ( mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL').and.
     &     mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG') ) THEN
         found_sos_coeff=.true.
         Call mh5_fetch_dset_array_real(fileid,'SOS_COEFFICIENTS_REAL',
     &                                          RR )
         Call mh5_fetch_dset_array_real(fileid,'SOS_COEFFICIENTS_IMAG',
     &                                          RI )
         ! assemble the complex matrix U:
         DO i=1,nss
            DO j=1,nss
               U(i,j)=cmplx( RR(i,j), RI(i,j), wp )
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dznrm2_(nss*nss,U,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'SOS-U matrix read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: U'
            DO i=1,nss
               DO j=1,nss
                  WRITE (6,'(2i4,A,2ES24.14)') i,j,' |',U(i,j)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'SO mixing coefficeints were not '//
     &                         'found on HDF5 file')
      END IF
      Call mma_deallocate(RR)
      Call mma_deallocate(RI)


!----------------------------------------------------------------------|
      ! read the angular momentum integrals (L):
      Call mma_allocate(AL,nstate,nstate,3,'AL')
      Call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
      IF (mh5_exists_dset(fileid,'SFS_ANGMOM')) THEN
         found_angmom=.true.
         Call mh5_fetch_dset_array_real(fileid,'SFS_ANGMOM',AL)
         DO i=1,nstate
            DO j=1,nstate
               DO l=1,3
                  ANGMOM(l,i,j)=AL(i,j,l)
               END DO
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dnrm2_(nstate*nstate*3,AL,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'SFS_ANGMOM read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: ANGMOM (x,y,z)'
            DO i=1,nstate
               DO j=1,nstate
                  WRITE (6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (ANGMOM(l,i,j),l=1,3)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'ANGMOM integrals were not found on '//
     &                         'HDF5 file')
      END IF
      Call mma_deallocate(AL)


!----------------------------------------------------------------------|
      ! read the electric dipole momentum integrals (EDMOM):
      Call mma_allocate(AL,nstate,nstate,3,'AL')
      Call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
      IF (mh5_exists_dset(fileid,'SFS_EDIPMOM')) THEN
         found_edmom=.true.
         Call mh5_fetch_dset_array_real(fileid,'SFS_EDIPMOM',AL)
         DO i=1,nstate
            DO j=1,nstate
               DO l=1,3
                  EDMOM(l,i,j)=AL(i,j,l)
               END DO
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dnrm2_(nstate*nstate*3,AL,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'SFS_EDIPMOM read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: SFS_EDIPMOM(x,y,z)'
            DO i=1,nstate
               DO j=1,nstate
                  WRITE (6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (EDMOM(l,i,j),l=1,3)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'EDMOM integrals were not found on '//
     &                         'HDF5 file')
      END IF
      Call mma_deallocate(AL)


!----------------------------------------------------------------------|
      ! read the spin-orbit integrals (AMFI):
      Call mma_allocate(AL,nstate,nstate,3,'AL')
      Call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
      IF (mh5_exists_dset(fileid,'SFS_AMFIINT')) THEN
         found_amfi=.true.
         Call mh5_fetch_dset_array_real(fileid,'SFS_AMFIINT',AL)
         DO i=1,nstate
            DO j=1,nstate
               DO l=1,3
                  AMFI(l,i,j)=AL(i,j,l)
               END DO
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dnrm2_(nstate*nstate*3,AL,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'SFS_AMFIINT read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: SFS_AMFIINT(x,y,z)'
            DO i=1,nstate
               DO j=1,nstate
                  WRITE (6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (AMFI(l,i,j),l=1,3)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'AMFI integrals were not found on '//
     &                         'HDF5 file')
      END IF
      Call mma_deallocate(AL)

!----------------------------------------------------------------------|
      ! read the RASSI SFS Hamiltonian (SFS_HAM):


!----------------------------------------------------------------------|
      ! read the RASSI SOS Hamiltonian (SOS_HAM):
      Call mma_allocate(RR,nss,nss,'RR')
      Call mma_allocate(RI,nss,nss,'RI')
      Call dcopy_(nss*nss,[0.0_wp],0,RR,1)
      Call dcopy_(nss*nss,[0.0_wp],0,RI,1)

      IF ( mh5_exists_dset(fileid,'HSO_MATRIX_REAL').and.
     &     mh5_exists_dset(fileid,'HSO_MATRIX_IMAG') ) THEN
         found_hso=.true.
         Call mh5_fetch_dset_array_real(fileid,'HSO_MATRIX_REAL',RR )
         Call mh5_fetch_dset_array_real(fileid,'HSO_MATRIX_IMAG',RI )
         ! assemble the complex matrix U:
         DO i=1,nss
            DO j=1,nss
               HSO(i,j)=cmplx( RR(i,j), RI(i,j), wp )
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dznrm2_(nss*nss,HSO,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'HSO matrix read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: HSO'
            DO i=1,nss
               DO j=1,nss
                  WRITE (6,'(2i4,A,2ES24.14)') i,j,' |',HSO(i,j)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'HSO matrix was not found on HDF5 file')
      END IF
      Call mma_deallocate(RR)
      Call mma_deallocate(RI)


!----------------------------------------------------------------------|
!  All info has been read
!  transform the SFS data to the SO basis
!----------------------------------------------------------------------|
      ! generate a local indexing table:
      iss=0
      ibas=0
      ipar=mod(multiplicity(1),2)
      DO Ist=1,nstate
         Mult=Multiplicity(Ist)
         DO I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            IF ( (Ipar == 0) .AND. (I == 0)) GO TO 310
               Iss=Iss+1
               Ibas(Ist,I)=Iss
  310       CONTINUE
         END DO ! i
      END DO ! ist



c----- expand the spin free basis to the spin-orbit basis:
      MM=(0.0_wp,0.0_wp)
      ML=(0.0_wp,0.0_wp)
      MS=(0.0_wp,0.0_wp)
      DO Ist=1,nstate
         Mult=Multiplicity(Ist)
         DO I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            IF ( (Ipar == 0) .AND. (I == 0) ) GO TO 301
            DO J=-(Mult-Ipar)/2,(Mult-Ipar)/2
               IF ( (Ipar == 0) .AND. (J == 0) ) GO TO 302
               DO l=1,3
                  i1=Ibas(Ist,I)
                  j1=Ibas(Ist,J)
                  MM(l,i1,j1)=-Spin(l,Mult,I,J)*g_e
                  MS(l,i1,j1)= Spin(l,Mult,I,J)
               END DO ! l
  302          CONTINUE
            END DO ! J
  301       CONTINUE
         END DO ! I
      END DO ! Ist


      DO Ist=1,nstate
         MultI=Multiplicity(Ist)
         DO Jst=1,nstate
            MultJ=Multiplicity(Jst)
            IF (MultI == MultJ) THEN
               DO I=-(MultI-Ipar)/2,(MultI-Ipar)/2
               IF ( (Ipar == 0) .AND. (I == 0) ) GO TO 303
                  DO l=1,3
                     i1=Ibas(Ist,I)
                     j1=Ibas(Jst,I)
                     MM(l,i1,j1)=MM(l,i1,j1)
     &                              -CMPLX(0.0_wp,Angmom(l,Ist,Jst),wp)
                     ML(l,i1,j1)=ML(l,i1,j1)
     &                              +CMPLX(0.0_wp,Angmom(l,Ist,Jst),wp)
                     DM(l,i1,j1)=DM(l,i1,j1)
     &                              +CMPLX(eDmom(l,Ist,Jst),0.0_wp,wp)
                  END DO   ! l
  303             CONTINUE
               END DO   ! I
            END IF
         END DO   ! Jst
      END DO   ! Ist

      ! calculate the matrix elements of the spin and magnetic moment
      ! in the spin-orbit basis:
      Call mma_allocate(tmp,nss,nss,'tmp')
      DO L=1,3
         TMP=(0.0_wp,0.0_wp)
         ! spin moment
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             MS(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             MS(L,:,:), nss )
         ! orbital moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             ML(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             ML(L,:,:), nss )
         ! magnetic moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             MM(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             MM(L,:,:), nss )

         IF (found_EDMOM) THEN
         ! electric dipole moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             DM(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             DM(L,:,:), nss )
         END IF
      END DO !L
      Call mma_deallocate(tmp)




      ! close the file
      Call mh5_close_file(fileid)

      RETURN
      END SUBROUTINE read_hdf5_all







      SUBROUTINE read_hdf5_poly( file_h5, nss, nstate,
     &                           eso, MM, MS,iReturn )
      Implicit None
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
#include "mh5.fh"
      Integer, intent(in)           :: nstate,nss
      Integer                       :: multiplicity(nstate)
      Integer                       :: iReturn

!      Real(kind=8), intent(out)    :: esfs(nstate)
      Real(kind=8), intent(out)    :: eso(nss)
!      Real(kind=8), intent(out)    ::  edmom(3,nstate,nstate)
!      Real(kind=8), intent(out)    ::   amfi(3,nstate,nstate)
      Real(kind=8)                 :: angmom(3,nstate,nstate)
      Complex(kind=8)  :: MM(3,nss,nss)
      Complex(kind=8)  :: MS(3,nss,nss)
!      Complex(kind=8), intent(out) :: ML(3,nss,nss)
!      Complex(kind=8), intent(out) :: DM(3,nss,nss) ! electric dipole moment
!      Complex(kind=8)              :: U(nLoc,nLoc)
!      Complex(kind=8), intent(out) :: HSO(nss,nss)

      Real(kind=8)                 :: AU2CM
      Real(kind=8), allocatable    :: etmp(:)
      Real(kind=8), allocatable    :: RR(:,:), RI(:,:)
      Real(kind=8), allocatable    :: AL(:,:,:)
      Complex(kind=8), allocatable :: U(:,:)
      Integer                       :: fileid,jend,INRM
      Character(Len=180)            :: file_h5
      Real(kind=8)                 :: RNRM
      Real(kind=8), external       :: dnrm2_, dznrm2_
      Complex(kind=8), external    :: spin
      ! local variables:
      Integer       :: iss, ibas(nstate,-50:50)
      Integer       :: i, j, i1, j1, ist, jst, mult, multI, multJ
      Integer       :: l, ipar
      Real(kind=8) :: g_e
      Complex(kind=8), allocatable    :: tmp(:,:)

!      Logical :: Exist
      Logical :: found_edmom, found_angmom, found_hso, found_amfi,
     &           found_sos_coeff, found_eso, found_esfs, found_mult
      Logical :: DBG

      DBG  =.false.
      AU2CM=219474.6313702_wp
      g_e=2.00231930437180_wp
      found_edmom=.false.
      found_angmom=.false.
      found_hso=.false.
      found_amfi=.false.
      found_sos_coeff=.false.
      found_eso=.false.
      found_esfs=.false.
      found_mult=.false.
      iReturn=0

      WRITE (6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)

      ! the presence of the RASSI-HDF5 has already been made in
      ! read_hdf5_init subroutine, executed earlier
      ! nss, nstate are already known
!----------------------------------------------------------------------|
      ! open the file
      fileid = mh5_open_file_r(trim(file_h5))
      IF (DBG) WRITE (6,'(A,I24)') 'read_hdf5_all:: fileid=',fileid
!----------------------------------------------------------------------|
      ! read spin multiplicity of each state:
c      IF (mh5_exists_dset(fileid,'STATE_SPINMULT')) THEN
         found_mult=.true.
         Call mh5_fetch_attr(fileid,'STATE_SPINMULT',
     &                       multiplicity(1:nstate))
         IF (DBG) WRITE (6,'(A)') 'read_hdf5_all:: multiplicity'
         IF (DBG) WRITE (6,'(20I4)') (multiplicity(i),i=1,nstate)
         INRM=0
         INRM=SUM(multiplicity(:))
         IF ( INRM == 0 ) THEN
            Call WarningMessage(2,'STATE_SPINMULT array read from '//
     &                          'HDF5 file has norm = zero')
            WRITE (6,*) 'Norm=',INRM
         END IF
c      ELSE
c         Call WarningMessage(2,'State multiplicity array was not '//
c     &                         'found on HDF5 file')
c      END IF

!----------------------------------------------------------------------|
!      ! read the spin free energies (cm-1)
!      Call mma_allocate(etmp,nstate,'etmp)')
!      Call dcopy_(nstate,[0.0_wp],0,etmp,1)
!      IF (mh5_exists_dset(fileid,'SFS_ENERGIES')) THEN
!         found_esfs=.true.
!         Call mh5_fetch_dset_array_real(fileid,'SFS_ENERGIES',etmp)
!         RNRM=0.0_wp
!         RNRM=dnrm2_(nstate,etmp,1)
!         IF ( RNRM.lt.1.0D-50 ) THEN
!            Call WarningMessage(2,'ESFS read from HDF5 file '//
!     &                          ' has norm = zero')
!            WRITE (6,*) 'Norm=',RNRM
!         END IF
!         ! compute the energeis in cm-1:
!         DO i=1,nstate
!            esfs(i)=(etmp(i)-etmp(1))*AU2CM
!         END DO
!         IF (DBG) THEN
!            WRITE (6,'(A)') 'read_hdf5_all:: esfs'
!            DO i=1,nstate,4
!               jEND=MIN(nstate,i+3)
!               WRITE (6,'(4ES24.14)') (esfs(j),j=i,jEnd)
!            END DO
!         END IF
!      ELSE
!         Call WarningMessage(2,'Spin-free energies were not '//
!     &                         'found on HDF5 file')
!      END IF
!      Call mma_deallocate(etmp)
!
!
!----------------------------------------------------------------------|
      ! read the spin-orbit energies (cm-1)
      Call mma_allocate(etmp,nss,'tmp)')
      Call dcopy_(nss,[0.0_wp],0,etmp,1)
      IF (mh5_exists_dset(fileid,'SOS_ENERGIES')) THEN
         found_eso=.true.
         Call mh5_fetch_dset_array_real(fileid,'SOS_ENERGIES',etmp)
         RNRM=0.0_wp
         RNRM=dnrm2_(nss,etmp,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'ESO read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         ! compute the energeis in cm-1:
         DO i=1,nss
            eso(i)=(etmp(i)-etmp(1))*AU2CM
         END DO
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: eso'
            DO i=1,nss,4
               jEnd=MIN(nss,i+3)
               WRITE (6,'(4ES24.14)') (eso(j),j=i,jEnd)
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'Spin-orbit energies were not '//
     &                         'found on HDF5 file')
      END IF
      Call mma_deallocate(etmp)


!----------------------------------------------------------------------|
      ! read the spin-orbit mixing coefficient matrix:
      Call mma_allocate(RR,nss,nss,'RR')
      Call mma_allocate(RI,nss,nss,'RI')
      Call mma_allocate(U,nss,nss,'U')
      Call dcopy_(nss*nss,[0.0_wp],0,RR,1)
      Call dcopy_(nss*nss,[0.0_wp],0,RI,1)
      Call zcopy_(nss*nss,[(0.0_wp,0.0_wp)],0,U,1)

      IF ( mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL').and.
     &     mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG') ) THEN
         found_sos_coeff=.true.
         Call mh5_fetch_dset_array_real(fileid,'SOS_COEFFICIENTS_REAL',
     &                                          RR )
         Call mh5_fetch_dset_array_real(fileid,'SOS_COEFFICIENTS_IMAG',
     &                                          RI )
         ! assemble the complex matrix U:
         DO i=1,nss
            DO j=1,nss
               U(i,j)=cmplx( RR(i,j), RI(i,j), wp )
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dznrm2_(nss*nss,U,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'SOS-U matrix read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: U'
            DO i=1,nss
               DO j=1,nss
                  WRITE (6,'(2i4,A,2ES24.14)') i,j,' |',U(i,j)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'SO mixing coefficeints were not '//
     &                         'found on HDF5 file')
      END IF
      Call mma_deallocate(RR)
      Call mma_deallocate(RI)


!----------------------------------------------------------------------|
      ! read the angular momentum integrals (L):
      Call mma_allocate(AL,nstate,nstate,3,'AL')
      Call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
      IF (mh5_exists_dset(fileid,'SFS_ANGMOM')) THEN
         found_angmom=.true.
         Call mh5_fetch_dset_array_real(fileid,'SFS_ANGMOM',AL)
         DO i=1,nstate
            DO j=1,nstate
               DO l=1,3
                  ANGMOM(l,i,j)=AL(i,j,l)
               END DO
            END DO
         END DO
         RNRM=0.0_wp
         RNRM=dnrm2_(nstate*nstate*3,AL,1)
         IF ( RNRM.lt.1.0D-50 ) THEN
            Call WarningMessage(2,'SFS_ANGMOM read from HDF5 file '//
     &                          ' has norm = zero')
            WRITE (6,*) 'Norm=',RNRM
         END IF
         IF (DBG) THEN
            WRITE (6,'(A)') 'read_hdf5_all:: ANGMOM (x,y,z)'
            DO i=1,nstate
               DO j=1,nstate
                  WRITE (6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (ANGMOM(l,i,j),l=1,3)
               END DO
            END DO
         END IF
      ELSE
         Call WarningMessage(2,'ANGMOM integrals were not found on '//
     &                         'HDF5 file')
      END IF
      Call mma_deallocate(AL)




!----------------------------------------------------------------------|
!      ! read the electric dipole momentum integrals (EDMOM):
!      Call mma_allocate(AL,nstate,nstate,3,'AL')
!      Call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
!      IF (mh5_exists_dset(fileid,'SFS_EDIPMOM')) THEN
!         found_edmom=.true.
!         Call mh5_fetch_dset_array_real(fileid,'SFS_EDIPMOM',AL)
!         DO i=1,nstate
!            DO j=1,nstate
!               DO l=1,3
!                  EDMOM(l,i,j)=AL(i,j,l)
!               END DO
!            END DO
!         END DO
!         RNRM=0.0_wp
!         RNRM=dnrm2_(nstate*nstate*3,AL,1)
!         IF ( RNRM.lt.1.0D-50 ) THEN
!            Call WarningMessage(2,'SFS_EDIPMOM read from HDF5 file '//
!     &                          ' has norm = zero')
!            WRITE (6,*) 'Norm=',RNRM
!         END IF
!         IF (DBG) THEN
!            WRITE (6,'(A)') 'read_hdf5_all:: SFS_EDIPMOM(x,y,z)'
!            DO i=1,nstate
!               DO j=1,nstate
!                  WRITE (6,'(2i4,A,3ES24.14)') i,j,' |',
!     &                                       (EDMOM(l,i,j),l=1,3)
!               END DO
!            END DO
!         END IF
!      ELSE
!         Call WarningMessage(2,'EDMOM integrals were not found on '//
!     &                         'HDF5 file')
!      END IF
!      Call mma_deallocate(AL)
!
!
!----------------------------------------------------------------------|
!      ! read the spin-orbit integrals (AMFI):
!      Call mma_allocate(AL,nstate,nstate,3,'AL')
!      Call dcopy_(nstate*nstate*3,[0.0_wp],0,AL,1)
!      IF (mh5_exists_dset(fileid,'SFS_AMFIINT')) THEN
!         found_amfi=.true.
!         Call mh5_fetch_dset_array_real(fileid,'SFS_AMFIINT',AL)
!         DO i=1,nstate
!            DO j=1,nstate
!               DO l=1,3
!                  AMFI(l,i,j)=AL(i,j,l)
!               END DO
!            END DO
!         END DO
!         RNRM=0.0_wp
!         RNRM=dnrm2_(nstate*nstate*3,AL,1)
!         IF ( RNRM.lt.1.0D-50 ) THEN
!            Call WarningMessage(2,'SFS_AMFIINT read from HDF5 file '//
!     &                          ' has norm = zero')
!            WRITE (6,*) 'Norm=',RNRM
!         END IF
!         IF (DBG) THEN
!            WRITE (6,'(A)') 'read_hdf5_all:: SFS_AMFIINT(x,y,z)'
!            DO i=1,nstate
!               DO j=1,nstate
!                  WRITE (6,'(2i4,A,3ES24.14)') i,j,' |',
!     &                                       (AMFI(l,i,j),l=1,3)
!               END DO
!            END DO
!         END IF
!      ELSE
!         Call WarningMessage(2,'AMFI integrals were not found on '//
!     &                         'HDF5 file')
!      END IF
!      Call mma_deallocate(AL)
!
!----------------------------------------------------------------------|
      ! read the RASSI SFS Hamiltonian (SFS_HAM):


!----------------------------------------------------------------------|
!      ! read the RASSI SOS Hamiltonian (SOS_HAM):
!      Call mma_allocate(RR,nss,nss,'RR')
!      Call mma_allocate(RI,nss,nss,'RI')
!      Call dcopy_(nss*nss,[0.0_wp],0,RR,1)
!      Call dcopy_(nss*nss,[0.0_wp],0,RI,1)
!
!      IF ( mh5_exists_dset(fileid,'HSO_MATRIX_REAL').and.
!     &     mh5_exists_dset(fileid,'HSO_MATRIX_IMAG') ) THEN
!         found_hso=.true.
!         Call mh5_fetch_dset_array_real(fileid,'HSO_MATRIX_REAL',RR )
!         Call mh5_fetch_dset_array_real(fileid,'HSO_MATRIX_IMAG',RI )
!         ! assemble the complex matrix U:
!         DO i=1,nss
!            DO j=1,nss
!               HSO(i,j)=cmplx( RR(i,j), RI(i,j), wp )
!            END DO
!         END DO
!         RNRM=0.0_wp
!         RNRM=dznrm2_(nss*nss,HSO,1)
!         IF ( RNRM.lt.1.0D-50 ) THEN
!            Call WarningMessage(2,'HSO matrix read from HDF5 file '//
!     &                          ' has norm = zero')
!            WRITE (6,*) 'Norm=',RNRM
!         END IF
!         IF (DBG) THEN
!            WRITE (6,'(A)') 'read_hdf5_all:: HSO'
!            DO i=1,nss
!               DO j=1,nss
!                  WRITE (6,'(2i4,A,2ES24.14)') i,j,' |',HSO(i,j)
!               END DO
!            END DO
!         END IF
!      ELSE
!         Call WarningMessage(2,'HSO matrix was not found on HDF5 file')
!      END IF
!      Call mma_deallocate(RR)
!      Call mma_deallocate(RI)
!
!
!----------------------------------------------------------------------|
!  All info has been read
!  transform the SFS data to the SO basis
!----------------------------------------------------------------------|
      ! generate a local indexing table:
      iss=0
      ibas=0
      ipar=mod(multiplicity(1),2)
      DO Ist=1,nstate
         Mult=Multiplicity(Ist)
         DO I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            IF ( (Ipar == 0) .AND. (I == 0)) GO TO 310
               Iss=Iss+1
               Ibas(Ist,I)=Iss
  310       CONTINUE
         END DO ! i
      END DO ! ist



c----- expand the spin free basis to the spin-orbit basis:
      MM=(0.0_wp,0.0_wp)
!      ML=(0.0_wp,0.0_wp)
      MS=(0.0_wp,0.0_wp)
      DO Ist=1,nstate
         Mult=Multiplicity(Ist)
         DO I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            IF ( (Ipar == 0) .AND. (I == 0) ) GO TO 301
            DO J=-(Mult-Ipar)/2,(Mult-Ipar)/2
               IF ( (Ipar == 0) .AND. (J == 0) ) GO TO 302
               DO l=1,3
                  i1=Ibas(Ist,I)
                  j1=Ibas(Ist,J)
                  MM(l,i1,j1)=-Spin(l,Mult,I,J)*g_e
                  MS(l,i1,j1)= Spin(l,Mult,I,J)
               END DO ! l
  302          CONTINUE
            END DO ! J
  301       CONTINUE
         END DO ! I
      END DO ! Ist


      DO Ist=1,nstate
         MultI=Multiplicity(Ist)
         DO Jst=1,nstate
            MultJ=Multiplicity(Jst)
            IF (MultI == MultJ) THEN
               DO I=-(MultI-Ipar)/2,(MultI-Ipar)/2
               IF ( (Ipar == 0) .AND. (I == 0) ) GO TO 303
                  DO l=1,3
                     i1=Ibas(Ist,I)
                     j1=Ibas(Jst,I)
                     MM(l,i1,j1)=MM(l,i1,j1)
     &                              -CMPLX(0.0_wp,Angmom(l,Ist,Jst),wp)
!                     ML(l,i1,j1)=ML(l,i1,j1)
!     &                              +CMPLX(0.0_wp,Angmom(l,Ist,Jst),wp)
!                     DM(l,i1,j1)=DM(l,i1,j1)
!     &                              +CMPLX(eDmom(l,Ist,Jst),0.0_wp,wp)
                  END DO   ! l
  303             CONTINUE
               END DO   ! I
            END IF
         END DO   ! Jst
      END DO   ! Ist



      ! calculate the matrix elements of the spin and magnetic moment
      ! in the spin-orbit basis:
      Call mma_allocate(tmp,nss,nss,'tmp')
      DO L=1,3
         TMP=(0.0_wp,0.0_wp)
         ! spin moment
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             MS(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             MS(L,:,:), nss )
!         ! orbital moment
!         TMP=(0.0_wp,0.0_wp)
!         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
!     &                     U, nss,
!     &             ML(L,:,:), nss,           (0.0_wp,0.0_wp),
!     &                   TMP, nss )
!         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
!     &                   TMP, nss,
!     &                     U, nss,           (0.0_wp,0.0_wp),
!     &             ML(L,:,:), nss )
         ! magnetic moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             MM(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             MM(L,:,:), nss )

!         IF (found_EDMOM) THEN
!         ! electric dipole moment
!         TMP=(0.0_wp,0.0_wp)
!         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
!     &                     U, nss,
!     &             DM(L,:,:), nss,           (0.0_wp,0.0_wp),
!     &                   TMP, nss )
!         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
!     &                   TMP, nss,
!     &                     U, nss,           (0.0_wp,0.0_wp),
!     &             DM(L,:,:), nss )
!         END IF
      END DO !L
      Call mma_deallocate(tmp)
      Call mma_deallocate(U)


      ! close the file
      Call mh5_close_file(fileid)

      RETURN
      END  SUBROUTINE read_hdf5_poly





#else
      !
      SUBROUTINE read_hdf5_init()
      RETURN
      END

      SUBROUTINE read_hdf5_all()
      RETURN
      END

      SUBROUTINE read_hdf5_poly()
      RETURN
      END

#endif
