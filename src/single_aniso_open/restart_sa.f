       Subroutine restart_sa( input_to_read, input_file_name,
     &                        nss, nstate )

       Implicit None
       Integer        :: nss, nstate, input_to_read, iDisk
       Integer        :: luaniso
       Character(180) :: input_file_name
       Integer        :: IsFreeUnit
       External       :: IsFreeUnit
       Logical        :: dbg
       Call qEnter('SA_restart')
       dbg=.false.


       If ( input_to_read .eq. 1 ) Then
          ! read the binary file "$Project.aniso":
          luaniso=8
          Call daname(luaniso,'POLYFILE')
          iDisk=0
          Call idafile(luaniso,2,nstate,1,iDisk)
          Call idafile(luaniso,2,nss,1,iDisk)
          Call daclos(luaniso)
          ! put them on RunFile:
          Call Put_iScalar('NSTATE_SINGLE   ',nstate)
          Call Put_iScalar('NSS_SINGLE      ',nss)
          Call Put_iScalar('MXJOB_SINGLE    ',1)
          Call Put_iScalar('NJOB_SINGLE     ',1)



       Else If ( (input_to_read .eq. 2) .OR.
     &           (input_to_read .eq. 4) ) Then
          ! read the ascii formatted "aniso.input" file:
          luaniso = IsFreeUnit(18)
          Call molcas_open(luaniso,input_file_name)
          READ(luaniso,*) nstate, nss
          Close(luaniso)
          ! put them on RunFile:
          Call Put_iScalar('NSTATE_SINGLE   ',nstate)
          Call Put_iScalar('NSS_SINGLE      ',nss)
          Call Put_iScalar('MXJOB_SINGLE    ',1)
          Call Put_iScalar('NJOB_SINGLE     ',1)




       Else If ( input_to_read .eq. 3 ) Then
          If(dbg) Write(6,*) 'restart_sa: file h5=',
     &                        trim(input_file_name)
#ifdef _HDF5_
          ! NSS and NSTATE are also placed on RunFile
          Call read_hdf5_init(input_file_name,nstate,nss)
          If(dbg) Write(6,*) 'restart_sa:    nss=',nss
          If(dbg) Write(6,*) 'restart_sa: nstate=',nstate
          Call Put_iScalar('NSTATE_SINGLE   ',nstate)
          Call Put_iScalar('NSS_SINGLE      ',nss)
          Call Put_iScalar('MXJOB_SINGLE    ',1)
          Call Put_iScalar('NJOB_SINGLE     ',1)
#else
          Write(6,'(A)') 'Warning:: restart option was set to 3:'//
     &                   'i.e. from an HDF5 file'
          Write(6,'(A,A)') 'file id =',trim(input_file_name)
          Call WarningMessage(2,'MOLCAS was compiled without '//
     &                          '_HDF5_ option.')
          Call Quit_OnUserError()
#endif


       Else
          Call WarningMessage(2,'SINGLE_ANISO:: RESTART  '//
     &                          'option is not known.')
          Write(6,'(A,I6)') 'restart_option =',input_to_read
          Write(6,'(A,I6)') 'restart_option can only take integer '//
     &                      'values:'
          Write(6,'(A,I6)') '1 - from binary $Project.aniso'
          Write(6,'(A,I6)') '2 - from formatted file "aniso.input" '//
     &                      '(filename can be given in the input)'
          Write(6,'(A,I6)') '3 - from an HDF5 type file generated '//
     &                      'by RASSI code (filename can be given '//
     &                      'in the input)'
          Write(6,'(A,I6)') '4 - from formatted file "aniso.input" '//
     &                      '(filename can be given in the input) '//
     &                      'in molcas-8.0 format'
          Call Quit_OnUserError()

       End If

       Call qExit('SA_restart')
      Return
      End
