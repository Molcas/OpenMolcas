************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2020, Jie J. Bao                                       *
************************************************************************
*       ***************************************************
*                          Calculating NTO
*       ****************************************************
*        Reference: J. Chem. Phys., 2014, 141, 024106
*        Notation used in the following of the code follows those in the
*        reference mentioned above, especially between equation 53 and
*        54 on page 024106-8.
*        NASHT is the number of active orbitals, originated from this
*        program.
*        Umat is the U matrix, which is the eigenvector  matrix
*        calculated by transition density matrix (TDM)
*        multiplied by its transpose. Vmat is the V matrix, the eigen-
*        vector matrix of a matrix calculated by the transpose
*        multiplied by the TDM.
*        Ueig is the eigenvalen matrix for the U matrix, Veig is
*        that
*        for the V matrix.
*        ONTO is the hole NTO, calculated by multiplying MO matrix with
*        the eigenvector matrix for U matrix. Note that the eigenvector
*        matrix is still named as U. Similar condition is for the
*        particle matirx.
*
*        However, the sets of particle and hole orbitals are switched
*        when I examined the results. So I put the data stored in LONTO
*        as the particle NTO. (Because JOB1 is for the second JobIph
*        file in the input and JOB2 is for the first.)
*
*                                         -------Jie Bao
*                  in Depart. of Chemistry, University of Minnesota, USA
*                                             2018/08/09

      SUBROUTINE NTOCalc(JOB1,JOB2,ISTATE,JSTATE,TRAD,TRASD,ISpin)

      use fortran_strings, only : str
      use stdalloc, only: mma_allocate, mma_deallocate
      use Symmetry_Info, only: nSym=>nIrrep
      use Constants, only: Two
      use rassi_data, only: NASHT,NBST,NISHT,NCMO,NASH,NBASF,NISH,NOSH

      Implicit None

      Integer ISpin,JOB1,JOB2
      Integer IState, jState
      Real*8,DIMENSION(NASHT**2)::TRAD,TRASD

      Character,DIMENSION(2) :: Spin
      INTEGER Iprint,Jprint,I,J,isym
! Printing or looping control
      INTEGER IUseSym,NUseSym,NSupBas,icactorb
      INTEGER,DIMENSION(NBST)  :: OrbUsedSym
      INTEGER,DIMENSION(NASHT+NISHT) :: OrbAct
! CMO Symmetry Contrl
      INTEGER,DIMENSION(NSym) :: NUsedBF,NUseBF,UsetoReal,RealtoUse
! Nr. of basis functions used prior to this symmetry (NUsedBF)
! and used in this symmetry (NUseBF) NSym >= NusedSym
      INTEGER   IOrb
!IOrb is the index  of orbitals.
      INTEGER NSUPCMO
      INTEGER NDge
      INTEGER NScrq
      REAL*8 WGRONK(2)
      INTEGER N_NTO,INFO, I_NTO
      REAL*8 PrintThres,SumEigVal
!     re-organizing orbitals
!     This is to convert active MO sets in any symmetry into a C1 symmetry
      INTEGER NAISHT
      INTEGER, DIMENSION(NISHT+NASHT) :: OrbBas,OrbSym
      !OrbBas() is the number of basis function for IOrb
      !OrbSym() is the index of symmetry/irrep  for IOrb
      ! The strings below should be converted to
      ! character(len=:), allocatable format, but currently
      ! gfortran has problems with this
      CHARACTER (len=128) FILENAME
      CHARACTER (len=8)  NTOType
      CHARACTER (len=9)  STATENAME
      Character*3 lIrrep(8)
      Logical DOTEST
      INTEGER LU
      INTEGER, External:: ISFREEUNIT
      EXTERNAL Molden_interface
      Real*8, allocatable:: SUPCMO1(:), SUPCMO2(:)
      Real*8, allocatable:: ONTO(:), UNTO(:)
      Real*8, allocatable:: CMO1(:), CMO2(:)
      Real*8, allocatable:: UMAT(:), VMAT(:)
      Real*8, allocatable:: UEig(:), VEig(:)
      Real*8, allocatable:: TDM(:), TDMT(:)
      Real*8, allocatable:: Scrq(:)
      Integer, allocatable:: Symfr(:), Symto(:)
      Integer, allocatable:: Indfr(:), Indto(:)

      LU=233

      statename=''
      DoTest=.false.
      PrintThres=1.0D-5
      CALL mma_allocate(CMO1,NCMO,Label='CMO1')
      CALL mma_allocate(CMO2,NCMO,Label='CMO2')
      CALL RDCMO_RASSI(JOB1,CMO1)
      CALL RDCMO_RASSI(JOB2,CMO2)

      if(dotest)then
       write(6,*) 'CMO1 '
       Do I=0,NCMO,5
       write(6,'(2X,5F10.6)')
     & (CMO1(I+IPrint),IPrint=1,MIN(5,NCMO-I))
       End Do
       write(6,*) 'CMO2 '
       Do I=0,NCMO,5
       write(6,'(2X,5F10.6)')
     & (CMO2(I+IPrint),IPrint=1,MIN(5,NCMO-I))
       End Do
      endif

C     Analyzing the symmetry of the wave function
      IUseSym=0
      NSupBas=0
      IOrb=0
      Do ISym=1,NSym
       IF (NASH(ISym).GT.0) Then
        IUseSym=IUseSym+1
        RealtoUse(ISym)=IUseSym
        UsetoReal(IUseSym)=ISym
        NSupBas=NSupBas+NBASF(ISym)
        If (IUseSym.gt.1) THEN
         NUsedBF(IUseSym)=NBASF(UsetoReal(IUseSym-1))+NUsedBF(IUseSym-1)
        Else
         NUsedBF(IUseSym)=0
        END IF
        NUseBF(IUseSym)=NBASF(ISym)
        Do I=1,NOSH(ISym)
         IOrb=IOrb+1
         If (I.gt.NISH(ISym)) Then
          OrbAct(IOrb)=1
         Else
          OrbAct(IOrb)=0
         End If
         OrbBas(IOrb)=NBASF(ISym)
         OrbSym(IOrb)=ISym
         OrbUsedSym(IOrb)=IUseSym
        End Do
       Else
        RealtoUse(ISym)=0
       End If
      End Do
      NSUPCMO=NASHT*NSupBas
      NUseSym=IUseSym
      NAISHT=NASHT+NISHT
      IF (DoTest) Then
       write(6,*) 'Reprinting MO information'
       write(6,*) 'Size of Super-CMO matrix',NSupCMO
       write(6,'(6X,A20,4X,16I4)')
     & 'MO Index',(IOrb,IOrb=1,NAISHT)
       write(6,'(6X,A20,4X,16I4)')
     & 'Irrep Belong to',(OrbSym(IOrb),IOrb=1,NAISHT)
       write(6,'(6X,A20,4X,16I4)')
     & 'Nr. of Basis F',(OrbBas(IOrb),IOrb=1,NAISHT)
       write(6,'(6X,A20,4X,16I4)')
     & 'Act Orbital?',(OrbAct(IOrb),IOrb=1,NAISHT)
       write(6,'(6X,A20,4X,16I4)')
     & 'used basis f',(NUsedBF(OrbUsedSym(IOrb)),IOrb=1,NAISHT)
      End If
C     End of analyzing wave function

C     building up a super-CMO matrix (to be C1-like)
      CALL mma_allocate (SUPCMO1,NSUPCMO,Label='SUPCMO1')
      CALL mma_allocate (SUPCMO2,NSUPCMO,Label='SUPCMO1')
      SUPCMO1(:)=0.0D0
      SUPCMO2(:)=0.0D0
      CALL mma_allocate(ONTO,NSUPCMO,Label='ONTO')
      CALL mma_allocate(UNTO,NSUPCMO,Label='UNTO')
      icactorb=0
      I=0
      Do IOrb=1,NAISHT
       IF (OrbAct(IOrb).eq.1) THEN
         icactorb=icactorb+1
         Do IPrint=1,(OrbBas(IOrb))
          JPrint=IPrint+NUsedBF(OrbUsedSym(IOrb))
          J=I+IPrint-1
          SUPCMO1(icactorb+(JPRINT-1)*NASHT)=CMO1(1+J)
          SUPCMO2(icactorb+(JPRINT-1)*NASHT)=CMO2(1+J)
        End DO
       End IF
       I=I+OrbBas(IOrb)
      End DO
      If (DoTest) Then
       write(6,*)'printing CMO1 in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(SUPCMO1(I+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do


       write(6,*)'printing CMO2 in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(SUPCMO2(I+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do
      End If
C     end of building up the super-CMO matrix
C     Start and initialize spaces
      statename = str(JSTATE)//'_'//str(ISTATE)
      NDge=NASHT**2
      CALL mma_allocate (Umat,NDge,Label='UMat')
      CALL mma_allocate (Vmat,NDge,Label='VMat')
      CALL mma_allocate (Ueig,NDge,Label='Ueig')
      CALL mma_allocate (Veig,NDge,Label='Veig')
      CALL mma_allocate (TDM,NDge,Label='TDM')
      CALL mma_allocate (TDMT,NDge,Label='TDMT')
       write(6,*)
       WRITE(6,'(6X,A)') repeat('*',100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A31,33X,A)')
     &     '*','NATURAL TRANSITION ORBITALS','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,38X,A14,I2,A12,I2 ,30X,A )')
     &'*','BETWEEN STATE ',JSTATE,' AND STATE ',ISTATE,'*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,A)') repeat('*',100)
       write(6,*)
       write(6,*)
      IF (ISpin.eq.1) Then
        N_NTO=1
        Spin(1)='a'
        write (6,'(10X,2a)')'NTO CALCULATION ONLY DONE FOR ALPHA SPIN ',
     &'BECAUSE THE WAVE FUNCTION IS A SINGLET,'
        write (6,'(10X,a)')'SO ALPHA NTOS ARE EQUAL TO BETA ONES'
       Else
        N_NTO=2
        Spin(1)='a'
        Spin(2)='b'
      End IF
      Do I_NTO=1,N_NTO
      Ueig(:)=0.0D0
      Veig(:)=0.0D0
      ONTO(:)=0.0D0
      UNTO(:)=0.0D0
      If (Spin(I_NTO).eq.'a') Then
      Do I=1,Ndge
       TDM(I)=(TRAD(I)+TRASD(I))/Two
      End DO
      else
      Do I=1,Ndge
       TDM(I)=(TRAD(I)-TRASD(I))/Two
      End DO
      End IF
      DO I=1,NASHT
       DO J=1,NASHT
        TDMT(I+NASHT*(J-1))=TDM(J+NASHT*(I-1))
       END DO
      END DO
C     Print out transition density matrix
      If(Spin(I_NTO).eq.'a') Then
      write (FILENAME,fmt='(a,a)')
     &"TDM",trim(adjustl(STATENAME))
      LU=ISFREEUNIT(LU)
      CALL Molcas_Open(LU,FILENAME)
      DO I=1,NASHT
        write (LU,'(5(1X,ES11.4E2))')
     &    (TRAD(NASHT*(I-1)+J),J=1,NASHT)
      END DO
      write (LU,*)
      write (LU,*)
      write (LU,*)
      DO I=1,NASHT
        write (LU,'(5(1X,ES11.4E2))')
     &    (TRASD(NASHT*(I-1)+J),J=1,NASHT)
      END DO
      close (LU)
      End If
C     Generalizing transpose of TDM, TDM_T

C     Calculating T_trans*T
      CALL DGEMM_('n','n',NASHT,NASHT,NASHT,1.0D0,TDMT,NASHT,
     &             TDM,NASHT,0.0D0,Vmat,NASHT)
C     Writing Particle Matrix
      write (FILENAME,fmt='(a,a,a)')
     &"Dhole.",trim(adjustl(STATENAME)),Spin(I_NTO)
      LU=ISFREEUNIT(LU)
      CALL Molcas_Open(LU,FILENAME)
      DO I=1,NASHT
        write (LU,'(10(1X,ES11.4E2))') (Vmat(NASHT*(I-1)+J),J=1,NASHT)
      END DO
      close (LU)
C     Calculating T*T_transpose
      CALL DGEMM_('n','n',NASHT,NASHT,NASHT,1.0D0,TDM,NASHT,
     &             TDMT,NASHT,0.0D0,Umat,NASHT)
      write (FILENAME,fmt='(a,a,a)')
     &"Dpart.",trim(adjustl(STATENAME)),Spin(I_NTO)
      LU=ISFREEUNIT(LU)
      CALL Molcas_Open(LU,FILENAME)
      DO I=1,NASHT
        write (LU,'(10(1X,ES11.4E2))') (Umat(NASHT*(I-1)+J),J=1,NASHT)
      END DO
      close (LU)
       CALL DSYEV_('V','U',NASHT,Vmat,NASHT,Veig,WGRONK,-1,INFO)
       NScrq=INT(WGRONK(1))
       If(Nscrq.eq.0) Then
        Nscrq=MAX(NDge,100)
        if(DoTest) Then
         write(6,*)'Size of scratch space is increased to max(NDge,100)'
        end if
       End If
       CALL mma_allocate (Scrq,NScrq,Label='Scrq')
C       Diagonalizing matrices
       CALL DSYEV_('V','U',NASHT,Umat,NASHT,Ueig,Scrq,NScrq,INFO)
       CALL DSYEV_('V','U',NASHT,Vmat,NASHT,Veig,Scrq,NScrq,INFO)
       CALL mma_deallocate (Scrq)
C     Printing some matrices
      If (DoTest) Then
       write (FILENAME,fmt='(a,a,a)')
     & "EigVecHole.",trim(adjustl(STATENAME)),Spin(I_NTO)
      LU=ISFREEUNIT(LU)
      CALL Molcas_Open(LU,FILENAME)
       DO I=1,NASHT
         write (LU,'(5(1X,ES11.4E2))') (Vmat(NASHT*(I-1)+J),J=1,NASHT)
       END DO
       close (LU)
       write (FILENAME,fmt='(a,a,a)')
     & "EigVecPart.",trim(adjustl(STATENAME)),Spin(I_NTO)
      LU=ISFREEUNIT(LU)
      CALL Molcas_Open(LU,FILENAME)
       DO I=1,NASHT
         write (LU,'(5(1X,ES11.4E2))') (Umat(NASHT*(I-1)+J),J=1,NASHT)
       END DO
       close (LU)
      End IF
C     End of Diagonlazing the mataces

C     Constructing hole and particle orbitals

      CALL DGEMM_('t','n',NASHT,NSupBas,NASHT,1.0D0,Umat,
     &      NASHT,SupCMO1,NASHT,0.0D0,ONTO,NASHT)
      CALL DGEMM_('t','n',NASHT,NSupBas,NASHT,1.0D0,Vmat,
     &      NASHT,SupCMO2,NASHT,0.0D0,UNTO,NASHT)

      If (DoTest) Then
       write(6,*)'printing Particle NTO in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(ONTO(I+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do

       write(6,*)'printing Hole     NTO in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(UNTO(I+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do
      End If

C     Printing NTOs
      CALL mma_allocate (Symto,NASHT,Label='Symto')
      CALL mma_allocate (Indto,NASHT,Label='Indto')
      CALL mma_allocate (Symfr,NASHT,Label='Symfr')
      CALL mma_allocate (Indfr,NASHT,Label='Indfr')
      NTOType='PART'
      CALL NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,ONTO,NTOType,
     &STATENAME,Ueig,UsetoReal,RealtoUse,Spin(I_NTO),Symto,Indto,
     &SumEigVal)
      NTOType='HOLE'
      CALL NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,UNTO,NTOType,
     &STATENAME,Veig,UsetoReal,RealtoUse,Spin(I_NTO),Symfr,Indfr,
     &SumEigVal)
C     End of Printing NTOs

      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         lIrrep(iSym) = adjustr(lIrrep(iSym))
      End Do

C     Putting particle-hole pairs in the output
        write(6,*)
      IF(I_NTO.eq.1) THEN
      write(6,'(10X,a)')
     &'NATURAL TRANSITION ORBTIAL INFORMATION FOR ALPHA SPIN'
      ELSE
      write(6,'(10X,a)')
     &'NATURAL TRANSITION ORBTIAL INFORMATION FOR BETA  SPIN'
      END IF
      WRITE(6,'(6X,A)') repeat('=',100)
      write(6,'(10X,5A18)')'EXCITATION','EIGENVALUE',
     &'EXCITATION','HOLE NTO','PARTICLE NTO'
      write(6,'(10X,A18,18X,3A18)')'AMPLITUDE','CONTRIBUTION(%)',
     &'SYMMETRY INDEX','SYMMETRY INDEX'
      WRITE(6,'(6X,A)') repeat('-',100)
      Do IOrb=NASHT,1,-1
       IF(Ueig(IOrb).lt.PrintThres)  EXIT
       write(6,'(10X,2(10X,F8.5),10X,F8.2,2(A9,I9))')
     & SQRT(Ueig(IOrb)),Ueig(IOrb),
     &                  Ueig(IOrb)/SumEigVal*1.0D2,
     & lIrrep(Symfr(IOrb)),Indfr(IOrb),
     & lIrrep(Symto(IOrb)),Indto(IOrb)
      End Do

      WRITE(6,'(6X,A)') repeat('-',100)
       write(6,'(6X,A,F8.5)')'SUM OF EIGENVALUES',SumEigVal
      WRITE(6,'(6X,A)') repeat('=',100)


      CALL mma_deallocate (Symto)
      CALL mma_deallocate (Indto)
      CALL mma_deallocate (Symfr)
      CALL mma_deallocate (Indfr)
      End DO
! End of loop over N_NTO (I_NTO=1 for alpha and 2 for beta)

       write(6,*)
       WRITE(6,'(6X,A)') repeat('*',100)
       WRITE(6,'(6X,A,33X,A34,31X,A)')
     &     '*','END OF NATURAL TRANSITION ORBITALS','*'
       WRITE(6,'(6X,A,33X,A14,I2,A12,I2 ,35X,A )')
     &'*','BETWEEN STATE ',JSTATE,' AND STATE ',ISTATE,'*'
      WRITE(6,'(6X,A)') repeat('*',100)

      CALL mma_deallocate(VMat)
      CALL mma_deallocate(UMat)
      CALL mma_deallocate(VEig)
      CALL mma_deallocate(UEig)
      CALL mma_deallocate(TDMT)
      CALL mma_deallocate(TDM)
      CALL mma_deallocate(CMO1)
      CALL mma_deallocate(CMO2)

      CALL mma_deallocate(SUPCMO2)
      CALL mma_deallocate(SUPCMO1)
      CALL mma_deallocate(UNTO)
      CALL mma_deallocate(ONTO)

      END SUBROUTINE NTOCalc



      SUBROUTINE  NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,NTO,
     &NTOType,STATENAME,EigVal,UsetoReal,RealtoUse,Spin,Sym,Ind,
     &SumEigVal)
      use Symmetry_Info, only: nSym=>nIrrep
      use Constants, only: Zero
      use rassi_data, only: NBST,NASHT,NASH,NBASF,NISH,NSSH

      Implicit None

C     input variables
      INTEGER NUseSym
      REAL*8 :: NTO(*)
      INTEGER,DIMENSION(NSym) :: NUseBF,NUsedBF,UsetoReal,RealtoUse
      CHARACTER (len=8) NTOType
      CHARACTER (len=1) Spin
      CHARACTER (len=5)  STATENAME
      Real*8 :: EigVal(*)
      Integer :: Sym(*), Ind(*)

      CHARACTER (len=128) FILENAME, molden_name
C     Loop control
      INTEGER I, J, ICount
C     Variables needed for judging the symmetry of a NTO
      INTEGER INTO,IUseSym,NNTO,ISym,IOrb
      REAL*8,DIMENSION(NUseSym) :: SquareSum
      REAL*8,DIMENSION(NBST) :: EigValArray
C     SquareSum=Sum over square of coefficients for certain symmetry
      INTEGER,DIMENSION(NUseSym) :: NOrbinSym
C     Total number of orbitals in IUseSym
      INTEGER,DIMENSION(NUseSym,NASHT) :: OrbSymIndex
C     OrbSymIndex gives the original orbital index for a orbital in iusesym
      REAL*8 Threshold,SumEigVal
C     If SquareSum(IUseSym) > Threshold, then print the coefficients in IUseSym symmetry
C     If there are more than one symmetry with SquareSum(IUseSym) > Threshold,
C     then give a warning message and print the one with the largest SquareSum
      INTEGER NPCMO,IPCMO
      Real*8,DIMENSION(:),allocatable::PCMO
C     Printing control
C
      INTEGER iPrintSym,OrbNum,IOrbinSym
      INTEGER LU
      Real*8,DIMENSION(2) :: vDum
      INTEGER,DIMENSION(7,8) :: v2Dum
      CHARACTER(len=72)Note
      Integer, External:: ISFREEUNIT

      Threshold=0.0D-10

      Do IUseSym=1,NUseSym
       NOrbinSym(IUseSym)=0
      End DO

      NPCMO=0
      Do ISym=1,NSym
       NPCMO=NPCMO+NBASF(ISym)**2
      End Do
      allocate(PCMO(NPCMO))


      Do INTO=1,NASHT
       iPrintSym=0
       Do IUseSym=1,NUseSym
        SquareSum(IUseSym)=Zero
        Do ICount=1,NUseBF(IUseSym)
         I=INTO
         J=ICount+NUsedBF(IUseSym)
         SquareSum(IUseSym)=SquareSum(IUseSym)
     &   +NTO(I+(J-1)*NASHT)**2
        End DO
        If (SquareSum(IUseSym).gt.Threshold) THEN
         If (iPrintSym.eq.0) Then
          iPrintSym=IUseSym
         Else
          write(6,'(a,a)')'There are at least two symmetries that have',
     &    ' a sum of coefficient**2 larger than the threshold'
       write(6,'(5A10)')'Threshold','Sum1','Sum2','Sym1','Sym2'
       write(6,'(3ES10.3E2,2I10)')Threshold,SquareSum(iPrintSym),
     & SquareSum(IUseSym),iPrintSym,IUseSym
          If (SquareSum(iPrintSym).lt.SquareSum(IUseSym)) THEN
           iPrintSym=IUseSym
          End If
         End IF
        End If
       End Do
       If(iPrintsym.eq.0) Then
        write(6,'(a,I2,a,a,a)') 'the symmetry of orbital ',INTO,
     & ' is not found. How is this possible?',
     & ' Change the value of DoTest in ntocalc.f to true ',
     & ' to print out the intermediate values'
       End If
        NOrbinSym(IPrintSym)=NOrbinSym(IPrintSym)+1
        OrbSymIndex(IPrintSym,NOrbinSym(IPrintSym))=INTO
      Sym(INTO)=UsetoReal(IPrintSym)
      Ind(INTO)=NOrbinSym(IPrintSym)+NISH(UsetoReal(IPrintSym))
      End Do

C     generating file in a similar way to other orbital files
      IOrb=0
      IPCMO=0
      SumEigVal=Zero
      Do ISym=1,NSym
       IUseSym=RealtoUse(ISym)
       If(IUseSym.ne.0) Then
C      If there are active orbitals in this symmetry
       NNTO=NOrbinSym(IUseSym)
C       write inactive part
       Do OrbNum=1,NISH(ISym)
        IOrb=IOrb+1
        EigValArray(IOrb)=Zero
       END DO
C Recording Printed NTO (PCMO)
       Do I=1,NISH(ISym)*NBASF(ISYM)
        IPCMO=IPCMO+1
        PCMO(IPCMO)=Zero
       End Do
C Recording Printed NTO (PCMO)
C       write active part
       Do IOrbinSym=1,NNTO
        OrbNum=IOrbinSym+NISH(ISym)
        IOrb=IOrb+1
        I=OrbSymIndex(IUseSym,IOrbinSym)
        EigValArray(IOrb)=EigVal(I)
        SumEigVal=SumEigVal+EigVal(I)
C Recording Printed NTO (PCMO)
        Do ICount=1,NUseBF(IUSeSym)
         IPCMO=IPCMO+1
         J=ICount+NUsedBF(IUseSym)
         PCMO(IPCMO)=NTO(I+(J-1)*NASHT)
        End Do
C Recording Printed NTO (PCMO)
       End Do
C       write virtual part
       Do OrBNum=NISH(ISym)+NASH(ISym)+1,NBASF(ISym)
        IOrb=IOrb+1
        EigValArray(IOrb)=Zero
       END DO
C Recording Printed NTO (PCMO)
       Do I=1,NSSH(ISym)*NBASF(ISYM)
        IPCMO=IPCMO+1
        PCMO(IPCMO)=Zero
       End Do
C Recording Printed NTO (PCMO)
       Else
C      If there is no active orbitals in this symmetry
       Do OrbNum=1,NBASF(ISym)
        IOrb=IOrb+1
        EigValArray(IOrb)=Zero
       END DO
C Recording Printed NTO (PCMO)
       Do I=1,NBASF(ISYM)**2
        IPCMO=IPCMO+1
        PCMO(IPCMO)=Zero
       End Do
C Recording Printed NTO (PCMO)
       End If
      End Do


      Do I=1,NBST
       EigValArray(I)=EigValArray(I)/SumEigVal
      End Do

      LU=50
      LU=ISFREEUNIT(LU)
      Note='*  Natural Transition Orbitals'
      WRITE(FILENAME,'(6(a))')
     & 'NTORB.SF.',trim(adjustl(STATENAME)),'.',Spin,'.',NTOType
      WRITE(molden_name,'(6(a))')
     & 'MD_NTO.SF.',trim(adjustl(STATENAME)),'.',Spin,'.',NTOType
      CALL WRVEC_(FILENAME,LU,'CO',0,NSYM,NBASF,NBASF,PCMO,vDum,
     & EigValArray,vDum,vDum,vDum,v2Dum,Note,0)
      CALL Molden_interface(0,trim(FILENAME),trim(molden_name))

      deallocate(PCMO)

      END SUBROUTINE  NTOSymAnalysis



