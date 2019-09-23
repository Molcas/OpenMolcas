C       ***************************************************
C                          Calculating NTO
C       ****************************************************
C        Reference: J. Chem. Phys., 2014, 141, 024106
C
C                                         -------Jie Bao 
C                  in Depart. of Chemistry, University of Minnesota, USA
C                                             2018/09/22
C       A document is available to explain the code on Minnesota
C       OpenMolcas Website:
C       https://comp.chem.umn.edu/openmolcas


      SUBROUTINE   NTOCalc(ISTATE,JSTATE,LTRAD,LTRASD,ISpin)
C Include every head file included in the higher level code, namely
C gtdmctl.f
#include "rasdim.fh"
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "rassiwfn.fh"

      Integer LTRAD, LTRASD,ISpin
      Character,DIMENSION(2) :: Spin
      INTEGER Iprint,Jprint,I,J,isym                        ! Printing or looping contral
      INTEGER IUseSym,NUseSym,NSupBas,icactorb              ! CMO Symmetry Contral
      INTEGER,DIMENSION(NBST)  :: OrbUsedSym                ! CMO Symmetry Contral
      INTEGER,DIMENSION(NASHT+NISHT) :: OrbAct              ! CMO Symmetry Contral
      INTEGER,DIMENSION(NSym) :: NUsedBF,NUseBF,UsetoReal,RealtoUse
! Nr. of basis functions used prior to this symmetry (NUsedBF) and used in this symmetry (NUseBF) NSym >= NusedSym
      INTEGER   IOrb
      !IOrb is the index  of orbitals.
      INTEGER LSUPCMO1,LSUPCMO2,NSUPCMO
      INTEGER NDge,LNTOUmat,LNTOVmat,LNTOVeig
      INTEGER LTDM,LTDMT,LScrq,NScrq
      DIMENSION WGRONK(2)

      INTEGER LONTO, LUNTO,N_NTO,INFO, LNTOUeig,I_NTO
      INTEGER LSymfr,LIndfr,LSymto,LIndto
      Double Precision Zero,Two,SumEigVal

C      DIMENSION OrbArray(NCMO),EigVArray(NASHT),TDMArray(NASHT**2)
C     re-organizing orbitals 
C     This is to convert active MO sets in any symmetry into a C1 symmetry
      INTEGER NAISHT
      INTEGER, DIMENSION(NISHT+NASHT) :: OrbBas,OrbSym
      !OrbBas() is the number of basis function for IOrb
      !OrbSym() is the index of symmetry/irrep  for IOrb
      CHARACTER (len=14) FILENAME
      CHARACTER (len=8)  NTOType
      CHARACTER (len=2)  STATENAME
      Character*3 lIrrep(8)
      Logical DOTEST

      COMMON SumEigVal
      DoTest=.false.
      Zero=0.0D0
      Two=2.0D0
      if(dotest)then
C       write (6,*) 'LTRad ',LTRad,WORK(LTRAD)
C       write(6,*) 'Transition density matrix '
C       Do IPrint=1,NASHT
C       write (6,'(10(2X,F10.7))') 
C     & (WORK(LTRAD+JPrint-1+NASHT*(IPrint-1)),JPrint=1,NASHT)
C       End Do
       write(6,*) 'LCMO1 '
       Do I=0,NCMO,5
       write(6,'(2X,5F10.6)') 
     & (WORK(LCMO1+I+IPrint-1),IPrint=1,MIN(5,NCMO-I))
       End Do
       write(6,*) 'LCMO2 '
       Do I=0,NCMO,5
       write(6,'(2X,5F10.6)') 
     & (WORK(LCMO2+I+IPrint-1),IPrint=1,MIN(5,NCMO-I))
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
      CALL GETMEM ('SupCMO1','Allo','Real',LSUPCMO1,NSUPCMO)
      CALL GETMEM ('SupCMO2','Allo','Real',LSUPCMO2,NSUPCMO)
      CALL GETMEM ('ONTO','Allo','Real',LONTO,NSUPCMO)
      CALL GETMEM ('UNTO','Allo','Real',LUNTO,NSUPCMO)
      CALL DCOPY_(NSUPCMO,Zero,0,WORK(LSUPCMO1),1)
      CALL DCOPY_(NSUPCMO,Zero,0,WORK(LSUPCMO2),1)
      icactorb=0
      I=0
      Do IOrb=1,NAISHT
       IF (OrbAct(IOrb).eq.1) THEN
         icactorb=icactorb+1
         Do IPrint=1,(OrbBas(IOrb))
          JPrint=IPrint+NUsedBF(OrbUsedSym(IOrb))
          J=I+IPrint-1
C          write(6,'(4X,5I4,2F10.6)') IOrb,icactorb,
C     &    NUsedBF(OrbUsedSym(IOrb)),I,J,WORK(LCMO1+J),WORK(LCMO2+J)
          WORK(LSUPCMO1+icactorb-1+(JPRINT-1)*NASHT)=WORK(LCMO1+J)
          WORK(LSUPCMO2+icactorb-1+(JPRINT-1)*NASHT)=WORK(LCMO2+J)
        End DO
       End IF
       I=I+OrbBas(IOrb)
      End DO
      If (DoTest) Then
      write (6,*) 'LSupCMO1=',LSupCMO1
      write (6,*) 'LSupCMO2=',LSupCMO2
       write(6,*)'printing CMO1 in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(WORK(LSUPCMO1+I-1+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do
     
 
       write(6,*)'printing CMO2 in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(WORK(LSUPCMO2+I-1+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do
      End If 
C     end of building up the super-CMO matrix
C     Start and initialize spaces
      write (STATENAME,'(I0)') ISTATE
      NDge=NASHT**2
      CALL GETMEM ('Umat','Allo','Real',LNTOUmat,NDge)
      CALL GETMEM ('Vmat','Allo','Real',LNTOVmat,NDge)
      CALL GETMEM ('Ueig','Allo','Real',LNTOUeig,NDge)
      CALL GETMEM ('Veig','Allo','Real',LNTOVeig,NDge)
      CALL GETMEM ('TDM' ,'Allo','Real',LTDM,NDge)
      CALL GETMEM ('TDMT','Allo','Real',LTDMT,NDge)
       write(6,*)
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A31,33X,A)')
     &     '*','NATURAL TRANSITION ORBITALS','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,38X,A14,I2,A12,I2 ,30X,A )')
     &'*','BETWEEN STATE ',JSTATE,' AND STATE ',ISTATE,'*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,100A1)') ('*',i=1,100)
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
      CALL DCOPY_(Ndge,Zero,0,WORK(LNTOUeig),1)
      CALL DCOPY_(Ndge,Zero,0,WORK(LNTOVeig),1)
      CALL DCOPY_(NSupCMO,Zero,0,WORK(LONTO),1)
      CALL DCOPY_(NSupCMO,Zero,0,WORK(LUNTO),1)
!      CALL DCOPY_(Ndge,WORK(LTRAD),1,WORK(LTDM),1)
      If (Spin(I_NTO).eq.'a') Then
C       WORK(LTDM-1+I)=WORK(LTRAD-1+I)
      Do I=1,Ndge
       WORK(LTDM-1+I)=(WORK(LTRAD-1+I)+WORK(LTRASD-1+I))/Two
      End DO
      else
      Do I=1,Ndge
       WORK(LTDM-1+I)=(WORK(LTRAD-1+I)-WORK(LTRASD-1+I))/Two
      End DO
      End IF
      DO I=1,NASHT
       DO J=1,NASHT
        WORK(LTDMT+(I-1)+NASHT*(J-1))=WORK(LTDM+(J-1)+NASHT*(I-1))
       END DO
      END DO
C     Print out transition density matrix
      If(Spin(I_NTO).eq.'a') Then
      write (FILENAME,fmt='(a,a)')
     &"TDM",trim(adjustl(STATENAME))
      open (unit=233,file=FILENAME)
      DO I=1,NASHT
        write (233,'(5(2X,E10.4E2))')
     &    (WORK(LTRAD-1+NASHT*(I-1)+J),J=1,NASHT)
      END DO
      write (233,*)
      write (233,*)
      write (233,*)
      DO I=1,NASHT
        write (233,'(5(2X,E10.4E2))')
     &    (WORK(LTRASD-1+NASHT*(I-1)+J),J=1,NASHT)
      END DO
      close (233)
      End If
C     Generalizing transpose of TDM, TDM_T

C     Calculating T_trans*T
      CALL DGEMM_('n','n',NASHT,NASHT,NASHT,1.0D0,WORK(LTDMT),NASHT,
     &             WORK(LTDM),NASHT,0.0D0,WORK(LNTOVmat),NASHT)
C     Writing Particle Matrix
      write (FILENAME,fmt='(a,a,a)')
     &"Dhole.",trim(adjustl(STATENAME)),Spin(I_NTO)
      open (unit=233,file=FILENAME)
      DO I=1,NASHT
        write (233,'(10(2X,E10.4E2))')
     &    (WORK(LNTOVmat-1+NASHT*(I-1)+J),J=1,NASHT)
      END DO
      close (233)
C     Calculating T*T_transpose
      CALL DGEMM_('n','n',NASHT,NASHT,NASHT,1.0D0,WORK(LTDM),NASHT,
     &             WORK(LTDMT),NASHT,0.0D0,WORK(LNTOUmat),NASHT)
      write (FILENAME,fmt='(a,a,a)')
     &"Dpart.",trim(adjustl(STATENAME)),Spin(I_NTO)
      open (unit=233,file=FILENAME)
      DO I=1,NASHT
        write (233,'(10(2X,E10.4E2))')
     &    (WORK(LNTOUmat-1+NASHT*(I-1)+J),J=1,NASHT)
      END DO
      close (233)
       CALL DSYEV_('V','U',NASHT,WORK(LNTOVmat),NASHT,WORK(LNTOVeig),
     &              WGRONK,-1,INFO)
       NScrq=INT(WGRONK(1))
       If(Nscrq.eq.0) Then
        Nscrq=MAX(NDge,100)
        if(DoTest) Then
         write(6,*)'Size of scratch space is increased to max(NDge,100)'
        end if
       End If
       CALL GETMEM ('Scrq','Allo','Real',LScrq,NScrq)
C       Diagonalizing matrices
       CALL DSYEV_('V','U',NASHT,WORK(LNTOUmat),NASHT,WORK(LNTOUeig),
     &              WORK(LScrq),NScrq,INFO)
       CALL DSYEV_('V','U',NASHT,WORK(LNTOVmat),NASHT,WORK(LNTOVeig),
     &              WORK(LScrq),NScrq,INFO)
       CALL GETMEM ('Scrq','Free','Real',LScrq,NScrq)
C     Printing some matrices
      If (DoTest) Then
       write (FILENAME,fmt='(a,a,a)')
     & "EigVecHole.",trim(adjustl(STATENAME)),Spin(I_NTO)
       open (unit=233,file=FILENAME)
       DO I=1,NASHT
         write (233,'(5(2X,E10.4E2))')
     &     (WORK(LNTOVmat-1+NASHT*(I-1)+J),J=1,NASHT)
       END DO
       close (233)
       write (FILENAME,fmt='(a,a,a)')
     & "EigVecPart.",trim(adjustl(STATENAME)),Spin(I_NTO)
       open (unit=233,file=FILENAME)
       DO I=1,NASHT
         write (233,'(5(2X,E10.4E2))')
     &     (WORK(LNTOUmat-1+NASHT*(I-1)+J),J=1,NASHT)
       END DO
       close (233)
      End IF
C     End of Diagonlazing the mataces

C     Constructing hole and particle orbitals

      CALL DGEMM_('t','n',NASHT,NSupBas,NASHT,1.0D0,WORK(LNTOUmat),
     &      NASHT,WORK(LSupCMO1),NASHT,0.0D0,WORK(LONTO),NASHT)
      CALL DGEMM_('t','n',NASHT,NSupBas,NASHT,1.0D0,WORK(LNTOVmat),
     &      NASHT,WORK(LSupCMO2),NASHT,0.0D0,WORK(LUNTO),NASHT)

      If (DoTest) Then
       write(6,*)'printing Particle NTO in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(WORK(LONTO+I-1+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do
      
       write(6,*)'printing Hole     NTO in a C1-like format'
       Do I=1,NASHT
        Do J=1,NSupBas,10
         write(6,'(4X,10F10.6)')(WORK(LUNTO+I-1+(JPrint-1)*NASHT),
     &   JPrint=J,MIN(J+9,NSupBas))
        End DO
       End Do
      End If

C     Printing NTOs
      CALL GETMEM ('PartNTOSyms','Allo','Inte',LSymto,NASHT)
      CALL GETMEM ('PartNTOIndx','Allo','Inte',LIndto,NASHT)
      CALL GETMEM ('PartNTOSyms','Allo','Inte',LSymfr,NASHT)
      CALL GETMEM ('PartNTOIndx','Allo','Inte',LIndfr,NASHT)
      NTOType='PartNTO'
      CALL NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,LONTO,NTOType,
     &STATENAME,LNTOUeig,UsetoReal,RealtoUse,Spin(I_NTO),LSymto,LIndto)
      NTOType='HoleNTO'
      CALL NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,LUNTO,NTOType,
     &STATENAME,LNTOVeig,UsetoReal,RealtoUse,Spin(I_NTO),LSymfr,LIndfr)
C     End of Printing NTOs     

      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         Call RightAd(lIrrep(iSym))
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
      WRITE(6,'(6X,100A1)') ('=',i=1,100)
      write(6,'(10X,5A18)')'EXCITATION','EIGENVALUE',
     &'EXCITATION','HOLE NTO','PARTICLE NTO'
      write(6,'(10X,A18,18X,3A18)')'AMPLITUDE','CONTRIBUTION(%)',
     &'SYMMETRY INDEX','SYMMETRY INDEX'
      WRITE(6,'(6X,100A1)') ('-',i=1,100)
      Do IOrb=NASHT,1,-1
       IF(WORK(LNTOUeig-1+IOrb).lt.Zero) THEN
        WORK(LNTOUeig-1+IOrb)=Zero
       End IF
       write(6,'(10X,2(10X,F8.5),10X,F8.2,2(A9,I9))')
     & SQRT(WORK(LNTOUeig-1+IOrb)),WORK(LNTOUeig-1+IOrb),
     &WORK(LNTOUeig-1+IOrb)/SumEigVal*1.0D2,
     & lIrrep(INT(WORK(LSymfr-1+IOrb))),INT(WORK(LIndfr-1+IOrb)),
     & lIrrep(INT(WORK(LSymto-1+IOrb))),INT(WORK(LIndto-1+IOrb))
      End Do
       
      WRITE(6,'(6X,100A1)') ('-',i=1,100)
       write(6,'(6X,A,F8.5)')'SUM OF EIGENVALUES',SumEigVal
      WRITE(6,'(6X,100A1)') ('=',i=1,100)
       
 
      CALL GETMEM ('PartNTOSyms','Free','Inte',LSymto,NASHT)
      CALL GETMEM ('PartNTOIndx','Free','Inte',LIndto,NASHT)
      CALL GETMEM ('PartNTOSyms','Free','Inte',LSymfr,NASHT)
      CALL GETMEM ('PartNTOIndx','Free','Inte',LIndfr,NASHT)
      End DO  ! End of loop over N_NTO (I_NTO=1 for alpha and 2 for beta)

       write(6,*)
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
       WRITE(6,'(6X,A,33X,A34,31X,A)')
     &     '*','END OF NATURAL TRANSITION ORBITALS','*'
       WRITE(6,'(6X,A,33X,A14,I2,A12,I2 ,35X,A )')
     &'*','BETWEEN STATE ',JSTATE,' AND STATE ',ISTATE,'*'
      WRITE(6,'(6X,100A1)') ('*',i=1,100)

      CALL GETMEM ('Umat','Free','Real',LNTOUmat,NDge)
      CALL GETMEM ('Vmat','Free','Real',LNTOVmat,NDge)
      CALL GETMEM ('Ueig','Free','Real',LNTOUeig,NDge)
      CALL GETMEM ('Veig','Free','Real',LNTOVeig,NDge)
      CALL GETMEM ('TDM' ,'Free','Real',LTDM,NDge)
      CALL GETMEM ('TDMT','Free','Real',LTDMT,NDge)

      CALL GETMEM ('SupCMO1','Free','Real',LSUPCMO1,NSUPCMO)
      CALL GETMEM ('SupCMO2','Free','Real',LSUPCMO2,NSUPCMO)
      CALL GETMEM ('ONTO','Free','Real',LONTO,NSUPCMO)
      CALL GETMEM ('UNTO','Free','Real',LUNTO,NSUPCMO)
      RETURN
      END 



      SUBROUTINE  NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,LNTO,
     &NTOType,STATENAME,LEigVal,UsetoReal,RealtoUse,Spin,LSym,LInd)
#include "rasdim.fh"
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "rassiwfn.fh"
     
C     input variables
      INTEGER NUseSym,LNTO,LEigVal
      INTEGER,DIMENSION(NSym) :: NUseBF,NUsedBF,UsetoReal,RealtoUse
      CHARACTER (len=8) NTOType
      CHARACTER (len=1) Spin
      CHARACTER (len=2)  STATENAME

C     Loop control
      INTEGER I, J, ICount, JCount

C     Variables needed for judging the symmetry of a NTO
      INTEGER INTO,IUseSym,NNTO,ISym,IOrb
      REAL*8,DIMENSION(NUseSym) :: SquareSum               
      REAL*8,DIMENSION(NBST) :: EigValArray 
C     SquareSum=Sum over square of coefficients for certain symmetry
      INTEGER,DIMENSION(NUseSym) :: NOrbinSym                
C     Total number of orbitals in IUseSym
      INTEGER,DIMENSION(NUseSym,NASHT) :: OrbSymIndex
C     OrbSymIndex gives the original orbital index for a orbital in iusesym
      REAL*8 Threshold,Zero,SumEigVal                      
C     If SquareSum(IUseSym) > Threshold, then print the coefficients in IUseSym symmetry
C     If there are more than one symmetry with SquareSum(IUseSym) > Threshold,
C     then give a warning message and print the one with the largest SquareSum
      COMMON SumEigVal


C     Printing control      
C     
      INTEGER iPrintSym,OrbNum,IOrbinSym,LSym,LInd
      CHARACTER (len=14) FILENAME

      Logical DoTest


      DoTest=.false.
      Threshold=0.0D-10
      Zero=0.0D0
c      Do IUseSym=1, NUseSym
c       write (6,'(2I5)') IUseSym,NUseBF(IUseSym)
c      End Do

      Do IUseSym=1,NUseSym
       NOrbinSym(IUseSym)=0
      End DO


      Do INTO=1,NASHT
C       write(6,*) 'INTO=',INTO
       iPrintSym=0
       Do IUseSym=1,NUseSym
C       write(6,*) 'IUseSym',IUseSym
        SquareSum(IUseSym)=Zero
        Do ICount=1,NUseBF(IUseSym)
C       write(6,*) 'ICount',ICount
         I=INTO
         J=ICount+NUsedBF(IUseSym)
         SquareSum(IUseSym)=SquareSum(IUseSym)
     &   +WORK(LNTO+I-1+(J-1)*NASHT)**2
C         write(6,*)'square sum=',SquareSum(IUseSym) 
        End DO
        If (SquareSum(IUseSym).gt.Threshold) THEN
         If (iPrintSym.eq.0) Then
          iPrintSym=IUseSym
         Else
          write(6,'(a,a)')'There are at least two symmetries that have',
     &    ' a sum of coefficient**2 larger than the threshold'
       write(6,'(5A10)')'Threshold','Sum1','Sum2','Sym1','Sym2'
       write(6,'(3E10.6E2,2I10)')Threshold,SquareSum(iPrintSym),IUseSym,
     &    SquareSum(iPrintSym),SquareSum(IUseSym)
          If (SquareSum(iPrintSym).lt.SquareSum(IUseSym)) THEN
           iPrintSym=IUseSym
          End If
         End IF
        End If
C        write (6,*) 'printsym=',printsym
       End Do
       If(iPrintsym.eq.0) Then
        write(6,'(a,I2,a,a,a)') 'the symmetry of orbital ',INTO, 
     & ' is not found. How is this possible?',
     & ' Change the value of DoTest in ntocalc.f to true ',
     & ' to print out the intermediate values'
       End If
        NOrbinSym(IPrintSym)=NOrbinSym(IPrintSym)+1
        OrbSymIndex(IPrintSym,NOrbinSym(IPrintSym))=INTO
      WORK(LSym-1+INTO)=UsetoReal(IPrintSym)
      WORK(LInd-1+INTO)=NOrbinSym(IPrintSym)+NISH(UsetoReal(IPrintSym))
      End Do


      WRITE(FILENAME,'(a,a,a,a)') 
     & trim(adjustl(STATENAME)),Spin,'.',NTOType
      OPEN  (unit=233,file=FILENAME)

C     generating file in a similar way to other orbital files
      write(233,'(a)') '#INPORB 2.2'
      write(233,'(a)') '#INFO'
      write(233,'(a)') '*  Natural Transition Orbitals'
      write(233,'(3i8)') 0,NSYM,0
      write(233,'(8i8)') (NBASF(I),I=1,NSYM)
      write(233,'(8i8)') (NBASF(I),I=1,NSYM)


      IOrb=0
      SumEigVal=Zero
      write(233,'(a)') '#ORB'
      Do ISym=1,NSym
       IUseSym=RealtoUse(ISym)
       If(IUseSym.ne.0) Then
C      If there are active orbitals in this symmetry
       NNTO=NOrbinSym(IUseSym)
C       write inactive part
       Do OrbNum=1,NISH(ISym)  
        IOrb=IOrb+1
        EigValArray(IOrb)=Zero
        write(233,'(A,2I5)') '* ORBITAL',ISym,OrbNum
        Do I=1,NBASF(ISym),5
        write(233,'(5E22.14E2)') (Zero,J=1,MIN(5,NBASF(ISym)+1-I)) 
        END DO
       END DO
C       write active part
       Do IOrbinSym=1,NNTO
        OrbNum=IOrbinSym+NISH(ISym)
        IOrb=IOrb+1
        I=OrbSymIndex(IUseSym,IOrbinSym)
        EigValArray(IOrb)=WORK(LEigVal-1+I)
        SumEigVal=SumEigVal+WORK(LEigVal-1+I)
        write(233,'(A,2I5)') '* ORBITAL',ISym,OrbNum
        Do ICount=1,NUseBF(IUseSym),5
         J=ICount+NUsedBF(IUseSym)
         write(233,'(5E22.14E2)')(WORK(LNTO+I-1+(JCount-1)*NASHT),
     &   JCount=J,MIN(J+4,NUseBF(IUseSym)+NUsedBF(IUseSym)))
        End DO
       End Do
C       write virtual part
       Do OrBNum=NISH(ISym)+NASH(ISym)+1,NBASF(ISym)
        IOrb=IOrb+1
        EigValArray(IOrb)=Zero
        write(233,'(A,2I5)') '* ORBITAL',ISym,OrbNum
        Do I=1,NBASF(ISym),5
        write(233,'(5E22.14E2)') (Zero,J=1,MIN(5,NBASF(ISym)+1-I)) 
        END DO
       END DO
       Else
C      If there is no active orbitals in this symmetry       
       Do OrbNum=1,NBASF(ISym)  
        IOrb=IOrb+1
        EigValArray(IOrb)=Zero
        write(233,'(A,2I5)') '* ORBITAL',ISym,OrbNum
        Do I=1,NBASF(ISym),5
        write(233,'(5E22.14E2)') (Zero,J=1,MIN(5,NBASF(ISym)+1-I)) 
        END DO
       END DO
       End If
      End Do
      
      write(233,'(a)') '#OCC'
      write(233,'(a)') '* EXCITATION CONTRIBUTION'
      IOrb=0
      DO ISym=1,NSym
       Do I=IOrb+1,IOrb+NBASF(ISym),5
      write(233,'(5E22.14E2)')
     & (EigValArray(I+J)/SumEigVal,J=0,MIN(4,NBASF(ISym)+IOrb-I))
       End Do
       IOrb=IOrb+NBASF(ISym)
      End Do

      write(233,'(a)') '#OCHR'
      write(233,'(a,a)') '* EXCITATION CONTRIBUTION'
     &, '(HUMAN-READABLE)'
      IOrb=0
      DO ISym=1,NSym
       Do I=IOrb+1,IOrb+NBASF(ISym),10
      write(233,'(10F8.4)')
     &(EigValArray(I+J)/SumEigVal,J=0,MIN(9,NBASF(ISym)+IOrb-I))
       End Do
       IOrb=IOrb+NBASF(ISym)
      End Do
      CLOSE (233)

      RETURN
      END



