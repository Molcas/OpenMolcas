      Subroutine BP_Driver()
      use Integral_interfaces, only: Int_Postprocess
      Implicit None

! use the Breit option computing 1/r^3 integralas but convert to
! conventional 1/r integrals
      Call set_Breit(1)

!     Prepare code to handle 2-particle densities
      Call PrepP()

!     port P to module, to be developed....

!setting up post-processing routines

      Int_PostProcess => Integral_WrOut_BP 
! to be developed.
!       Integral_WrOut_BP calls PGet0!
!       Modify Integral_Wrout2 to Integral_wrout_bp
!       Modify plf2 to plt2_BP
!       Modify IndSft2 to IndSft2_BP

!     Compute the BP integrals and contract with P
      ThrAO=1.0D-16
      Call Drv2el(ThrAO)
    
!     Reset the environment
      Call CloseP()
      Call set_Breit(0)
      End Subroutine BP_Driver

