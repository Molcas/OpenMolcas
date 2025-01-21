Subroutine BP_Driver()
use Integral_interfaces, only: int_wrout, Int_Postprocess
use k2_arrays, only: DeDe
use definitions, only: wp
use stdalloc, only: mma_allocate, mma_deallocate
Implicit None
real(kind=wp) :: ThrAO
procedure(int_wrout) :: Integral_WrOut2_BP

! use the Breit option computing 1/r^3 integralas but convert to
! conventional 1/r integrals
Call set_Breit(1)

!     Prepare code to handle 2-particle densities
Call PrepP()

!     port P to module, to be developed....

!setting up post-processing routines

Int_PostProcess => Integral_WrOut2_BP
!     to be developed.
!     Modify Integral_Wrout2 to Integral_wrout_bp
!     Modify plf2 to plt2_BP
!     Modify IndSft2 to IndSft2_BP

call mma_allocate(DeDe,[-1,-1],label='DeDe') ! Dummy allocation

!     Compute the BP integrals and contract with P
ThrAO=1.0D-16
Call Drv2el(ThrAO)

!     Reset the environment
Call CloseP()
Call set_Breit(0)
Int_PostProcess => Null()
call mma_deallocate(DeDe)

End Subroutine BP_Driver

