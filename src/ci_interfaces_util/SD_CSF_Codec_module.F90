module SD_CSF_Codec_Module

  use Definitions, only: wp, iwp

  implicit none
  private

  public :: SD_CSF_Codec
  public :: CSF_TO_SD
  public :: SD_TO_CSF

  integer(iwp), parameter :: CSF_TO_SD = 1
  integer(iwp), parameter :: SD_TO_CSF = 2

contains

  subroutine SD_CSF_Codec(CSF,SD,Direction,IRefSm)

    real(wp), intent(inout) :: CSF(:)
    real(wp), intent(inout) :: SD(:)

    integer(iwp), intent(in) :: Direction
    integer(iwp), intent(in) :: IRefSm

    select case(Direction)

    case(CSF_TO_SD)

      call CSDTVC(CSF,SD,1,IRefSm,0)

    case(SD_TO_CSF)

      call CSDTVC(CSF,SD,2,IRefSm,0)

    case default

      call Abend()

    end select

  end subroutine SD_CSF_Codec

end module SD_CSF_Codec_Module
