!> @file reaction_sandbox_aw_inhibit.F90
!! Reaction Sandbox that inhibits a reaction below water activity 0.50
!! Implements: 4 H2(aq) + HCO3- + H+ -> CH4(aq) + 3 H2O

module Reaction_Sandbox_AWInhibit_class

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use Option_module
  use Input_Aux_module
  use String_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_module, only: material_auxvar_type
  use Reactive_Transport_Aux_module

  implicit none
  private

  public :: AWInhibitCreate

  type, extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibit_type
     integer :: id_H2, id_HCO3, id_Hplus, id_CH4, id_H2O
     real(8) :: kmax
  contains
     procedure :: ReadInput => AWInhibitReadInput
     procedure :: Setup => AWInhibitSetup
     procedure :: Evaluate => AWInhibitEvaluate
     procedure :: Destroy => AWInhibitDestroy
  end type reaction_sandbox_awinhibit_type

contains

function AWInhibitCreate() result(this)
  type(reaction_sandbox_awinhibit_type), pointer :: this
  allocate(this)
  this%kmax = 1.0d-6
end function AWInhibitCreate

subroutine AWInhibitReadInput(this, input, option)
  class(reaction_sandbox_awinhibit_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  call InputReadUntilEnd(input, option)
end subroutine AWInhibitReadInput

subroutine AWInhibitSetup(this, reaction, option)
  class(reaction_sandbox_awinhibit_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%id_H2    = ReactionGetSpeciesId(reaction, 'H2(aq)', option)
  this%id_HCO3  = ReactionGetSpeciesId(reaction, 'HCO3-', option)
  this%id_Hplus = ReactionGetSpeciesId(reaction, 'H+', option)
  this%id_CH4   = ReactionGetSpeciesId(reaction, 'CH4(aq)', option)
  this%id_H2O   = ReactionGetSpeciesId(reaction, 'H2O', option)
end subroutine AWInhibitSetup

subroutine AWInhibitEvaluate(this, Residual, Jacobian, compute_derivative, &
                              rt_auxvar, global_auxvar, material_auxvar, &
                              reaction, option)

  class(reaction_sandbox_awinhibit_type) :: this
  real(8), dimension(:) :: Residual
  real(8), dimension(:,:) :: Jacobian
  logical :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  real(8) :: rate, aw
  real(8) :: conc_H2, conc_HCO3, conc_Hplus

  conc_H2    = rt_auxvar%aqueous_concentration(this%id_H2)
  conc_HCO3  = rt_auxvar%aqueous_concentration(this%id_HCO3)
  conc_Hplus = rt_auxvar%aqueous_concentration(this%id_Hplus)
  aw         = global_auxvar%water_activity_coefficient

  if (aw <= 0.5d0) then
    rate = 0.0d0
  else
    rate = this%kmax * conc_H2**4 * conc_HCO3 * conc_Hplus
  end if

  Residual(this%id_H2)    = Residual(this%id_H2)    - 4.0d0 * rate
  Residual(this%id_HCO3)  = Residual(this%id_HCO3)  - 1.0d0 * rate
  Residual(this%id_Hplus) = Residual(this%id_Hplus) - 1.0d0 * rate
  Residual(this%id_CH4)   = Residual(this%id_CH4)   + 1.0d0 * rate
  Residual(this%id_H2O)   = Residual(this%id_H2O)   + 3.0d0 * rate

end subroutine AWInhibitEvaluate

subroutine AWInhibitDestroy(this)
  class(reaction_sandbox_awinhibit_type) :: this
  ! No dynamic memory to deallocate
end subroutine AWInhibitDestroy

end module Reaction_Sandbox_AWInhibit_class
