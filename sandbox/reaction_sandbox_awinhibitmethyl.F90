module Reaction_Sandbox_AWInhibitMethyl_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  PetscInt, parameter :: AWINHIBITMETHYL_THRESHOLD_INHIBITION = 1
  PetscInt, parameter :: AWINHIBITMETHYL_SMOOTHSTEP_INHIBITION = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibitmethyl_type
    
    ! Water activity inhibition parameters
    PetscReal :: aw_threshold
    PetscInt :: inhibition_type

    ! Reaction parameters for methylotrophic methanogenesis: CH3OH + H2(aq) -> CH4(aq) + H2O
    PetscReal :: rate_constant
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    ! Species names and indices
    character(len=MAXWORDLENGTH) :: species_names(3)
    PetscInt :: nspecies
    PetscInt :: i_ch3oh
    PetscInt :: i_h2
    PetscInt :: i_ch4

  contains
    procedure, public :: ReadInput  => AWInhibitMethylRead
    procedure, public :: Setup      => AWInhibitMethylSetup
    procedure, public :: Evaluate   => AWInhibitMethylEvaluate
    procedure, public :: Destroy    => AWInhibitMethylDestroy
  end type reaction_sandbox_awinhibitmethyl_type

  public :: AWInhibitMethylCreate

contains

! ************************************************************************** !

function AWInhibitMethylCreate()
  !
  ! Allocates AWInhibitMethyl reaction object for water activity inhibited methylotrophic methanogenesis
  !

  implicit none

  class(reaction_sandbox_awinhibitmethyl_type), pointer :: AWInhibitMethylCreate

  allocate(AWInhibitMethylCreate)
  
  ! Default water activity threshold of 0.5
  AWInhibitMethylCreate%aw_threshold = 0.5d0
  AWInhibitMethylCreate%inhibition_type = AWINHIBITMETHYL_THRESHOLD_INHIBITION
  
  ! Reaction parameters
  AWInhibitMethylCreate%rate_constant = UNINITIALIZED_DOUBLE
  AWInhibitMethylCreate%activation_energy = UNINITIALIZED_DOUBLE
  AWInhibitMethylCreate%reference_temperature = UNINITIALIZED_DOUBLE
  
  ! Species indices
  AWInhibitMethylCreate%i_ch3oh = UNINITIALIZED_INTEGER
  AWInhibitMethylCreate%i_h2 = UNINITIALIZED_INTEGER
  AWInhibitMethylCreate%i_ch4 = UNINITIALIZED_INTEGER

  nullify(AWInhibitMethylCreate%next)

end function AWInhibitMethylCreate

! ************************************************************************** !

subroutine AWInhibitMethylRead(this,input,option)
  !
  ! Reads input deck for water activity inhibited methylotrophic methanogenesis parameters
  !

  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_awinhibitmethyl_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  error_string = 'AWINHIBITMETHYL'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       trim(error_string))
    call StringToUpper(word)

    select case(trim(word))
      case('WATER_ACTIVITY_THRESHOLD')
        call InputReadDouble(input,option,this%aw_threshold)
        call InputErrorMsg(input,option,'water_activity_threshold',error_string)
        if (this%aw_threshold < 0.d0 .or. this%aw_threshold > 1.d0) then
          option%io_buffer = 'WATER_ACTIVITY_THRESHOLD must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant',error_string)
        call InputReadAndConvertUnits(input,this%rate_constant,'mol/m^3-sec',&
                        trim(error_string)//',rate_constant',option)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = AWINHIBITMETHYL_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = AWINHIBITMETHYL_SMOOTHSTEP_INHIBITION
          case default
            error_string = trim(error_string) // ',INHIBITION_TYPE'
            call InputKeywordUnrecognized(input,word,error_string ,option)
        end select

      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,this%activation_energy)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%activation_energy,'j/mol',&
                          trim(error_string)//',activation energy',option)

      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,this%reference_temperature)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%reference_temperature,'C',&
                          trim(error_string)//',reference temperature',option)

      case default
        call InputKeywordUnrecognized(input,word,error_string ,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine AWInhibitMethylRead

! ************************************************************************** !

subroutine AWInhibitMethylSetup(this,reaction,option)
  !
  ! Sets up the water activity inhibited methylotrophic methanogenesis reaction
  !

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_awinhibitmethyl_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  ! Check that rate constant is provided
  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided for AWInhibitMethyl reaction'
    call PrintErrMsg(option)
  endif

  ! Get species indices for methylotrophic methanogenesis reaction: CH3OH + H2(aq) -> CH4(aq) + H2O
  word = 'CH3OH'
  this%i_ch3oh = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'H2(aq)'
  this%i_h2 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'CH4(aq)'
  this%i_ch4 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  if (Initialized(this%activation_energy) .and. &
      UnInitialized(this%reference_temperature)) then
    option%io_buffer = 'A REFERENCE_TEMPERATURE must be provided when an &
      &ACTIVATION_ENERGY is defined in the AWInhibitMethyl Reaction Sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitMethylSetup

! ************************************************************************** !

subroutine AWInhibitMethylEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Evaluates the water activity inhibited methylotrophic methanogenesis reaction
  ! Reaction: CH3OH + H2(aq) -> CH4(aq) + H2O
  ! Inhibited below water activity threshold
  !

  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  class(reaction_sandbox_awinhibitmethyl_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: molality_to_molarity
  PetscReal :: water_activity, aw_inhibition, tempreal
  PetscReal :: rate_constant
  PetscReal :: reaction_rate

  ! Concentrations
  PetscReal :: C_ch3oh, C_h2, C_ch4

  ! Calculate water volume in liters
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! m^3 -> L

  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3

  ! Get water activity from rt_auxvar
  water_activity = exp(rt_auxvar%ln_act_h2o)

  ! Apply temperature correction if specified
  rate_constant = this%rate_constant
  if (Initialized(this%activation_energy)) then
    rate_constant = rate_constant * Arrhenius(this%activation_energy, &
                                            global_auxvar%temp, &
                                            this%reference_temperature)
  endif

  ! Get aqueous concentrations (convert from molality to molarity)
  C_ch3oh = rt_auxvar%pri_molal(this%i_ch3oh) * &
            rt_auxvar%pri_act_coef(this%i_ch3oh) * molality_to_molarity
  C_h2 = rt_auxvar%pri_molal(this%i_h2) * &
         rt_auxvar%pri_act_coef(this%i_h2) * molality_to_molarity
  C_ch4 = rt_auxvar%pri_molal(this%i_ch4) * &
          rt_auxvar%pri_act_coef(this%i_ch4) * molality_to_molarity

  ! Calculate water activity inhibition
  select case(this%inhibition_type)
    case(AWINHIBITMETHYL_SMOOTHSTEP_INHIBITION)
      ! Smooth transition around threshold
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.05d0, aw_inhibition, tempreal)
      ! For smoothstep, we want inhibition when aw < threshold, so invert
      aw_inhibition = 1.d0 - aw_inhibition
    case(AWINHIBITMETHYL_THRESHOLD_INHIBITION)
      ! Sharp threshold
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  ! Calculate reaction rate: CH3OH + H2(aq) -> CH4(aq) + H2O
  ! Rate law assumes first order in each reactant
  reaction_rate = rate_constant * C_ch3oh * C_h2 * aw_inhibition
  
  ! Apply water activity inhibition
  reaction_rate = reaction_rate * L_water  ! Convert to mol/sec

  ! Update residuals (negative stoichiometry for reactants, positive for products)
  ! 1 CH3OH consumed
  Residual(this%i_ch3oh) = Residual(this%i_ch3oh) + reaction_rate
  
  ! 1 H2(aq) consumed
  Residual(this%i_h2) = Residual(this%i_h2) + reaction_rate
  
  ! 1 CH4(aq) produced
  Residual(this%i_ch4) = Residual(this%i_ch4) - reaction_rate

  ! Note: H2O is not tracked as it's the solvent

  if (compute_derivative) then
    option%io_buffer = 'REACTION_SANDBOX AWINHIBITMETHYL must be run with &
      &NUMERICAL_JACOBIAN listed in the NUMERICAL_METHODS TRANSPORT &
      &NEWTON_SOLVER block as analytical derivatives are not calculated &
      &in the sandbox evaluate routine.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitMethylEvaluate

! ************************************************************************** !

subroutine AWInhibitMethylDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !

  implicit none
  class(reaction_sandbox_awinhibitmethyl_type) :: this

  ! Nothing to deallocate in this simple version

end subroutine AWInhibitMethylDestroy

end module Reaction_Sandbox_AWInhibitMethyl_class