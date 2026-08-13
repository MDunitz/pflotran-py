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
  PetscInt, parameter :: AWINHIBITMETHYL_ONE_MINUS_AW_INHIBITION = 3

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibitmethyl_type

    PetscReal :: aw_threshold
    PetscInt :: inhibition_type
    PetscReal :: fixed_water_activity

    ! Network-matching Monod kinetics for methylotrophic methanogenesis:
    !   CH3OH + H2(aq) -> CH4(aq) + H2O
    PetscReal :: rate_constant
    PetscReal :: half_saturation_ch3oh
    PetscReal :: half_saturation_h2
    PetscReal :: threshold_ch3oh
    PetscReal :: threshold_h2
    PetscReal :: o2_inhibition
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    PetscInt :: i_ch3oh
    PetscInt :: i_h2
    PetscInt :: i_ch4
    PetscInt :: i_o2

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

  implicit none

  class(reaction_sandbox_awinhibitmethyl_type), pointer :: AWInhibitMethylCreate

  allocate(AWInhibitMethylCreate)

  AWInhibitMethylCreate%aw_threshold = 0.95d0
  AWInhibitMethylCreate%inhibition_type = AWINHIBITMETHYL_SMOOTHSTEP_INHIBITION
  AWInhibitMethylCreate%fixed_water_activity = UNINITIALIZED_DOUBLE

  AWInhibitMethylCreate%rate_constant = UNINITIALIZED_DOUBLE
  AWInhibitMethylCreate%half_saturation_ch3oh = 1.0d-1
  AWInhibitMethylCreate%half_saturation_h2 = 1.0d-1
  AWInhibitMethylCreate%threshold_ch3oh = 1.1d-15
  AWInhibitMethylCreate%threshold_h2 = 1.1d-15
  AWInhibitMethylCreate%o2_inhibition = 1.0d-6
  AWInhibitMethylCreate%activation_energy = UNINITIALIZED_DOUBLE
  AWInhibitMethylCreate%reference_temperature = UNINITIALIZED_DOUBLE

  AWInhibitMethylCreate%i_ch3oh = UNINITIALIZED_INTEGER
  AWInhibitMethylCreate%i_h2 = UNINITIALIZED_INTEGER
  AWInhibitMethylCreate%i_ch4 = UNINITIALIZED_INTEGER
  AWInhibitMethylCreate%i_o2 = UNINITIALIZED_INTEGER

  nullify(AWInhibitMethylCreate%next)

end function AWInhibitMethylCreate

! ************************************************************************** !

subroutine AWInhibitMethylRead(this,input,option)

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
    call InputErrorMsg(input,option,'keyword',trim(error_string))
    call StringToUpper(word)

    select case(trim(word))
      case('WATER_ACTIVITY_THRESHOLD')
        call InputReadDouble(input,option,this%aw_threshold)
        call InputErrorMsg(input,option,'water_activity_threshold',error_string)
        if (this%aw_threshold < 0.d0 .or. this%aw_threshold > 1.d0) then
          option%io_buffer = 'WATER_ACTIVITY_THRESHOLD must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('FIXED_WATER_ACTIVITY')
        call InputReadDouble(input,option,this%fixed_water_activity)
        call InputErrorMsg(input,option,'fixed_water_activity',error_string)
        if (this%fixed_water_activity < 0.d0 .or. &
            this%fixed_water_activity > 1.d0) then
          option%io_buffer = 'FIXED_WATER_ACTIVITY must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant',error_string)
        call InputReadAndConvertUnits(input,this%rate_constant,'mol/L-sec',&
                        trim(error_string)//',rate_constant',option)

      case('HALF_SATURATION_CH3OH')
        call InputReadDouble(input,option,this%half_saturation_ch3oh)
        call InputErrorMsg(input,option,'half_saturation_ch3oh',error_string)

      case('HALF_SATURATION_H2')
        call InputReadDouble(input,option,this%half_saturation_h2)
        call InputErrorMsg(input,option,'half_saturation_h2',error_string)

      case('THRESHOLD_CH3OH')
        call InputReadDouble(input,option,this%threshold_ch3oh)
        call InputErrorMsg(input,option,'threshold_ch3oh',error_string)

      case('THRESHOLD_H2')
        call InputReadDouble(input,option,this%threshold_h2)
        call InputErrorMsg(input,option,'threshold_h2',error_string)

      case('O2_INHIBITION')
        call InputReadDouble(input,option,this%o2_inhibition)
        call InputErrorMsg(input,option,'o2_inhibition',error_string)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = AWINHIBITMETHYL_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = AWINHIBITMETHYL_SMOOTHSTEP_INHIBITION
          case('ONE_MINUS_AW')
            this%inhibition_type = AWINHIBITMETHYL_ONE_MINUS_AW_INHIBITION
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

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_awinhibitmethyl_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided for AWInhibitMethyl'
    call PrintErrMsg(option)
  endif

  word = 'CH3OH'
  this%i_ch3oh = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'H2(aq)'
  this%i_h2 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'CH4(aq)'
  this%i_ch4 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%i_o2 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)

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
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water, molality_to_molarity
  PetscReal :: water_activity, aw_inhibition, tempreal
  PetscReal :: rate_constant, reaction_rate
  PetscReal :: C_ch3oh, C_h2, C_o2
  PetscReal :: monod_ch3oh, monod_h2, inhib_o2

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  if (Initialized(this%fixed_water_activity)) then
    water_activity = this%fixed_water_activity
  else
    water_activity = exp(rt_auxvar%ln_act_h2o)
  endif

  rate_constant = this%rate_constant
  if (Initialized(this%activation_energy)) then
    rate_constant = rate_constant * Arrhenius(this%activation_energy, &
                                            global_auxvar%temp, &
                                            this%reference_temperature)
  endif

  C_ch3oh = rt_auxvar%pri_molal(this%i_ch3oh) * &
            rt_auxvar%pri_act_coef(this%i_ch3oh) * molality_to_molarity
  C_h2 = rt_auxvar%pri_molal(this%i_h2) * &
         rt_auxvar%pri_act_coef(this%i_h2) * molality_to_molarity
  C_o2 = rt_auxvar%pri_molal(this%i_o2) * &
         rt_auxvar%pri_act_coef(this%i_o2) * molality_to_molarity

  if (C_ch3oh < this%threshold_ch3oh .or. C_h2 < this%threshold_h2) return

  monod_ch3oh = C_ch3oh / (this%half_saturation_ch3oh + C_ch3oh)
  monod_h2 = C_h2 / (this%half_saturation_h2 + C_h2)
  inhib_o2 = this%o2_inhibition / (this%o2_inhibition + C_o2)

  select case(this%inhibition_type)
    case(AWINHIBITMETHYL_SMOOTHSTEP_INHIBITION)
      ! Positive threshold => INHIBIT_BELOW: factor -> 1 as a_w rises (wet on).
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.20d0, aw_inhibition, tempreal)
    case(AWINHIBITMETHYL_ONE_MINUS_AW_INHIBITION)
      if (this%aw_threshold >= 1.d0) then
        aw_inhibition = 1.d0
      else if (water_activity <= this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = (water_activity - this%aw_threshold) / &
                        (1.d0 - this%aw_threshold)
      endif
    case(AWINHIBITMETHYL_THRESHOLD_INHIBITION)
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  reaction_rate = rate_constant * monod_ch3oh * monod_h2 * inhib_o2 * &
                  aw_inhibition
  reaction_rate = reaction_rate * L_water

  Residual(this%i_ch3oh) = Residual(this%i_ch3oh) + reaction_rate
  Residual(this%i_h2) = Residual(this%i_h2) + reaction_rate
  Residual(this%i_ch4) = Residual(this%i_ch4) - reaction_rate

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

  implicit none
  class(reaction_sandbox_awinhibitmethyl_type) :: this

end subroutine AWInhibitMethylDestroy

end module Reaction_Sandbox_AWInhibitMethyl_class
