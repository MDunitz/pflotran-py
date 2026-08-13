module Reaction_Sandbox_AWInhibit_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  PetscInt, parameter :: AWINHIBIT_THRESHOLD_INHIBITION = 1
  PetscInt, parameter :: AWINHIBIT_SMOOTHSTEP_INHIBITION = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibit_type

    ! Water activity inhibition parameters
    PetscReal :: aw_threshold
    PetscInt :: inhibition_type

    ! Network-matching Monod kinetics for hydrogenotrophic methanogenesis:
    !   4 H2(aq) + HCO3- + H+ -> CH4(aq) + 3 H2O
    ! Rate = k * Monod(H2) * Monod(HCO3) * I_O2 * I_Fe * I_H * f(a_w)
    PetscReal :: rate_constant
    PetscReal :: half_saturation_h2
    PetscReal :: half_saturation_hco3
    PetscReal :: threshold_h2
    PetscReal :: threshold_hco3
    PetscReal :: o2_inhibition
    PetscReal :: fe_inhibition
    PetscReal :: h_inhibition
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    PetscInt :: i_h2
    PetscInt :: i_hco3
    PetscInt :: i_h
    PetscInt :: i_ch4
    PetscInt :: i_o2
    PetscInt :: i_fe

  contains
    procedure, public :: ReadInput  => AWInhibitRead
    procedure, public :: Setup      => AWInhibitSetup
    procedure, public :: Evaluate   => AWInhibitEvaluate
    procedure, public :: Destroy    => AWInhibitDestroy
  end type reaction_sandbox_awinhibit_type

  public :: AWInhibitCreate

contains

! ************************************************************************** !

function AWInhibitCreate()

  implicit none

  class(reaction_sandbox_awinhibit_type), pointer :: AWInhibitCreate

  allocate(AWInhibitCreate)

  AWInhibitCreate%aw_threshold = 0.95d0
  AWInhibitCreate%inhibition_type = AWINHIBIT_SMOOTHSTEP_INHIBITION

  AWInhibitCreate%rate_constant = UNINITIALIZED_DOUBLE
  AWInhibitCreate%half_saturation_h2 = 1.0d-1
  AWInhibitCreate%half_saturation_hco3 = 1.0d-1
  AWInhibitCreate%threshold_h2 = 1.1d-15
  AWInhibitCreate%threshold_hco3 = 1.1d-15
  AWInhibitCreate%o2_inhibition = 1.0d-6
  AWInhibitCreate%fe_inhibition = 1.0d-9
  AWInhibitCreate%h_inhibition = 1.78d-6
  AWInhibitCreate%activation_energy = UNINITIALIZED_DOUBLE
  AWInhibitCreate%reference_temperature = UNINITIALIZED_DOUBLE

  AWInhibitCreate%i_h2 = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_hco3 = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_h = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_ch4 = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_o2 = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_fe = UNINITIALIZED_INTEGER

  nullify(AWInhibitCreate%next)

end function AWInhibitCreate

! ************************************************************************** !

subroutine AWInhibitRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_awinhibit_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  error_string = 'CHEMISTRY,RXN_SANDBOX,AWINHIBIT'

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
        ! mol/L-sec matches MICROBIAL_REACTION RATE_CONSTANT; residual is
        ! assembled as rate * L_water -> mol/sec.
        call InputReadAndConvertUnits(input,this%rate_constant,'mol/L-sec',&
                        trim(error_string)//',rate_constant',option)

      case('HALF_SATURATION_H2')
        call InputReadDouble(input,option,this%half_saturation_h2)
        call InputErrorMsg(input,option,'half_saturation_h2',error_string)

      case('HALF_SATURATION_HCO3')
        call InputReadDouble(input,option,this%half_saturation_hco3)
        call InputErrorMsg(input,option,'half_saturation_hco3',error_string)

      case('THRESHOLD_H2')
        call InputReadDouble(input,option,this%threshold_h2)
        call InputErrorMsg(input,option,'threshold_h2',error_string)

      case('THRESHOLD_HCO3')
        call InputReadDouble(input,option,this%threshold_hco3)
        call InputErrorMsg(input,option,'threshold_hco3',error_string)

      case('O2_INHIBITION')
        call InputReadDouble(input,option,this%o2_inhibition)
        call InputErrorMsg(input,option,'o2_inhibition',error_string)

      case('FE_INHIBITION')
        call InputReadDouble(input,option,this%fe_inhibition)
        call InputErrorMsg(input,option,'fe_inhibition',error_string)

      case('H_INHIBITION')
        call InputReadDouble(input,option,this%h_inhibition)
        call InputErrorMsg(input,option,'h_inhibition',error_string)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = AWINHIBIT_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = AWINHIBIT_SMOOTHSTEP_INHIBITION
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

end subroutine AWInhibitRead

! ************************************************************************** !

subroutine AWInhibitSetup(this,reaction,option)

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_awinhibit_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided for AWInhibit reaction'
    call PrintErrMsg(option)
  endif

  word = 'H2(aq)'
  this%i_h2 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%i_hco3 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'H+'
  this%i_h = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'CH4(aq)'
  this%i_ch4 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%i_o2 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Fe+++'
  this%i_fe = ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  if (Initialized(this%activation_energy) .and. &
      UnInitialized(this%reference_temperature)) then
    option%io_buffer = 'A REFERENCE_TEMPERATURE must be provided when an &
      &ACTIVATION_ENERGY is defined in the AWInhibit Reaction Sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitSetup

! ************************************************************************** !

subroutine AWInhibitEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Dual-Monod hydrogenotrophic methanogenesis, inhibited by water activity.
  ! Matches the network MICROBIAL_REACTION rate law so this sandbox can replace
  ! that reaction rather than running as a dead parallel pathway.
  !

  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  class(reaction_sandbox_awinhibit_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
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
  PetscReal :: C_h2, C_hco3, C_h, C_o2, C_fe
  PetscReal :: monod_h2, monod_hco3, inhib_o2, inhib_fe, inhib_h

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  water_activity = exp(rt_auxvar%ln_act_h2o)

  rate_constant = this%rate_constant
  if (Initialized(this%activation_energy)) then
    rate_constant = rate_constant * Arrhenius(this%activation_energy, &
                                            global_auxvar%temp, &
                                            this%reference_temperature)
  endif

  C_h2 = rt_auxvar%pri_molal(this%i_h2) * &
         rt_auxvar%pri_act_coef(this%i_h2) * molality_to_molarity
  C_hco3 = rt_auxvar%pri_molal(this%i_hco3) * &
           rt_auxvar%pri_act_coef(this%i_hco3) * molality_to_molarity
  C_h = rt_auxvar%pri_molal(this%i_h) * &
        rt_auxvar%pri_act_coef(this%i_h) * molality_to_molarity
  C_o2 = rt_auxvar%pri_molal(this%i_o2) * &
         rt_auxvar%pri_act_coef(this%i_o2) * molality_to_molarity
  C_fe = rt_auxvar%pri_molal(this%i_fe) * &
         rt_auxvar%pri_act_coef(this%i_fe) * molality_to_molarity

  if (C_h2 < this%threshold_h2 .or. C_hco3 < this%threshold_hco3) then
    return
  endif

  monod_h2 = C_h2 / (this%half_saturation_h2 + C_h2)
  monod_hco3 = C_hco3 / (this%half_saturation_hco3 + C_hco3)
  ! Monod inhibition ABOVE threshold: Ki / (Ki + C)
  inhib_o2 = this%o2_inhibition / (this%o2_inhibition + C_o2)
  inhib_fe = this%fe_inhibition / (this%fe_inhibition + C_fe)
  inhib_h = this%h_inhibition / (this%h_inhibition + C_h)

  select case(this%inhibition_type)
    case(AWINHIBIT_SMOOTHSTEP_INHIBITION)
      ! Positive threshold => INHIBIT_BELOW polarity: factor -> 1 as a_w rises
      ! above the threshold (rate on when wet). Do not invert -- the old
      ! 1-factor flipped that and made salt *increase* methane.
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.05d0, aw_inhibition, tempreal)
    case(AWINHIBIT_THRESHOLD_INHIBITION)
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  reaction_rate = rate_constant * monod_h2 * monod_hco3 * &
                  inhib_o2 * inhib_fe * inhib_h * aw_inhibition
  reaction_rate = reaction_rate * L_water

  Residual(this%i_h2) = Residual(this%i_h2) + 4.d0 * reaction_rate
  Residual(this%i_hco3) = Residual(this%i_hco3) + reaction_rate
  Residual(this%i_h) = Residual(this%i_h) + reaction_rate
  Residual(this%i_ch4) = Residual(this%i_ch4) - reaction_rate

  if (compute_derivative) then
    option%io_buffer = 'REACTION_SANDBOX AWINHIBIT must be run with &
      &NUMERICAL_JACOBIAN listed in the NUMERICAL_METHODS TRANSPORT &
      &NEWTON_SOLVER block as analytical derivatives are not calculated &
      &in the sandbox evaluate routine.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitEvaluate

! ************************************************************************** !

subroutine AWInhibitDestroy(this)

  implicit none
  class(reaction_sandbox_awinhibit_type) :: this

end subroutine AWInhibitDestroy

end module Reaction_Sandbox_AWInhibit_class
