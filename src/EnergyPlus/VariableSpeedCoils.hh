// EnergyPlus, Copyright (c) 1996-2020, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without the U.S. Department of Energy's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef VariableSpeedCoils_hh_INCLUDED
#define VariableSpeedCoils_hh_INCLUDED

// ObjexxFCL Headers
#include <ObjexxFCL/Array1D.fwd.hh>
#include <ObjexxFCL/Optional.fwd.hh>

// EnergyPlus Headers
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/EnergyPlus.hh>

namespace EnergyPlus {

// Forward declarations
struct EnergyPlusData;

namespace VariableSpeedCoils {

    struct VariableSpeedCoilData // variable speed coil
    {
        int MaxSpedLevels = 10;
        // Members
        std::string Name;              // Name of the  Coil
        std::string VarSpeedCoilType;  // type of coil
        int NumOfSpeeds = 2;               // Number of speeds
        int NormSpedLevel = MaxSpedLevels;          // Nominal speed level
        Real64 RatedWaterVolFlowRate = EnergyPlus::DataSizing::AutoSize;  // Rated/Ref Water Volumetric Flow Rate [m3/s]
        Real64 RatedWaterMassFlowRate = EnergyPlus::DataSizing::AutoSize; // Rated/Ref Water Volumetric Flow Rate [m3/s]
        Real64 RatedAirVolFlowRate = EnergyPlus::DataSizing::AutoSize;    // Rated/Ref Air Volumetric Flow Rate [m3/s]
        Real64 RatedCapHeat = EnergyPlus::DataSizing::AutoSize;           // Rated/Ref Heating Capacity [W]
        Real64 RatedCapCoolTotal = EnergyPlus::DataSizing::AutoSize;      // Rated/Ref Total Cooling Capacity [W]
        Real64 MaxONOFFCyclesperHour = 0.0;  // Maximum ON/OFF cycles per hour for the compressor (cycles/hour)
        Real64 Twet_Rated = 0.0;             // Nominal time for condensate to begin leaving the coil's
        // condensate drain line (sec)
        Real64 Gamma_Rated = 0.0; // Initial moisture evaporation rate divided by steady-state
        // AC latent capacity (dimensionless)
        int HOTGASREHEATFLG = 0;            // whether to use hot gas reheat
        Real64 HPTimeConstant = 0.0;          // Heat pump time constant [s]
        int PLFFPLR = 0;                    // index of part load curve as a function of part load ratio
        std::string CoolHeatType;       // Type of WatertoAirHP ie. Heating or Cooling
        int VSCoilTypeOfNum = 0;            // type of component in plant
        bool SimFlag = false;                   // Heat Pump Simulation Flag
        Real64 DesignWaterMassFlowRate = 0.0; // design water mass flow rate [kg/s]
        Real64 DesignWaterVolFlowRate = 0.0;  // design water volumetric flow rate [m3/s]
        Real64 DesignAirMassFlowRate = 0.0;   // Design Air Mass Flow Rate [kg/s]
        Real64 DesignAirVolFlowRate = 0.0;    // Design Air Volumetric Flow Rate [m3/s]
        Real64 AirVolFlowRate = 0.0;          // Air Volumetric Flow Rate[m3/s], real time
        Real64 AirMassFlowRate = 0.0;         // Air Mass Flow Rate[kg/s], real time
        Real64 InletAirPressure = 0.0;        // air inlet pressure [pa]
        Real64 InletAirDBTemp = 0.0;          // Inlet Air Dry Bulb Temperature [C], real time
        Real64 InletAirHumRat = 0.0;          // Inlet Air Humidity Ratio [kg/kg], real time
        Real64 InletAirEnthalpy = 0.0;        // Inlet Air Enthalpy [J/kg], real time
        Real64 OutletAirDBTemp = 0.0;         // Outlet Air Dry Bulb Temperature [C], real time
        Real64 OutletAirHumRat = 0.0;         // Outlet Air Humidity Ratio [kg/kg], real time
        Real64 OutletAirEnthalpy = 0.0;       // Outlet Air Enthalpy [J/kg], real time
        Real64 WaterVolFlowRate = 0.0;        // Water Volumetric Flow Rate [m3/s], real time
        Real64 WaterMassFlowRate = 0.0;       // Water Mass Flow Rate [kg/s], real time
        Real64 InletWaterTemp = 0.0;          // Inlet Water Temperature [C]
        Real64 InletWaterEnthalpy = 0.0;      // Inlet Water Enthalpy [J/kg]
        Real64 OutletWaterTemp = 0.0;         // Outlet Water Temperature [C]
        Real64 OutletWaterEnthalpy = 0.0;     // Outlet Water Enthalpy [J/kg]
        Real64 Power = 0.0;                   // Power Consumption [W]
        Real64 QLoadTotal = 0.0;              // Load Side Total Heat Transfer Rate [W]
        Real64 QSensible = 0.0;               // Sensible Load Side Heat Transfer Rate [W]
        Real64 QLatent = 0.0;                 // Latent Load Side Heat Transfer Rate [W]
        Real64 QSource = 0.0;                 // Source Side Heat Transfer Rate [W]
        Real64 QWasteHeat = 0.0;              // Recoverable waste Heat Transfer Rate [W]
        Real64 Energy = 0.0;                  // Energy Consumption [J]
        Real64 EnergyLoadTotal = 0.0;         // Load Side Total Heat Transferred [J]
        Real64 EnergySensible = 0.0;          // Sensible Load Side Heat Transferred [J]
        Real64 EnergyLatent = 0.0;            // Latent Load Side Heat Transferred [J]
        Real64 EnergySource = 0.0;            // Source Side Heat Transferred [J]
        Real64 COP = 0.0;                     // Heat Pump Coefficient of Performance [-]
        Real64 RunFrac = 0.0;                 // Duty Factor
        Real64 PartLoadRatio = 0.0;           // Part Load Ratio
        Real64 RatedPowerHeat = 0.0;          // Rated/Ref Heating Power Consumption[W]
        Real64 RatedCOPHeat = 0.0;            // Rated/Ref Heating COP [W/W]
        Real64 RatedCapCoolSens = 0.0;        // Rated/Ref Sensible Cooling Capacity [W]
        Real64 RatedPowerCool = 0.0;          // Rated/Ref Cooling Power Consumption[W]
        Real64 RatedCOPCool = 0.0;            // Rated/Ref Cooling COP [W/W]
        int AirInletNodeNum = 0;            // Node Number of the Air Inlet
        int AirOutletNodeNum = 0;           // Node Number of the Air Outlet
        int WaterInletNodeNum = 0;          // Node Number of the Water Onlet
        int WaterOutletNodeNum = 0;         // Node Number of the Water Outlet
        int LoopNum = 0;                    // plant loop index for water side
        int LoopSide = 0;                   // plant loop side index
        int BranchNum = 0;                  // plant branch index
        int CompNum = 0;                    // plant component index
        // set by parent object and "pushed" to this structure in SetVSWSHPData subroutine
        bool FindCompanionUpStreamCoil = true; // Flag to get the companion coil in Init
        bool IsDXCoilInZone = false;            // true means dx coil is in zone instead of outside
        int CompanionCoolingCoilNum = 0;    // Heating coil companion cooling coil index
        int CompanionHeatingCoilNum = 0;    // Cooling coil companion heating coil index
        Real64 FanDelayTime = 0.0;            // Fan delay time, time delay for the HP's fan to
        // beginning for multispeed coil type
        int MSHPDesignSpecIndex = -1;              // index to UnitarySystemPerformance:Multispeed object
        Array1D_int MSErrIndex;               // index flag for num speeds/recurring messages
        Array1D<Real64> MSRatedPercentTotCap; // Percentage to the total cooling capacity for MS heat pump at the highest speed [dimensionless]
        Array1D<Real64> MSRatedTotCap;        // Rated cooling capacity for MS heat pump [W]
        Array1D<Real64> MSRatedSHR;           // Rated SHR for MS heat pump [dimensionless]
        Array1D<Real64> MSRatedCOP;           // Rated COP for MS heat pump [dimensionless]
        Array1D<Real64> MSRatedAirVolFlowPerRatedTotCap;
        // Rated Air volume flow rate per total capacity through unit at rated conditions [m^3/w]
        Array1D<Real64> MSRatedAirVolFlowRate;
        // Air volume flow rate through unit at rated conditions [m3/s]
        Array1D<Real64> MSRatedAirMassFlowRate;
        // Air mass flow rate through unit at rated conditions [kg/s]
        Array1D<Real64> MSRatedWaterVolFlowPerRatedTotCap;
        // Rated water volume flow rate per total  capacity through unit at rated conditions [m^3/w]
        Array1D<Real64> MSRatedWaterVolFlowRate;
        // Water volume flow rate through unit at rated conditions [m3/s]
        Array1D<Real64> MSRatedWaterMassFlowRate;
        // Water mass flow rate through unit at rated conditions [kg/s]
        Array1D<Real64> MSRatedCBF;
        // rated coil bypass factor
        Array1D<Real64> MSEffectiveAo;
        // effective heat transfer surface at each speed
        Array1D_int MSCCapFTemp;
        // index of total capacity modifier curve
        Array1D_int MSCCapAirFFlow;
        // index of total capacity modifier curve as a function of air flow
        Array1D_int MSCCapWaterFFlow;
        // index of total capacity modifier curve as a function of water flow
        Array1D_int MSEIRFTemp;
        // index of energy input ratio modifier curve as a function of temperature
        Array1D_int MSEIRAirFFlow;
        // index of energy input ratio modifier curve as a function of air flow fraction
        Array1D_int MSEIRWaterFFlow;
        // index of energy input ratio modifier curve as a function of water flow fraction
        Array1D_int MSWasteHeat;
        // index of waste heat as a function of temperature
        Array1D<Real64> MSWasteHeatFrac;
        // water heating coil pump power at various speeds
        Array1D<Real64> MSWHPumpPower;
        Array1D<Real64> MSWHPumpPowerPerRatedTotCap;
        // Waste heat fraction
        Real64 SpeedNumReport = 0.0;
        // speed number for output
        Real64 SpeedRatioReport = 0.0;
        // speed ratio for output between two neighboring speeds
        // End of multispeed water source coil input
        //----------------------------------------------------------------
        // added variables and arrays for variable speed air-source heat pump
        // defrosting
        int DefrostStrategy = 0;       // defrost strategy; 1=reverse-cycle, 2=resistive
        int DefrostControl = 0;        // defrost control; 1=timed, 2=on-demand
        int EIRFPLR = 0;               // index of energy input ratio vs part-load ratio curve
        int DefrostEIRFT = 0;          // index of defrost mode total cooling capacity for reverse cycle heat pump
        Real64 MinOATCompressor = 0.0;   // Minimum OAT for heat pump compressor operation
        Real64 OATempCompressorOn = 0.0; // The outdoor tempearture when the compressor is automatically turned back on,
        // if applicable, following automatic shut off. This field is used only for
        // HSPF calculation.
        Real64 MaxOATDefrost = 0.0;           // Maximum OAT for defrost operation
        Real64 DefrostTime = 0.0;             // Defrost time period in hours
        Real64 DefrostCapacity = 0.0;         // Resistive defrost to nominal capacity (at 21.11C/8.33C) ratio
        Real64 HPCompressorRuntime = 0.0;     // keep track of compressor runtime
        Real64 HPCompressorRuntimeLast = 0.0; // keep track of last time step compressor runtime (if simulation downshifts)
        Real64 TimeLeftToDefrost = 0.0;       // keep track of time left to defrost heat pump
        Real64 DefrostPower = 0.0;            // power used during defrost
        Real64 DefrostConsumption = 0.0;      // energy used during defrost
        // crankcase heater
        bool ReportCoolingCoilCrankcasePower = true; // logical determines if the cooling coil crankcase heater power is reported
        Real64 CrankcaseHeaterCapacity = 0.0;       // total crankcase heater capacity [W]
        Real64 CrankcaseHeaterPower = 0.0;          // report variable for average crankcase heater power [W]
        Real64 MaxOATCrankcaseHeater = 0.0;         // maximum OAT for crankcase heater operation [C]
        Real64 CrankcaseHeaterConsumption = 0.0;    // report variable for total crankcase heater energy consumption [J]
        // condenser evaporative precooling
        int CondenserInletNodeNum= 0;       // Node number of outdoor condenser
        int CondenserType = EnergyPlus::DataHVACGlobals::AirCooled;               // Type of condenser for DX cooling coil: AIR COOLED or EVAP COOLED
        bool ReportEvapCondVars = false;         // true if any performance mode includes an evap condenser
        Real64 EvapCondPumpElecNomPower = 0.0; // Nominal power input to the evap condenser water circulation pump [W]
        Real64 EvapCondPumpElecPower = 0.0;    // Average power consumed by the evap condenser water circulation pump over
        // the time step [W]
        Real64 EvapWaterConsumpRate = 0.0;        // Evap condenser water consumption rate [m3/s]
        Real64 EvapCondPumpElecConsumption = 0.0; // Electric energy consumed by the evap condenser water circulation pump [J]
        Real64 EvapWaterConsump = 0.0;            // Evap condenser water consumption [m3]
        Real64 BasinHeaterConsumption = 0.0;      // Basin heater energy consumption (J)
        Real64 BasinHeaterPowerFTempDiff = 0.0;   // Basin heater capacity per degree C below setpoint (W/C)
        Real64 BasinHeaterSetPointTemp = 0.0;     // setpoint temperature for basin heater operation (C)
        Real64 BasinHeaterPower = 0.0;            // Basin heater power (W)
        int BasinHeaterSchedulePtr = 0;         // Pointer to basin heater schedule
        Array1D<Real64> EvapCondAirFlow;    // Air flow rate through the evap condenser at high speed, volumetric flow rate
        // for water use calcs [m3/s]
        Array1D<Real64> EvapCondEffect; // effectiveness of the evaporatively cooled condenser
        // [high speed for multi-speed unit] (-)
        Array1D<Real64> MSRatedEvapCondVolFlowPerRatedTotCap; // evap condenser air flow ratio to capacity
        // begin variables for Water System interactions
        int EvapWaterSupplyMode = 101;         // where does water come from
        std::string EvapWaterSupplyName; // name of water source e.g. water storage tank
        int EvapWaterSupTankID = 0;
        int EvapWaterTankDemandARRID = 0;
        int CondensateCollectMode = 1001;         // where does water come from
        std::string CondensateCollectName; // name of water source e.g. water storage tank
        int CondensateTankID = 0;
        int CondensateTankSupplyARRID = 0;
        Real64 CondensateVdot = 0.0;         // rate of water condensation from air stream [m3/s]
        Real64 CondensateVol = 0.0;          // amount of water condensed from air stream [m3]
        Real64 CondInletTemp = 0.0;          // Evap condenser inlet temperature [C], report variable
        int SupplyFanIndex = 0;            // index of this fan in fan array or vector
        int SupplyFan_TypeNum = 0;         // type of fan, in DataHVACGlobals
        std::string SupplyFanName;     // name of fan associated with this dx coil
        Real64 SourceAirMassFlowRate = 0.0;  // source air mass flow rate [kg/s]
        Real64 InletSourceAirTemp = 0.0;     // source air temperature entering the outdoor coil [C]
        Real64 InletSourceAirEnthalpy = 0.0; // source air enthalpy entering the outdoor coil [J/kg]
        // end variables for water system interactions

        // begin varibles for HPWH
        Real64 RatedCapWH = 0.0;                  // Rated water heating Capacity [W]
        int InletAirTemperatureType = 0;        // Specifies to use either air wet-bulb or dry-bulb temp for curve objects
        Real64 WHRatedInletDBTemp = 0.0;          // Rated inlet air dry-bulb temperature [C]
        Real64 WHRatedInletWBTemp = 0.0;          // Rated inlet air wet-bulb temperature [C]
        Real64 WHRatedInletWaterTemp = 0.0;       // Rated condenser water inlet temperature [C]
        Real64 HPWHCondPumpElecNomPower = 0.0;    // Nominal power input to the condenser water circulation pump [W]
        Real64 HPWHCondPumpFracToWater = 1.0;     // Nominal power fraction to water for the condenser water circulation pump
        Real64 RatedHPWHCondWaterFlow = 0.0;      // Rated water flow rate through the condenser of the HPWH DX coil [m3/s]
        Real64 ElecWaterHeatingPower = 0.0;       // Total electric power consumed by compressor and condenser pump [W]
        Real64 ElecWaterHeatingConsumption = 0.0; // Total electric consumption by compressor and condenser pump [J]
        bool FanPowerIncludedInCOP = false;         // Indicates that fan heat is included in heating capacity and COP
        bool CondPumpHeatInCapacity = false;        // Indicates that condenser pump heat is included in heating capacity
        bool CondPumpPowerInCOP = false;            // Indicates that condenser pump power is included in heating COP
        bool AirVolFlowAutoSized = false;           // Used to report autosizing info for the HPWH DX coil
        bool WaterVolFlowAutoSized = false;         // Used to report autosizing info for the HPWH DX coil
        Real64 TotalHeatingEnergy = 0.0;          // total water heating energy
        Real64 TotalHeatingEnergyRate = 0.0;      // total WH energy rate
        bool bIsDesuperheater = false;              // whether the coil is used for a desuperheater, i.e. zero all the cooling capacity and power
        // end variables for HPWH
        bool reportCoilFinalSizes = true; // one time report of sizes to coil selection report
        Real64 capModFacTotal = 0.0;     // coil  TotCapTempModFac * TotCapAirFFModFac * TotCapWaterFFModFac, for result for simulation peak reporting

        // Default Constructor
        VariableSpeedCoilData() 
          : MSErrIndex(MaxSpedLevels, 0),  MSRatedPercentTotCap(MaxSpedLevels, 0.0), MSRatedTotCap(MaxSpedLevels, 0.0),
            MSRatedSHR(MaxSpedLevels, 0.0), MSRatedCOP(MaxSpedLevels, 0.0), MSRatedAirVolFlowPerRatedTotCap(MaxSpedLevels, 0.0),
            MSRatedAirVolFlowRate(MaxSpedLevels, 0.0), MSRatedAirMassFlowRate(MaxSpedLevels, 0.0),
            MSRatedWaterVolFlowPerRatedTotCap(MaxSpedLevels, 0.0), MSRatedWaterVolFlowRate(MaxSpedLevels, 0.0),
            MSRatedWaterMassFlowRate(MaxSpedLevels, 0.0), MSRatedCBF(MaxSpedLevels, 0.0), MSEffectiveAo(MaxSpedLevels, 0.0),
            MSCCapFTemp(MaxSpedLevels, 0), MSCCapAirFFlow(MaxSpedLevels, 0), MSCCapWaterFFlow(MaxSpedLevels, 0), MSEIRFTemp(MaxSpedLevels, 0),
            MSEIRAirFFlow(MaxSpedLevels, 0), MSEIRWaterFFlow(MaxSpedLevels, 0), MSWasteHeat(MaxSpedLevels, 0), MSWasteHeatFrac(MaxSpedLevels, 0.0),
            MSWHPumpPower(MaxSpedLevels, 0.0), MSWHPumpPowerPerRatedTotCap(MaxSpedLevels, 0.0), EvapCondAirFlow(MaxSpedLevels, 0.0), EvapCondEffect(MaxSpedLevels, 0.0),
            MSRatedEvapCondVolFlowPerRatedTotCap(MaxSpedLevels, 0.0)
        {
        }
    };

    void SimVariableSpeedCoils(EnergyPlusData &state, std::string const &CompName,   // Coil Name
                               int &CompIndex,                // Index for Component name
                               int const CyclingScheme,       // Continuous fan OR cycling compressor
                               Real64 &MaxONOFFCyclesperHour, // Maximum cycling rate of heat pump [cycles/hr]
                               Real64 &HPTimeConstant,        // Heat pump time constant [s]
                               Real64 &FanDelayTime,          // Fan delay time, time delay for the HP's fan to
                               int const CompOp,              // compressor on/off. 0 = off; 1= on
                               Real64 const PartLoadFrac,
                               int const SpeedNum,                        // compressor speed number
                               Real64 const SpeedRatio,                   // compressor speed ratio
                               Real64 const SensLoad,                     // Sensible demand load [W]
                               Real64 const LatentLoad,                   // Latent demand load [W]
                               Optional<Real64 const> OnOffAirFlowRat = _ // ratio of comp on to comp off air flow rate
    );

    void GetVarSpeedCoilInput(EnergyPlusData &state);

    // Beginning Initialization Section of the Module
    //******************************************************************************

    void InitVarSpeedCoil(EnergyPlusData &state, int const DXCoilNum,                // Current DXCoilNum under simulation
                          Real64 const MaxONOFFCyclesperHour, // Maximum cycling rate of heat pump [cycles/hr]
                          Real64 const HPTimeConstant,        // Heat pump time constant [s]
                          Real64 const FanDelayTime,          // Fan delay time, time delay for the HP's fan to
                          Real64 const SensLoad,              // Control zone sensible load[W]
                          Real64 const LatentLoad,            // Control zone latent load[W]
                          int const CyclingScheme,            // fan operating mode
                          Real64 const OnOffAirFlowRatio,     // ratio of compressor on flow to average flow over time step
                          Real64 const SpeedRatio,            // compressor speed ratio
                          int const SpeedNum                  // compressor speed number
    );

    void SizeVarSpeedCoil(EnergyPlusData &state, int const DXCoilNum);

    void CalcVarSpeedCoilCooling(EnergyPlusData &state,
                                 int const DXCoilNum,            // Heat Pump Number
                                 int const CyclingScheme,        // Fan/Compressor cycling scheme indicator
                                 Real64 &RuntimeFrac,            // Runtime Fraction of compressor or percent on time (on-time/cycle time)
                                 Real64 const SensDemand,        // Cooling Sensible Demand [W] !unused1208
                                 Real64 const LatentDemand,      // Cooling Latent Demand [W]
                                 int const CompOp,               // compressor operation flag
                                 Real64 const PartLoadRatio,     // compressor part load ratio
                                 Real64 const OnOffAirFlowRatio, // ratio of compressor on flow to average flow over time step
                                 Real64 const SpeedRatio,        // SpeedRatio varies between 1.0 (higher speed) and 0.0 (lower speed)
                                 int const SpeedNum              // Speed number, high bound
    );

    void CalcVarSpeedCoilHeating(EnergyPlusData &state,
                                 int const DXCoilNum,            // Heat Pump Number
                                 int const CyclingScheme,        // Fan/Compressor cycling scheme indicator
                                 Real64 &RuntimeFrac,            // Runtime Fraction of compressor or percent on time (on-time/cycle time)
                                 Real64 const SensDemand,        // Cooling Sensible Demand [W] !unused1208
                                 int const CompOp,               // compressor operation flag
                                 Real64 const PartLoadRatio,     // compressor part load ratio
                                 Real64 const OnOffAirFlowRatio, // ratio of compressor on flow to average flow over time step
                                 Real64 const SpeedRatio,        // SpeedRatio varies between 1.0 (higher speed) and 0.0 (lower speed)
                                 int const SpeedNum              // Speed number, high bound, i.e. SpeedNum - 1 is the other side
    );

    Real64 GetCoilCapacityVariableSpeed(EnergyPlusData &state,
                                        std::string const &CoilType, // must match coil types in this module
                                        std::string const &CoilName, // must match coil names for the coil type
                                        bool &ErrorsFound            // set to true if problem
    );

    int GetCoilIndexVariableSpeed(EnergyPlusData &state,
                                  std::string const &CoilType, // must match coil types in this module
                                  std::string const &CoilName, // must match coil names for the coil type
                                  bool &ErrorsFound            // set to true if problem
    );

    Real64 GetCoilAirFlowRateVariableSpeed(EnergyPlusData &state,
                                           std::string const &CoilType, // must match coil types in this module
                                           std::string const &CoilName, // must match coil names for the coil type
                                           bool &ErrorsFound            // set to true if problem
    );

    int GetCoilInletNodeVariableSpeed(EnergyPlusData &state,
                                      std::string const &CoilType, // must match coil types in this module
                                      std::string const &CoilName, // must match coil names for the coil type
                                      bool &ErrorsFound            // set to true if problem
    );

    int GetCoilOutletNodeVariableSpeed(EnergyPlusData &state,
                                       std::string const &CoilType, // must match coil types in this module
                                       std::string const &CoilName, // must match coil names for the coil type
                                       bool &ErrorsFound            // set to true if problem
    );

    int GetVSCoilCondenserInletNode(EnergyPlusData &state,
                                    std::string const &CoilName, // must match coil names for the coil type
                                    bool &ErrorsFound            // set to true if problem
    );

    int GetVSCoilPLFFPLR(EnergyPlusData &state,
                         std::string const &CoilType, // must match coil types in this module
                         std::string const &CoilName, // must match coil names for the coil type
                         bool &ErrorsFound            // set to true if problem
    );

    int GetVSCoilCapFTCurveIndex(EnergyPlusData &state,
                                 int const &CoilIndex, // must match coil names for the coil type
                                 bool &ErrorsFound     // set to true if problem
    );

    Real64 GetVSCoilMinOATCompressor(EnergyPlusData &state,
                                     std::string const &CoilName, // must match coil names for the coil type
                                     bool &ErrorsFound            // set to true if problem
    );

    Real64 GetVSCoilMinOATCompressorUsingIndex(EnergyPlusData &state,
                                               int const CoilIndex, // index to cooling coil
                                               bool &ErrorsFound    // set to true if problem
    );

    int GetVSCoilNumOfSpeeds(EnergyPlusData &state,
                             std::string const &CoilName, // must match coil names for the coil type
                             bool &ErrorsFound            // set to true if problem
    );

    void SetVarSpeedCoilData(EnergyPlusData &state,
                             int const WSHPNum,                        // Number of OA Controller
                             bool &ErrorsFound,                        // Set to true if certain errors found
                             Optional_int CompanionCoolingCoilNum = _, // Index to cooling coil for heating coil = SimpleWSHPNum
                             Optional_int CompanionHeatingCoilNum = _, // Index to heating coil for cooling coil = SimpleWSHPNum
                             Optional_int MSHPDesignSpecIndex = _      // index to UnitarySystemPerformance:Multispeed object
    );

    void UpdateVarSpeedCoil(EnergyPlusData &state, int const DXCoilNum);

    Real64 CalcEffectiveSHR(EnergyPlusData &state, int const DXCoilNum,     // Index number for cooling coil
                            Real64 const SHRss,      // Steady-state sensible heat ratio
                            int const CyclingScheme, // Fan/compressor cycling scheme indicator
                            Real64 const RTF,        // Compressor run-time fraction
                            Real64 const QLatRated,  // Rated latent capacity
                            Real64 const QLatActual, // Actual latent capacity
                            Real64 const EnteringDB, // Entering air dry-bulb temperature
                            Real64 const EnteringWB  // Entering air wet-bulb temperature
    );

    void CalcTotCapSHR_VSWSHP(EnergyPlusData &state,
                              Real64 const InletDryBulb,       // inlet air dry bulb temperature [C]
                              Real64 const InletHumRat,        // inlet air humidity ratio [kg water / kg dry air]
                              Real64 const InletEnthalpy,      // inlet air specific enthalpy [J/kg]
                              Real64 &InletWetBulb,            // inlet air wet bulb temperature [C]
                              Real64 const AirMassFlowRatio,   // Ratio of actual air mass flow to nominal air mass flow
                              Real64 const WaterMassFlowRatio, // Ratio of actual water mass flow to nominal water mass flow
                              Real64 const AirMassFlow,        // actual mass flow for capacity and SHR calculation
                              Real64 const CBF,                // coil bypass factor
                              Real64 const TotCapNom1,         // nominal total capacity at low speed [W]
                              int const CCapFTemp1,            // capacity modifier curve index, function of entering wetbulb at low speed
                              int const CCapAirFFlow1,         // capacity modifier curve, function of actual air flow vs rated flow at low speed
                              int const CCapWaterFFlow1,       // capacity modifier curve, function of actual water flow vs rated flow at low speed
                              Real64 const TotCapNom2,         // nominal total capacity at high speed [W]
                              int const CCapFTemp2,            // capacity modifier curve index, function of entering wetbulb at high speed
                              int const CCapAirFFlow2,         // capacity modifier curve, function of actual air flow vs rated flow at high speed
                              int const CCapWaterFFlow2,       // capacity modifier curve, function of actual water flow vs rated flow at high speed
                              Real64 &TotCap1,                 // total capacity at the given conditions [W] at low speed
                              Real64 &TotCap2,                 // total capacity at the given conditions [W] at high speed
                              Real64 &TotCapSpeed,             // integrated total capacity corresponding to the speed ratio
                              Real64 &SHR,                     // sensible heat ratio at the given conditions
                              Real64 const CondInletTemp,      // Condenser inlet temperature [C]
                              Real64 const Pressure,           // air pressure [Pa]
                              Real64 const SpeedRatio,         // from 0.0 to 1.0
                              int const NumSpeeds,             // number of speeds for input
                              Real64 &TotCapModFac             // capacity modification factor, func of temp and func of flow
    );

    void CalcVarSpeedHPWH(EnergyPlusData &state,
                          int const DXCoilNum,        // the number of the DX coil to be simulated
                          Real64 &RuntimeFrac,        // Runtime Fraction of compressor or percent on time (on-time/cycle time)
                          Real64 const PartLoadRatio, // sensible water heating load / full load sensible water heating capacity
                          Real64 const SpeedRatio,    // SpeedRatio varies between 1.0 (higher speed) and 0.0 (lower speed)
                          int const SpeedNum,         // Speed number, high bound capacity
                          int const CyclingScheme     // Continuous fan OR cycling compressor
    );

    Real64 getVarSpeedPartLoadRatio(EnergyPlusData &state, int const DXCoilNum); // the number of the DX coil to mined for current PLR

    void setVarSpeedHPWHFanTypeNum(EnergyPlusData &state, int const dXCoilNum, int const fanTypeNum);

    void setVarSpeedHPWHFanIndex(EnergyPlusData &state, int const dXCoilNum, int const fanIndex);

    void setVarSpeedFanInfo(EnergyPlusData &state, int const dXCoilNum, std::string const fanName, int const fanIndex, int const fanTypeNum);

} // namespace VariableSpeedCoils

struct VariableSpeedCoilsData : BaseGlobalStruct {

    Real64 const RatedInletAirTemp = 26.6667;       // 26.6667C or 80F
    Real64 const RatedInletWetBulbTemp = 19.4444;   // 19.44 or 67F, cooling mode
    Real64 const RatedInletAirHumRat = 0.0111847;   // Humidity ratio corresponding to 80F dry bulb/67F wet bulb
    Real64 const RatedInletWaterTemp = 29.4444;     // 85 F cooling mode
    Real64 const RatedAmbAirTemp = 35.0;            // 95 F cooling mode
    Real64 const RatedInletAirTempHeat = 21.1111;   // 21.11C or 70F, heating mode
    Real64 const RatedInletWaterTempHeat = 21.1111; // 21.11C or 70F, heating mode
    Real64 const RatedAmbAirTempHeat = 8.3333;      // 8.33 or 47F, heating mode
    Real64 const RatedAmbAirWBHeat = 6.1111;        // 8.33 or 43F, heating mode, rated wet bulb temperature
                                                       // Water Systems
    int const CondensateDiscarded = 1001; // default mode where water is "lost"
    int const CondensateToTank = 1002;    // collect coil condensate from air and store in water storage tank

    int const WaterSupplyFromMains = 101;
    int const WaterSupplyFromTank = 102;

    // Curve Types
    int const Linear = 1;
    int const BiLinear = 2;
    int const Quadratic = 3;
    int const BiQuadratic = 4;
    int const Cubic = 5;

    // Defrost strategy (heat pump only)
    int const ReverseCycle = 1; // uses reverse cycle defrost strategy
    int const Resistive = 2;    // uses electric resistance heater for defrost
                               // Defrost control  (heat pump only)
    int const Timed = 1;    // defrost cycle is timed
    int const OnDemand = 2; // defrost cycle occurs only when required

    int const MaxSpedLevels = 10; // Maximum number of speed that supports  

    int NumVarSpeedCoils = 0; // The Number of Water to Air Heat Pumps found in the Input

    bool MyOneTimeFlag = true;     // one time allocation flag
    bool GetCoilsInputFlag = true; // Flag set to make sure you get input once
                                  // LOGICAL, ALLOCATABLE, DIMENSION(:) :: MySizeFlag

    Real64 SourceSideMassFlowRate = 0.0; // Source Side Mass flow rate [Kg/s]
    Real64 SourceSideInletTemp = 0.0;    // Source Side Inlet Temperature [C]
    Real64 SourceSideInletEnth = 0.0;    // Source Side Inlet Enthalpy [J/kg]
    Real64 LoadSideMassFlowRate = 0.0;   // Load Side Mass flow rate [Kg/s]
    Real64 LoadSideInletDBTemp = 0.0;    // Load Side Inlet Dry Bulb Temp [C]
    Real64 LoadSideInletWBTemp = 0.0;    // Load Side Inlet Wet Bulb Temp [C]
    Real64 LoadSideInletHumRat = 0.0;    // Load Side Outlet Humidity ratio
    Real64 LoadSideInletEnth = 0.0;      // Load Side Inlet Enthalpy [J/kg]
    Real64 LoadSideOutletDBTemp = 0.0;   // Load Side Outlet Dry Bulb Temp [C]
    Real64 LoadSideOutletHumRat = 0.0;   // Load Side Outlet Humidity ratio
    Real64 LoadSideOutletEnth = 0.0;     // Load Side Outlet Enthalpy [J/kg]
    Real64 QSensible = 0.0;              // Load side sensible heat transfer rate [W]
    Real64 QLoadTotal = 0.0;             // Load side total heat transfer rate [W]
    Real64 QLatRated = 0.0;              // Latent Capacity [W] rated at entering air conditions [Tdb=26.7C Twb=19.4C]
    Real64 QLatActual = 0.0;             // Actual Latent Capacity [W]
    Real64 QSource = 0.0;                // Source side heat transfer rate [W]
    Real64 Winput = 0.0;                 // Power Consumption [W]
    Real64 PLRCorrLoadSideMdot = 0.0;    // Load Side Mdot corrected for Part Load Ratio of the unit

    Real64 VSHPWHHeatingCapacity = 0.0; // Used by Heat Pump:Water Heater object as total water heating capacity [W]
    Real64 VSHPWHHeatingCOP = 0.0;      // Used by Heat Pump:Water Heater object as water heating COP [W/W]
    
    Array1D<VariableSpeedCoils::VariableSpeedCoilData> VarSpeedCoil;
    bool firstTime = true;

    void clear_state() override
    {
        NumVarSpeedCoils = 0;
        MyOneTimeFlag = true;
        GetCoilsInputFlag = true;
        SourceSideMassFlowRate = 0.0;
        SourceSideInletTemp = 0.0;
        SourceSideInletEnth = 0.0;
        LoadSideMassFlowRate = 0.0;
        LoadSideInletDBTemp = 0.0;
        LoadSideInletWBTemp = 0.0;
        LoadSideInletHumRat = 0.0;
        LoadSideInletEnth = 0.0;
        LoadSideOutletDBTemp = 0.0;
        LoadSideOutletHumRat = 0.0;
        LoadSideOutletEnth = 0.0;
        QSensible = 0.0;
        QLoadTotal = 0.0;
        QLatRated = 0.0;
        QLatActual = 0.0;
        QSource = 0.0;
        Winput = 0.0;
        PLRCorrLoadSideMdot = 0.0;
        VSHPWHHeatingCapacity = 0.0;
        VSHPWHHeatingCOP = 0.0;
        VarSpeedCoil.deallocate();
        firstTime = true;
    }

    // Default Constructor
    VariableSpeedCoilsData()
    {
    }
};
} // namespace EnergyPlus

#endif
