// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Contains the classes required to extend the black-oil model by polymer.
 */
#ifndef EWOMS_BLACK_OIL_POLYMER_MODULE_HH
#define EWOMS_BLACK_OIL_POLYMER_MODULE_HH

#include "blackoilproperties.hh"
#include <opm/models/io/vtkblackoilpolymermodule.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/IntervalTabulated2DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyadsTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlymaxTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyrockTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyshlogTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyviscTable.hpp>
#endif

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by polymer.
 */
template <class TypeTag, bool enablePolymerV = getPropValue<TypeTag, Properties::EnablePolymer>()>
class BlackOilPolymerModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Toolbox = Opm::MathToolbox<Evaluation>;

    using TabulatedFunction = typename Opm::Tabulated1DFunction<Scalar>;
    using TabulatedTwoDFunction = typename Opm::IntervalTabulated2DFunction<Scalar>;

    static constexpr unsigned polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr unsigned polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
    static constexpr unsigned contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr unsigned contiPolymerMolarWeightEqIdx = Indices::contiPolymerMWEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;


    static constexpr unsigned enablePolymer = enablePolymerV;
    static constexpr bool enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr bool enablePolymerMechanicalDegradation = getPropValue<TypeTag, Properties::EnablePolymerMechanicalDegradation>();

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;

    struct SkprpolyTable {
        double refConcentration;
        TabulatedTwoDFunction table_func;
    };

public:
    enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

    // a struct containing the constants to calculate polymer viscosity
    // based on Mark-Houwink equation and Huggins equation, the constants are provided
    // by the keyword PLYVMH
    struct PlyvmhCoefficients {
        Scalar k_mh;
        Scalar a_mh;
        Scalar gamma;
        Scalar kappa;
    };
    
    struct StashedValues {
        Dune::FieldVector<Scalar, 6> velocity; // at faces
        //Dune::FieldVector<Scalar, 6> distance; // cell centre --> face
        Scalar average_cell_velocity; // m/sec (seepage velocity)
        Scalar cell_shear_rate;
        
        Scalar debug_value;
    };
    
    static std::map<int, StashedValues> stashedValues_;

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the polymer module
     */
    static void initFromState(const Opm::EclipseState& eclState)
    {
        // some sanity checks: if polymers are enabled, the POLYMER keyword must be
        // present, if polymers are disabled the keyword must not be present.
        if (enablePolymer && !eclState.runspec().phases().active(Phase::POLYMER)) {
            throw std::runtime_error("Non-trivial polymer treatment requested at compile time, but "
                                     "the deck does not contain the POLYMER keyword");
        }
        else if (!enablePolymer && eclState.runspec().phases().active(Phase::POLYMER)) {
            throw std::runtime_error("Polymer treatment disabled at compile time, but the deck "
                                     "contains the POLYMER keyword");
        }

        if (enablePolymerMolarWeight && !eclState.runspec().phases().active(Phase::POLYMW)) {
            throw std::runtime_error("Polymer molecular weight tracking is enabled at compile time, but "
                                     "the deck does not contain the POLYMW keyword");
        }
        else if (!enablePolymerMolarWeight && eclState.runspec().phases().active(Phase::POLYMW)) {
            throw std::runtime_error("Polymer molecular weight tracking is disabled at compile time, but the deck "
                                     "contains the POLYMW keyword");
        }
        
        if (enablePolymerMechanicalDegradation && !eclState.runspec().phases().active(Phase::POLYMW)) {
            throw std::runtime_error("Polymer mechanical degradation model is enabled at compile time, but "
                                     "the deck does not contain the POLYMW keyword");
        }

        if (enablePolymerMolarWeight && !enablePolymer) {
            throw std::runtime_error("Polymer molecular weight tracking is enabled while polymer treatment "
                                     "is disabled at compile time");
        }
        
        if (enablePolymerMechanicalDegradation && !enablePolymer) {
            throw std::runtime_error("Polymer mechanical degradation model is enabled while polymer treatment "
                                     "is disabled at compile time");
        }

        if (!eclState.runspec().phases().active(Phase::POLYMER))
            return; // polymer treatment is supposed to be disabled

        const auto& tableManager = eclState.getTableManager();

        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        setNumSatRegions(numSatRegions);

        // initialize the objects which deal with the PLYROCK keyword
        const auto& plyrockTables = tableManager.getPlyrockTables();
        if (!plyrockTables.empty()) {
            assert(numSatRegions == plyrockTables.size());
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                const auto& plyrockTable = plyrockTables.template getTable<Opm::PlyrockTable>(satRegionIdx);
                setPlyrock(satRegionIdx,
                           plyrockTable.getDeadPoreVolumeColumn()[0],
                           plyrockTable.getResidualResistanceFactorColumn()[0],
                           plyrockTable.getRockDensityFactorColumn()[0],
                           static_cast<AdsorptionBehaviour>(plyrockTable.getAdsorbtionIndexColumn()[0]),
                           plyrockTable.getMaxAdsorbtionColumn()[0]);
            }
        }
        else {
            throw std::runtime_error("PLYROCK must be specified in POLYMER runs\n");
        }

        // initialize the objects which deal with the PLYADS keyword
        const auto& plyadsTables = tableManager.getPlyadsTables();
        if (!plyadsTables.empty()) {
            assert(numSatRegions == plyadsTables.size());
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                const auto& plyadsTable = plyadsTables.template getTable<Opm::PlyadsTable>(satRegionIdx);
                // Copy data
                const auto& c = plyadsTable.getPolymerConcentrationColumn();
                const auto& ads = plyadsTable.getAdsorbedPolymerColumn();
                plyadsAdsorbedPolymer_[satRegionIdx].setXYContainers(c, ads);
            }
        }
        else {
            throw std::runtime_error("PLYADS must be specified in POLYMER runs\n");
        }


        unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
        setNumPvtRegions(numPvtRegions);

        // initialize the objects which deal with the PLYVISC keyword
        const auto& plyviscTables = tableManager.getPlyviscTables();
        if (!plyviscTables.empty()) {
            // different viscosity model is used for POLYMW
            if (enablePolymerMolarWeight && !enablePolymerMechanicalDegradation) {
                Opm::OpmLog::warning("PLYVISC should not be used in POLYMW runs, "
                                     "it will have no effect. A viscosity model based on PLYVMH is used instead.\n");
            }
            else {
                assert(numPvtRegions == plyviscTables.size());
                for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                    const auto& plyadsTable = plyviscTables.template getTable<Opm::PlyviscTable>(pvtRegionIdx);
                    // Copy data
                    const auto& c = plyadsTable.getPolymerConcentrationColumn();
                    const auto& visc = plyadsTable.getViscosityMultiplierColumn();
                    plyviscViscosityMultiplierTable_[pvtRegionIdx].setXYContainers(c, visc);
                }
            }
            
            if(enablePolymerMolarWeight && enablePolymerMechanicalDegradation)
            {
                Opm::OpmLog::warning("Using experimental mechanical degradation model...\n");
            }
            
        }
        else if (!enablePolymerMolarWeight) {
            throw std::runtime_error("PLYVISC must be specified in POLYMER runs\n");
        }

        // initialize the objects which deal with the PLYMAX keyword
        const auto& plymaxTables = tableManager.getPlymaxTables();
        const unsigned numMixRegions = plymaxTables.size();
        setNumMixRegions(numMixRegions);
        if (!plymaxTables.empty()) {
            for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                const auto& plymaxTable = plymaxTables.template getTable<Opm::PlymaxTable>(mixRegionIdx);
                setPlymax(mixRegionIdx, plymaxTable.getPolymerConcentrationColumn()[0]);
            }
        }
        else {
            throw std::runtime_error("PLYMAX must be specified in POLYMER runs\n");
        }

        if (!eclState.getTableManager().getPlmixparTable().empty()) {
            if (enablePolymerMolarWeight) {
                Opm::OpmLog::warning("PLMIXPAR should not be used in POLYMW runs, it will have no effect.\n");
            }
            else {
                const auto& plmixparTable = eclState.getTableManager().getPlmixparTable();
                // initialize the objects which deal with the PLMIXPAR keyword
                for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                    setPlmixpar(mixRegionIdx, plmixparTable[mixRegionIdx].todd_langstaff);
                }
            }
        }
        else if (!enablePolymerMolarWeight) {
            throw std::runtime_error("PLMIXPAR must be specified in POLYMER runs\n");
        }

        hasPlyshlog_ = eclState.getTableManager().hasTables("PLYSHLOG");
        hasShrate_ = eclState.getTableManager().useShrate();

        if ((hasPlyshlog_ || hasShrate_) && enablePolymerMolarWeight) {
            Opm::OpmLog::warning("PLYSHLOG and SHRATE should not be used in POLYMW runs, they will have no effect.\n");
        }

        if (hasPlyshlog_ && !enablePolymerMolarWeight) {
            const auto& plyshlogTables = tableManager.getPlyshlogTables();
            assert(numPvtRegions == plyshlogTables.size());
            plyshlogShearEffectRefMultiplier_.resize(numPvtRegions);
            plyshlogShearEffectRefLogVelocity_.resize(numPvtRegions);
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                const auto& plyshlogTable = plyshlogTables.template getTable<Opm::PlyshlogTable>(pvtRegionIdx);

                Scalar plyshlogRefPolymerConcentration = plyshlogTable.getRefPolymerConcentration();
                auto waterVelocity = plyshlogTable.getWaterVelocityColumn().vectorCopy();
                auto shearMultiplier = plyshlogTable.getShearMultiplierColumn().vectorCopy();

                // do the unit version here for the waterVelocity
                Opm::UnitSystem unitSystem = eclState.getDeckUnitSystem();
                double siFactor = hasShrate_? unitSystem.parse("1/Time").getSIScaling() : unitSystem.parse("Length/Time").getSIScaling();
                for (size_t i = 0; i < waterVelocity.size(); ++i) {
                    waterVelocity[i] *= siFactor;
                    // for plyshlog the input must be stored as logarithms
                    // the interpolation is then done the log-space.
                    waterVelocity[i] = std::log(waterVelocity[i]);
                }

                Scalar refViscMult = plyviscViscosityMultiplierTable_[pvtRegionIdx].eval(plyshlogRefPolymerConcentration, /*extrapolate=*/true);
                // convert the table using referece conditions
                for (size_t i = 0; i < waterVelocity.size(); ++i) {
                    shearMultiplier[i] *= refViscMult;
                    shearMultiplier[i] -= 1;
                    shearMultiplier[i] /= (refViscMult - 1);
                    shearMultiplier[i] = shearMultiplier[i];
                }
                plyshlogShearEffectRefMultiplier_[pvtRegionIdx].resize(waterVelocity.size());
                plyshlogShearEffectRefLogVelocity_[pvtRegionIdx].resize(waterVelocity.size());

                for (size_t i = 0; i < waterVelocity.size(); ++i) {
                    plyshlogShearEffectRefMultiplier_[pvtRegionIdx][i] = shearMultiplier[i];
                    plyshlogShearEffectRefLogVelocity_[pvtRegionIdx][i] = waterVelocity[i];
                }
            }
        }

        if (hasShrate_ && !enablePolymerMolarWeight) {
            if(!hasPlyshlog_) {
                throw std::runtime_error("PLYSHLOG must be specified if SHRATE is used in POLYMER runs\n");
            }
            const auto& shrateTable = eclState.getTableManager().getShrateTable();
            shrate_.resize(numPvtRegions);
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                if (shrateTable.empty()) {
                    shrate_[pvtRegionIdx] = 4.8; //default;
                }
                else if (shrateTable.size() == numPvtRegions) {
                    shrate_[pvtRegionIdx] = shrateTable[pvtRegionIdx].rate;
                }
                else {
                    throw std::runtime_error("SHRATE must either have 0 or number of NUMPVT entries\n");
                }
            }
        }

        if (enablePolymerMolarWeight) { // include for both kinds of POLYMW models
            const auto& plyvmhTable = eclState.getTableManager().getPlyvmhTable();
            if (!plyvmhTable.empty()) {
                assert(plyvmhTable.size() == numMixRegions);
                for (size_t regionIdx = 0; regionIdx < numMixRegions; ++regionIdx) {
                    plyvmhCoefficients_[regionIdx].k_mh = plyvmhTable[regionIdx].k_mh;
                    plyvmhCoefficients_[regionIdx].a_mh = plyvmhTable[regionIdx].a_mh;
                    plyvmhCoefficients_[regionIdx].gamma = plyvmhTable[regionIdx].gamma;
                    plyvmhCoefficients_[regionIdx].kappa = plyvmhTable[regionIdx].kappa;
                }
            }
            else {
                throw std::runtime_error("PLYVMH keyword must be specified in POLYMW runs \n");
            }
        }
        
        if (enablePolymerMolarWeight){
            
            // handling PLYMWINJ keyword
            const auto& plymwinjTables = tableManager.getPlymwinjTables();
            for (const auto& table : plymwinjTables) {
                const int tableNumber = table.first;
                const auto& plymwinjtable = table.second;
                const std::vector<double>& throughput = plymwinjtable.getThroughputs();
                const std::vector<double>& watervelocity = plymwinjtable.getVelocities();
                const std::vector<std::vector<double>>& molecularweight = plymwinjtable.getMoleWeights();
                TabulatedTwoDFunction tablefunc(throughput, watervelocity, molecularweight, true, false);
                plymwinjTables_[tableNumber] = std::move(tablefunc);
            }

            // handling SKPRWAT keyword
            const auto& skprwatTables = tableManager.getSkprwatTables();
            for (const auto& table : skprwatTables) {
                const int tableNumber = table.first;
                const auto& skprwattable = table.second;
                const std::vector<double>& throughput = skprwattable.getThroughputs();
                const std::vector<double>& watervelocity = skprwattable.getVelocities();
                const std::vector<std::vector<double>>& skinpressure = skprwattable.getSkinPressures();
                TabulatedTwoDFunction tablefunc(throughput, watervelocity, skinpressure, true, false);
                skprwatTables_[tableNumber] = std::move(tablefunc);
            }

            // handling SKPRPOLY keyword
            const auto& skprpolyTables = tableManager.getSkprpolyTables();
            for (const auto& table : skprpolyTables) {
                const int tableNumber = table.first;
                const auto& skprpolytable = table.second;
                const std::vector<double>& throughput = skprpolytable.getThroughputs();
                const std::vector<double>& watervelocity = skprpolytable.getVelocities();
                const std::vector<std::vector<double>>& skinpressure = skprpolytable.getSkinPressures();
                const double refPolymerConcentration = skprpolytable.referenceConcentration();
                SkprpolyTable tablefunc = {refPolymerConcentration,
                                           TabulatedTwoDFunction(throughput, watervelocity, skinpressure, true, false)};
                skprpolyTables_[tableNumber] = std::move(tablefunc);
            }
        }
        /*
        else if(enablePolymerMolarWeight && enablePolymerMechanicalDegradation)
        {
            // TO DO: Read needed data from input file; need to implement new keyword(s)!
        }
        */
    }
#endif

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the PLYROCK and PLYADS of any region.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        plyrockDeadPoreVolume_.resize(numRegions);
        plyrockResidualResistanceFactor_.resize(numRegions);
        plyrockRockDensityFactor_.resize(numRegions);
        plyrockAdsorbtionIndex_.resize(numRegions);
        plyrockMaxAdsorbtion_.resize(numRegions);
        plyadsAdsorbedPolymer_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer rock properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setPlyrock(unsigned satRegionIdx,
                           const Scalar& plyrockDeadPoreVolume,
                           const Scalar& plyrockResidualResistanceFactor,
                           const Scalar& plyrockRockDensityFactor,
                           const Scalar& plyrockAdsorbtionIndex,
                           const Scalar& plyrockMaxAdsorbtion)
    {
        plyrockDeadPoreVolume_[satRegionIdx] = plyrockDeadPoreVolume;
        plyrockResidualResistanceFactor_[satRegionIdx] = plyrockResidualResistanceFactor;
        plyrockRockDensityFactor_[satRegionIdx] = plyrockRockDensityFactor;
        plyrockAdsorbtionIndex_[satRegionIdx] = plyrockAdsorbtionIndex;
        plyrockMaxAdsorbtion_[satRegionIdx] = plyrockMaxAdsorbtion;
    }

    /*!
     * \brief Specify the number of pvt regions.
     *
     * This must be called before setting the PLYVISC of any region.
     */
    static void setNumPvtRegions(unsigned numRegions)
    {
        plyviscViscosityMultiplierTable_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer viscosity a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setPlyvisc(unsigned satRegionIdx,
                           const TabulatedFunction& plyviscViscosityMultiplierTable)
    {
        plyviscViscosityMultiplierTable_[satRegionIdx] = plyviscViscosityMultiplierTable;
    }

    /*!
     * \brief Specify the number of mix regions.
     *
     * This must be called before setting the PLYMAC and PLMIXPAR of any region.
     */
    static void setNumMixRegions(unsigned numRegions)
    {
        plymaxMaxConcentration_.resize(numRegions);
        plymixparToddLongstaff_.resize(numRegions);

        if (enablePolymerMolarWeight) {
            plyvmhCoefficients_.resize(numRegions);
        }
    }

    /*!
     * \brief Specify the maximum polymer concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    static void setPlymax(unsigned mixRegionIdx,
                          const Scalar& plymaxMaxConcentration)
    {
        plymaxMaxConcentration_[mixRegionIdx] = plymaxMaxConcentration;
    }

    /*!
     * \brief Specify the maximum polymer concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    static void setPlmixpar(unsigned mixRegionIdx,
                            const Scalar& plymixparToddLongstaff)
    {
        plymixparToddLongstaff_[mixRegionIdx] = plymixparToddLongstaff;
    }

    /*!
    * \brief get the PLYMWINJ table
    */
    static TabulatedTwoDFunction& getPlymwinjTable(const int tableNumber)
    {
        const auto iterTable = plymwinjTables_.find(tableNumber);
        if (iterTable != plymwinjTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the PLYMWINJ table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
    * \brief get the SKPRWAT table
    */
    static TabulatedTwoDFunction& getSkprwatTable(const int tableNumber)
    {
        const auto iterTable = skprwatTables_.find(tableNumber);
        if (iterTable != skprwatTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the SKPRWAT table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
    * \brief get the SKPRPOLY table
    */
    static SkprpolyTable& getSkprpolyTable(const int tableNumber)
    {
        const auto iterTable = skprpolyTables_.find(tableNumber);
        if (iterTable != skprpolyTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the SKPRPOLY table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
     * \brief Register all run-time parameters for the black-oil polymer module.
     */
    static void registerParameters()
    {
        if (!enablePolymer)
            // polymers have been disabled at compile time
            return;

        Opm::VtkBlackOilPolymerModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all polymer specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enablePolymer)
            // polymers have been disabled at compile time
            return;

        model.addOutputModule(new Opm::VtkBlackOilPolymerModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enablePolymer)
            // polymers have been disabled at compile time
            return false;

        if (!enablePolymerMolarWeight)
           return pvIdx == polymerConcentrationIdx;

        // both enablePolymer and enablePolymerMolarWeight are true here
        return pvIdx == polymerConcentrationIdx || pvIdx == polymerMoleWeightIdx;
    }

    static std::string primaryVarName(unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        if (pvIdx == polymerConcentrationIdx) {
            return "polymer_waterconcentration";
        }
        else {
            return "polymer_molecularweight";
        }
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enablePolymer)
            return false;

        if (!enablePolymerMolarWeight)
            return eqIdx == contiPolymerEqIdx;

        // both enablePolymer and enablePolymerMolarWeight are true here
        return eqIdx == contiPolymerEqIdx || eqIdx == contiPolymerMolarWeightEqIdx;
    }

    static std::string eqName(unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        if (eqIdx == contiPolymerEqIdx)
            return "conti^polymer";
        else
            return "conti^polymer_molecularweight";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants,
                           unsigned globalCellIdx)
    {
        if (!enablePolymer)
            return;

        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeWater =
                Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        surfaceVolumeWater = Opm::max(surfaceVolumeWater, 1e-10);

        // polymer in water phase
        const LhsEval massPolymer = surfaceVolumeWater
                * Toolbox::template decay<LhsEval>(intQuants.polymerConcentration())
                * (1.0 - Toolbox::template decay<LhsEval>(intQuants.polymerDeadPoreVolume()));

        // polymer in solid phase
        const LhsEval adsorptionPolymer =
                Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity())
                * Toolbox::template decay<LhsEval>(intQuants.polymerRockDensity())
                * Toolbox::template decay<LhsEval>(intQuants.polymerAdsorption());

        LhsEval accumulationPolymer = massPolymer + adsorptionPolymer;

        storage[contiPolymerEqIdx] += accumulationPolymer;

        // tracking the polymer molecular weight
        if (enablePolymerMolarWeight) {
            accumulationPolymer = Opm::max(accumulationPolymer, 1e-10);

            storage[contiPolymerMolarWeightEqIdx]  += accumulationPolymer
                                         * Toolbox::template decay<LhsEval> (intQuants.polymerMoleWeight());
        }
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enablePolymer)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const unsigned upIdx = extQuants.upstreamIndex(FluidSystem::waterPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const unsigned contiWaterEqIdx = Indices::conti0EqIdx + Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);

        unsigned globalSpaceIdxUp = elemCtx.globalSpaceIndex(upIdx, timeIdx);
        unsigned globalSpaceIdxInt = elemCtx.globalSpaceIndex(inIdx, timeIdx);

        if (upIdx == inIdx) {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *up.fluidState().invB(waterPhaseIdx)
                    *up.polymerViscosityCorrection()
                    /extQuants.polymerShearFactor()
                    *up.polymerConcentration();

            // modify water
            flux[contiWaterEqIdx] /=
                    extQuants.waterShearFactor();
        }
        else {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *Opm::decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                    *Opm::decay<Scalar>(up.polymerViscosityCorrection())
                    /Opm::decay<Scalar>(extQuants.polymerShearFactor())
                    *Opm::decay<Scalar>(up.polymerConcentration());

            // modify water
            flux[contiWaterEqIdx] /=
                    Opm::decay<Scalar>(extQuants.waterShearFactor());
        }

        // flux related to transport of polymer molecular weight
        if (enablePolymerMolarWeight) {
            if (upIdx == inIdx)
                flux[contiPolymerMolarWeightEqIdx] =
                    flux[contiPolymerEqIdx]*up.polymerMoleWeight();
            else
                flux[contiPolymerMolarWeightEqIdx] =
                    flux[contiPolymerEqIdx]*Opm::decay<Scalar>(up.polymerMoleWeight());
        }

    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the cange of polymer primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enablePolymer)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[polymerConcentrationIdx];
        outstream << priVars[polymerMoleWeightIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enablePolymer)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[polymerConcentrationIdx];
        instream >> priVars0[polymerMoleWeightIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1[polymerConcentrationIdx] = priVars0[polymerConcentrationIdx];
        priVars1[polymerMoleWeightIdx] = priVars0[polymerMoleWeightIdx];
    }

    static const Scalar plyrockDeadPoreVolume(const ElementContext& elemCtx,
                                              unsigned scvIdx,
                                              unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockDeadPoreVolume_[satnumRegionIdx];
    }

    static const Scalar plyrockResidualResistanceFactor(const ElementContext& elemCtx,
                                                        unsigned scvIdx,
                                                        unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockResidualResistanceFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockRockDensityFactor(const ElementContext& elemCtx,
                                                 unsigned scvIdx,
                                                 unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockRockDensityFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockAdsorbtionIndex(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockAdsorbtionIndex_[satnumRegionIdx];
    }

    static const Scalar plyrockMaxAdsorbtion(const ElementContext& elemCtx,
                                             unsigned scvIdx,
                                             unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockMaxAdsorbtion_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyadsAdsorbedPolymer(const ElementContext& elemCtx,
                                                          unsigned scvIdx,
                                                          unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyadsAdsorbedPolymer_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyviscViscosityMultiplierTable(const ElementContext& elemCtx,
                                                                    unsigned scvIdx,
                                                                    unsigned timeIdx)
    {
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const TabulatedFunction& plyviscViscosityMultiplierTable(unsigned pvtnumRegionIdx)
    {
        return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const Scalar plymaxMaxConcentration(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymaxMaxConcentration_[polymerMixRegionIdx];
    }

    static const Scalar plymixparToddLongstaff(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymixparToddLongstaff_[polymerMixRegionIdx];
    }

    static const PlyvmhCoefficients& plyvmhCoefficients(const ElementContext& elemCtx,
                                                        const unsigned scvIdx,
                                                        const unsigned timeIdx)
    {
        const unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyvmhCoefficients_[polymerMixRegionIdx];
    }

    static bool hasPlyshlog()
    {
        return hasPlyshlog_;
    }

    static bool hasShrate()
    {
        return hasShrate_;
    }

    static const Scalar shrate(unsigned pvtnumRegionIdx)
    {
        return shrate_[pvtnumRegionIdx];
    }

    /*!
     * \brief Computes the shear factor
     *
     * Input is polymer concentration and either the water velocity or the shrate if hasShrate_ is true.
     * The pvtnumRegionIdx is needed to make sure the right table is used.
     */
    template <class Evaluation>
    static Evaluation computeShearFactor(const Evaluation& polymerConcentration,
                                         unsigned pvtnumRegionIdx,
                                         const Evaluation& v0)
    {
        using ToolboxLocal = Opm::MathToolbox<Evaluation>;

        const auto& viscosityMultiplierTable = plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
        Scalar viscosityMultiplier = viscosityMultiplierTable.eval(Opm::scalarValue(polymerConcentration), /*extrapolate=*/true);

        const Scalar eps = 1e-14;
        // return 1.0 if the polymer has no effect on the water.
        if (std::abs((viscosityMultiplier - 1.0)) < eps)
            return ToolboxLocal::createConstant(v0, 1.0);

        const std::vector<Scalar>& shearEffectRefLogVelocity = plyshlogShearEffectRefLogVelocity_[pvtnumRegionIdx];
        auto v0AbsLog = Opm::log(Opm::abs(v0));
        // return 1.0 if the velocity /sharte is smaller than the first velocity entry.
        if (v0AbsLog < shearEffectRefLogVelocity[0])
            return ToolboxLocal::createConstant(v0, 1.0);

        // compute shear factor from input
        // Z = (1 + (P - 1) * M(v)) / P
        // where M(v) is computed from user input
        // and P = viscosityMultiplier
        const std::vector<Scalar>& shearEffectRefMultiplier = plyshlogShearEffectRefMultiplier_[pvtnumRegionIdx];
        size_t numTableEntries = shearEffectRefLogVelocity.size();
        assert(shearEffectRefMultiplier.size() == numTableEntries);

        std::vector<Scalar> shearEffectMultiplier(numTableEntries, 1.0);
        for (size_t i = 0; i < numTableEntries; ++i) {
            shearEffectMultiplier[i] = (1.0 + (viscosityMultiplier - 1.0)*shearEffectRefMultiplier[i]) / viscosityMultiplier;
            shearEffectMultiplier[i] = Opm::log(shearEffectMultiplier[i]);
        }
        // store the logarithmic velocity and logarithmic multipliers in a table for easy look up and
        // linear interpolation in the logarithmic space.
        TabulatedFunction logShearEffectMultiplier = TabulatedFunction(numTableEntries, shearEffectRefLogVelocity, shearEffectMultiplier, /*bool sortInputs =*/ false);

        // Find sheared velocity (v) that satisfies
        // F = log(v) + log (Z) - log(v0) = 0;

        // Set up the function
        // u = log(v)
        auto F = [&logShearEffectMultiplier, &v0AbsLog](const Evaluation& u) {
            return u + logShearEffectMultiplier.eval(u, true) - v0AbsLog;
        };
        // and its derivative
        auto dF = [&logShearEffectMultiplier](const Evaluation& u) {
            return 1 + logShearEffectMultiplier.evalDerivative(u, true);
        };

        // Solve F = 0 using Newton
        // Use log(v0) as initial value for u
        auto u = v0AbsLog;
        bool converged = false;
        // TODO make this into parameters
        for (int i = 0; i < 20; ++i) {
            auto f = F(u);
            auto df = dF(u);
            u -= f/df;
            if (std::abs(Opm::scalarValue(f)) < 1e-12) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error("Not able to compute shear velocity. \n");
        }

        // return the shear factor
        return Opm::exp(logShearEffectMultiplier.eval(u, /*extrapolate=*/true));

    }
    
    // TO DO: How many, and which, of the input arguments should be Evaluations ? (EvalWell)
    template <class Evaluation>
    static Evaluation computeWellShearRate(const Evaluation& darcyVelocityWell,
                                           const Scalar absolutePermeability,
                                           const Scalar porosity,
                                           const Evaluation& Sw,
                                           const Evaluation& relWater,
                                           const Scalar alphaFactor)
    {
        return  (4.0*alphaFactor*darcyVelocityWell)/(Opm::sqrt(8.0*relWater*absolutePermeability*porosity*Sw));
    }
    
    // TO DO: How many, and which, of the input arguments should be Evaluations ?
    template <class Evaluation>
    static Evaluation calcPolymerApparentViscosityForMechanicalDegradationModel(const Evaluation& polymerConcentration,
                                                                                const Evaluation& Mw,
                                                                                unsigned polymerMixRegionIdx,
                                                                                const Evaluation& shearRate,
                                                                                const Evaluation& muWater,
                                                                                const Scalar& porosity,
                                                                                const Scalar& temperature,
                                                                                const Evaluation& Sw,
                                                                                const Evaluation& relWater,
                                                                                bool elongation=false)
    {
        const bool isPolymer = (polymerConcentration > 1.0e-12) ? true : false;
        if(!isPolymer){
            return Evaluation(muWater);
        }
        const auto& plyvmhCoefficients = plyvmhCoefficients_[polymerMixRegionIdx];
        const Scalar k_mh = plyvmhCoefficients.k_mh;
        const Scalar a_mh = plyvmhCoefficients.a_mh;
        
        // hack: currently this must be manually updated in TWO different places..!!!
        const bool useCarreauYasudaModel = true;
        
        Scalar huggins1 = 0.4;
        Scalar huggins2 = 0.1;
        if(useCarreauYasudaModel){
            huggins1 = 0.01;
            huggins2 = 0.077;
        }
        const Scalar pow_x = 1.0;
        const Evaluation intrinsicViscosity(k_mh * Opm::pow(Mw, a_mh));
        const Evaluation x(polymerConcentration * intrinsicViscosity);
        const Evaluation specificViscosityAtZeroShearRate(x*(1 + x*(huggins1 + huggins2*x))); // x + k1*x^2 + k2*x^3
        
        const Scalar lambda_a = 0.350;
        const Scalar a_n = 0.063;
        const Scalar b_n = 0.804;
        
        // Calculate Carreau-Yasuda parameters:
        
        const Evaluation pow_n = 1.0 - Opm::pow(1.0 + Opm::pow(a_n * x, b_n), -1.0); // NOTE: pow([], -1) == 1/[]
        //const Evaluation pow_n(1.0 - 1.0/(1.0+pow(a_n*x, b_n)));
        
        // note: lam1_ICS (Lohne et al. 2017) = (lam1_OPM * 1.07)*(1.0+2.5*0.0015) when using muWater=1.0 and cpol=1.5 kg/m3
        const Evaluation lambda1((lambda_a*muWater*specificViscosityAtZeroShearRate*Mw)/(polymerConcentration*temperature));
        
        const Scalar gamma_hf = 20.0; // 1/sec
        const Scalar p_alpha = 1.75;

        // hack. Should really set to zero here, but need it to work for both Evaluation & EvalWell...
        Evaluation shearThinningSpecificViscosity = intrinsicViscosity;
        if(useCarreauYasudaModel){
            shearThinningSpecificViscosity = specificViscosityAtZeroShearRate*Opm::pow(1+Opm::pow(lambda1*shearRate, pow_x), -pow_n/pow_x);
        }
        else{
            // Meter's model:
            shearThinningSpecificViscosity = specificViscosityAtZeroShearRate / (1.0 + Opm::pow(shearRate/gamma_hf, p_alpha - 1.0));
        }
        Evaluation specificViscosity = shearThinningSpecificViscosity;
        
        if(elongation){
            const Scalar critical_Deborah_number = 0.5; // el_lamf = 2.0 / N_De
            const Scalar el_m2 = 1.5;
            const Scalar el_x2 = 3.0;
            // Used in som old code: R= 8.314459848
            const Scalar gas_constant = 8.31446261815324; // J/K/mol
            const Scalar pre_factor = 3.0/(5.0*gas_constant*critical_Deborah_number);
            // Note (9/5-2019): pre_factor includes effect of oil being present on lambda2.
            Evaluation fac_k_Sw = relWater;
            if(Sw > 0){
                fac_k_Sw /= Sw;
            }
            // TO DO: Add ipv and/or further Sw correction here??
            const Evaluation volume_factor(Opm::sqrt(fac_k_Sw)* porosity / (1. - porosity));
            // note: lam2_ICS (Lohne et al. 2017) = (lam2_OPM * 1.07)*(1.0+2.5*0.0015) when using muWater=1.0 and cpol=1.5 kg/m3
            const Evaluation lambda2(pre_factor*volume_factor*muWater*intrinsicViscosity*Mw/temperature);
            const Evaluation t_exp((el_m2 + pow_n) /  el_x2);
            const Evaluation elongationFactor(Opm::pow(1.0 + Opm::pow(lambda2 * shearRate, el_x2), t_exp));
            specificViscosity *= elongationFactor;
        }
        
        Evaluation relativeViscosity(1.0 + specificViscosity);
        return muWater*relativeViscosity;
        
    }

    const Scalar molarMass() const
    {
        return 0.25; // kg/mol
    }
    

private:
    static std::vector<Scalar> plyrockDeadPoreVolume_;
    static std::vector<Scalar> plyrockResidualResistanceFactor_;
    static std::vector<Scalar> plyrockRockDensityFactor_;
    static std::vector<Scalar> plyrockAdsorbtionIndex_;
    static std::vector<Scalar> plyrockMaxAdsorbtion_;
    static std::vector<TabulatedFunction> plyadsAdsorbedPolymer_;
    static std::vector<TabulatedFunction> plyviscViscosityMultiplierTable_;
    static std::vector<Scalar> plymaxMaxConcentration_;
    static std::vector<Scalar> plymixparToddLongstaff_;
    static std::vector<std::vector<Scalar>> plyshlogShearEffectRefMultiplier_;
    static std::vector<std::vector<Scalar>> plyshlogShearEffectRefLogVelocity_;
    static std::vector<Scalar> shrate_;
    static bool hasShrate_;
    static bool hasPlyshlog_;

    static std::vector<PlyvmhCoefficients> plyvmhCoefficients_;
    static std::map<int, TabulatedTwoDFunction> plymwinjTables_;
    static std::map<int, TabulatedTwoDFunction> skprwatTables_;

    static std::map<int, SkprpolyTable> skprpolyTables_;
};



template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockDeadPoreVolume_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockResidualResistanceFactor_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockRockDensityFactor_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockAdsorbtionIndex_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockMaxAdsorbtion_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyadsAdsorbedPolymer_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyviscViscosityMultiplierTable_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plymaxMaxConcentration_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plymixparToddLongstaff_;

template <class TypeTag, bool enablePolymerV>
std::vector<std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyshlogShearEffectRefMultiplier_;

template <class TypeTag, bool enablePolymerV>
std::vector<std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyshlogShearEffectRefLogVelocity_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::shrate_;

template <class TypeTag, bool enablePolymerV>
bool
BlackOilPolymerModule<TypeTag, enablePolymerV>::hasShrate_;

template <class TypeTag, bool enablePolymerV>
bool
BlackOilPolymerModule<TypeTag, enablePolymerV>::hasPlyshlog_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::PlyvmhCoefficients>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyvmhCoefficients_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedTwoDFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plymwinjTables_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedTwoDFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::skprwatTables_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::SkprpolyTable>
BlackOilPolymerModule<TypeTag, enablePolymerV>::skprpolyTables_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::StashedValues>
BlackOilPolymerModule<TypeTag, enablePolymerV>::stashedValues_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilPolymerIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag, bool enablePolymerV = getPropValue<TypeTag, Properties::EnablePolymer>()>
class BlackOilPolymerIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using StashedValues = typename PolymerModule::StashedValues;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr bool enablePolymerMechanicalDegradation = getPropValue<TypeTag, Properties::EnablePolymerMechanicalDegradation>();
    static constexpr int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;


public:

    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    void polymerPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        
        emvis_ = 17.0; // question: what is this quantity ???
        epvis_ = 1.0;
        ewv_pol_ = 1.0;
        eta_rel_ = 1.0;
        shear_rate_ = 0.0;
        
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        polymerConcentration_ = priVars.makeEvaluation(polymerConcentrationIdx, timeIdx, linearizationType);
        if (enablePolymerMolarWeight) {
            polymerMoleWeight_ = priVars.makeEvaluation(polymerMoleWeightIdx, timeIdx, linearizationType);
        }
        const Scalar cmax = PolymerModule::plymaxMaxConcentration(elemCtx, dofIdx, timeIdx);

        // permeability reduction due to polymer
        const Scalar& maxAdsorbtion = PolymerModule::plyrockMaxAdsorbtion(elemCtx, dofIdx, timeIdx);
        const auto& plyadsAdsorbedPolymer = PolymerModule::plyadsAdsorbedPolymer(elemCtx, dofIdx, timeIdx);
        polymerAdsorption_ = plyadsAdsorbedPolymer.eval(polymerConcentration_, /*extrapolate=*/true);
        if (PolymerModule::plyrockAdsorbtionIndex(elemCtx, dofIdx, timeIdx) == PolymerModule::NoDesorption) {
            const Scalar& maxPolymerAdsorption = elemCtx.problem().maxPolymerAdsorption(elemCtx, dofIdx, timeIdx);
            polymerAdsorption_ = std::max(Evaluation(maxPolymerAdsorption) , polymerAdsorption_);
        }

        // compute resitanceFactor
        const Scalar& residualResistanceFactor = PolymerModule::plyrockResidualResistanceFactor(elemCtx, dofIdx, timeIdx);
        const Evaluation resistanceFactor = 1.0 + (residualResistanceFactor - 1.0) * polymerAdsorption_ / maxAdsorbtion;

        // compute effective viscosities
        if (!enablePolymerMolarWeight) {
            const auto& fs = asImp_().fluidState_;
            const Evaluation& muWater = fs.viscosity(waterPhaseIdx);
            const auto& viscosityMultiplier = PolymerModule::plyviscViscosityMultiplierTable(elemCtx, dofIdx, timeIdx);
            const Evaluation viscosityMixture = viscosityMultiplier.eval(polymerConcentration_, /*extrapolate=*/true) * muWater;

            // Do the Todd-Longstaff mixing
            const Scalar plymixparToddLongstaff = PolymerModule::plymixparToddLongstaff(elemCtx, dofIdx, timeIdx);
            const Evaluation viscosityPolymer = viscosityMultiplier.eval(cmax, /*extrapolate=*/true) * muWater;
            const Evaluation viscosityPolymerEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(viscosityPolymer, 1.0 - plymixparToddLongstaff);
            const Evaluation viscosityWaterEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(muWater, 1.0 - plymixparToddLongstaff);

            const Evaluation cbar = polymerConcentration_ / cmax;
            // waterViscosity / effectiveWaterViscosity
            waterViscosityCorrection_ = muWater * ((1.0 - cbar) / viscosityWaterEffective + cbar / viscosityPolymerEffective);
            // effectiveWaterViscosity / effectivePolymerViscosity
            polymerViscosityCorrection_ =  (muWater / waterViscosityCorrection_) / viscosityPolymerEffective;
            
            epvis_ = viscosityPolymerEffective;
            ewv_pol_ = viscosityWaterEffective;
            eta_rel_ = 1.0; // not used together with eclipse model
            
            /*  NOTES (odn):
            *
            *   viscosityMixture = mu_m(C_p) (== viscosity of fully mixed polymer solution)
            *   viscosityPolymerEffective = mu_p,eff
            *   viscosityWaterEffective = mu_w,e (== partially mixed water viscosity)
            *
            *   waterViscosityCorrection_ = [(1-C_bar)/mu_w,e + C_bar/mu_p,eff]*muWater = mu_water / mu_w,eff
            *
            *       ---> mult. with 1/mu_water to get 1/mu_w,eff in the water flow equation
            *            [adjust water phase mobility inside polymerPropertiesUpdate_]
            *
            *   polymerViscosityCorrection_ = ((muWater)*(mu_w,eff/mu_water)) / mu_p,eff = mu_w,eff / mu_p,eff
            *
            *       ---> First, waterVolumeFlux_, which already includes waterViscosityCorrection_ (??), 
            *            is recovered from extQuants.volumeFlux(waterPhaseIdx).
            *       ---> Next, the polymer flux is computed from waterVolumeFlux_ by further multiplying with 
            *            polymerViscosityCorrection_. 
            *            This is done for both the regular polymer component and the Mw component.
            *       ---> Since the water flux has already been corrected by mu_water / mu_w,eff, the net effect 
            *            is to correct by mu_water / mu_p,eff in the polymer equations.
            */
            
        }
        else if(enablePolymerMechanicalDegradation) // based on hard-coded parameters + PLYVMH (for now)
        {
            const auto& fs = asImp_().fluidState_;
            
            const Evaluation relWater = 1.0; // TO DO: Get rel.perm. of water also?
            const Evaluation& muWater = fs.viscosity(waterPhaseIdx);
            const Evaluation& Sw = fs.saturation(waterPhaseIdx);
            /*
            const Scalar porosity = elemCtx.problem().porosity(elemCtx, dofIdx, timeIdx);
            const Scalar temperature = elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx);
            
            const auto& plyvmhCoefficients = PolymerModule::plyvmhCoefficients(elemCtx, dofIdx, timeIdx);
            const Scalar k_mh = plyvmhCoefficients.k_mh;
            const Scalar a_mh = plyvmhCoefficients.a_mh;
            const Scalar gamma = plyvmhCoefficients.gamma;
            const Scalar kappa = plyvmhCoefficients.kappa;

            Scalar shear_rate = 0.0;
            unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
            if(!PolymerModule::stashedValues_.empty()){
                StashedValues& stashedValues = PolymerModule::stashedValues_[globalSpaceIdx];
                unsigned cellGlobalIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
                shear_rate = stashedValues.cell_shear_rate;
                
                // TO DO: just for testing....
                epvis_ = stashedValues.debug_value;
                
            }
            
            const bool isPolymer = (polymerConcentration_.value() > 1.0e-12) ? true : false;
            const Scalar cPol = isPolymer ? polymerConcentration_.value() : 0.0;
            
            // hack: currently this must be manually updated in TWO different places..!!!
            const bool useCarreauYasudaModel = true;
            
            Scalar huggins1 = 0.4;
            Scalar huggins2 = 0.1;
            if(useCarreauYasudaModel){
                huggins1 = 0.01;
                huggins2 = 0.077;
            }
            const Scalar pow_x = 1.0;
            
            const Evaluation Mw = polymerMoleWeight_ * 1000.;
            const Evaluation intrinsicViscosity = k_mh * pow(Mw, a_mh);
            const Evaluation x = cPol * intrinsicViscosity;
            const Evaluation specificViscosityAtZeroShearRate = x*(1 + x*(huggins1 + huggins2*x)); // x + k1*x^2 + k2*x^3
            
            const Scalar lambda_a = 0.350;
            const Scalar a_n = 0.063;
            const Scalar b_n = 0.804;
            
            // Calculate Carreau-Yasuda parameters:
            const Evaluation pow_n = 1.0 - 1.0/(1.0+pow(a_n*x, b_n));
            
            // note: lam1_ICS (Lohne et al. 2017) = (lam1_OPM * 1.07)*(1.0+2.5*0.0015) when using muWater=1.0 and cpol=1.5 kg/m3
            Evaluation lambda1 = lambda_a*muWater*specificViscosityAtZeroShearRate*Mw; // unit of mu = Pa*s ??
            if(cPol > 0.0){
                const Scalar denominator = cPol*temperature;
                lambda1 /= denominator;
            }
            
            const Scalar gamma_hf = 20.0; // 1/sec
            const Scalar p_alpha = 1.75;
            
            Evaluation shearThinningSpecificViscosity = 0.0;
            if(useCarreauYasudaModel){
                shearThinningSpecificViscosity = specificViscosityAtZeroShearRate*pow(1+pow(lambda1*shear_rate, pow_x), -pow_n/pow_x);
            }
            else{
                // Meter's model:
                shearThinningSpecificViscosity = specificViscosityAtZeroShearRate / (1.0 + pow(shear_rate/gamma_hf, p_alpha - 1.0));
            }
        
            // Elongation & degradation parameters:
            const Scalar critical_Deborah_number = 0.5; // el_lamf = 2.0 / N_De
            const Scalar el_m2 = 1.5;
            const Scalar el_x2 = 3.0;
            
            const Scalar gas_constant = 8.31446261815324; // J/K/mol
            const Scalar pre_factor = 3.0/(5.0*gas_constant*critical_Deborah_number);
            // Note (9/5-2019): pre_factor includes effect of oil being present on lambda2.
            const Evaluation fac_k_Sw = (Sw > 0) ? (relWater / Sw) : relWater;
            // TO DO: Add ipv and/or further Sw correction here??
            const Evaluation volume_factor = Opm::sqrt(fac_k_Sw)* porosity / (1. - porosity);
            // note: lam2_ICS (Lohne et al. 2017) = (lam2_OPM * 1.07)*(1.0+2.5*0.0015) when using muWater=1.0 and cpol=1.5 kg/m3
            const Evaluation lambda2 = pre_factor*volume_factor*muWater*intrinsicViscosity*Mw/temperature;
            const Evaluation t_exp = (el_m2 + pow_n) / el_x2;
            const Evaluation elongationFactor = Opm::pow(1.0 + Opm::pow(lambda2 * shear_rate, el_x2), t_exp);
            
            // TO DO: when to select shear thickening??
            const Evaluation specificViscosity = shearThinningSpecificViscosity; // *elongationFactor;
            const Evaluation relativeViscosity = 1.0 + specificViscosity;
            
            waterViscosityCorrection_ = 1.0 / relativeViscosity;
            polymerViscosityCorrection_ =  1.0;
            
            // TO DO: Not all quantities printed to summary file match their names (FOR DEBUGGING)
            emvis_ = lambda2;
            epvis_ = shearThinningSpecificViscosity*elongationFactor;
            eta_rel_ = relativeViscosity;
            shear_rate_ = shear_rate;
            */
            
            // DEBUG ONLY:
            Evaluation shear_rate = 0.0;
            unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
            if(!PolymerModule::stashedValues_.empty()){
                StashedValues& stashedValues = PolymerModule::stashedValues_[globalSpaceIdx];
                unsigned cellGlobalIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
                shear_rate = stashedValues.cell_shear_rate;
                
                // TO DO: just for testing....
                epvis_ = stashedValues.debug_value;
            }
            
            const Scalar porosity = elemCtx.problem().porosity(elemCtx, dofIdx, timeIdx);
            const Scalar temperature = elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx);
            const unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, dofIdx, timeIdx);
            const Evaluation& apparentViscosity = PolymerModule::calcPolymerApparentViscosityForMechanicalDegradationModel(polymerConcentration_,
                                                                                                                           polymerMoleWeight_*1.0e3,  // --> kg/mol
                                                                                                                           polymerMixRegionIdx,
                                                                                                                           shear_rate,
                                                                                                                           muWater,
                                                                                                                           porosity,
                                                                                                                           temperature,
                                                                                                                           Sw,
                                                                                                                           relWater);
                                                                                                                           
            const Evaluation relativeViscosity = apparentViscosity / muWater;                           
            waterViscosityCorrection_ = 1.0 / relativeViscosity;
            polymerViscosityCorrection_ = 1.0;

            //emvis_ = lambda2;
            //epvis_ = shearThinningSpecificViscosity*elongationFactor;
            eta_rel_ = relativeViscosity;
            shear_rate_ = shear_rate;
            //std::cout << "Viscosity=" << 1.0e3*visc << " mPas\n.";
        }
        else { // based on PLYVMH
            const auto& plyvmhCoefficients = PolymerModule::plyvmhCoefficients(elemCtx, dofIdx, timeIdx);
            const Scalar k_mh = plyvmhCoefficients.k_mh;
            const Scalar a_mh = plyvmhCoefficients.a_mh;
            const Scalar gamma = plyvmhCoefficients.gamma;
            const Scalar kappa = plyvmhCoefficients.kappa;

            // viscosity model based on Mark-Houwink equation and Huggins equation
            // 1000 is a emperical constant, most likely related to unit conversion
            const Evaluation intrinsicViscosity = k_mh * pow(polymerMoleWeight_ * 1000., a_mh);
            const Evaluation x = polymerConcentration_ * intrinsicViscosity;
            const Evaluation relativeViscosity = 1.0 + gamma * (x + kappa * x * x);
            
            //waterViscosityCorrection_ = 1.0 / (1.0 + gamma * (x + kappa * x * x));
            waterViscosityCorrection_ = 1.0 / relativeViscosity;
            polymerViscosityCorrection_ = 1.0;
            
            // TO DO: Are EPVIS, EWVPOL, EMVIS relevant here??
            eta_rel_ = relativeViscosity;
        }

        // adjust water mobility
        asImp_().mobility_[waterPhaseIdx] *= waterViscosityCorrection_ / resistanceFactor;

        // update rock properties
        polymerDeadPoreVolume_ = PolymerModule::plyrockDeadPoreVolume(elemCtx, dofIdx, timeIdx);
        polymerRockDensity_ = PolymerModule::plyrockRockDensityFactor(elemCtx, dofIdx, timeIdx);
    }

    const Evaluation& polymerConcentration() const
    { return polymerConcentration_; }

    const Evaluation& polymerMoleWeight() const
    {
        if (!enablePolymerMolarWeight)
            throw std::logic_error("polymerMoleWeight() is called but polymer milecular weight is disabled");

        return polymerMoleWeight_;
    }

    const Scalar& polymerDeadPoreVolume() const
    { return polymerDeadPoreVolume_; }

    const Evaluation& polymerAdsorption() const
    { return polymerAdsorption_; }

    const Scalar& polymerRockDensity() const
    { return polymerRockDensity_; }

    // effectiveWaterViscosity / effectivePolymerViscosity
    const Evaluation& polymerViscosityCorrection() const
    { return polymerViscosityCorrection_; }

    // waterViscosity / effectiveWaterViscosity
    const Evaluation& waterViscosityCorrection() const
    { return waterViscosityCorrection_; }
    
    const Evaluation& polymerRelativeViscosity() const
    { 
    if (!enablePolymerMolarWeight)
        throw std::logic_error("polymerRelativeViscosity() is called but polymer molecular weight is disabled");
    
    return eta_rel_;
    }
    
    const Evaluation& epvis() const
    { 
    if (!enablePolymerMolarWeight)
        throw std::logic_error("epvis() is called but polymer molecular weight is disabled");
    
    return epvis_; 
    }
    
    const Evaluation& emvis() const
    { 
    if (!enablePolymerMolarWeight)
        throw std::logic_error("emvis() is called but polymer molecular weight is disabled");
    
    return emvis_; 
    }
    
    const Evaluation& ewvpol() const
    { 
    if (!enablePolymerMolarWeight)
        throw std::logic_error("ewvpol() called but but polymer molecular weight is disabled"); 
    
    return ewv_pol_;
    }
    
    const Evaluation& shear() const
    { 
    if (!enablePolymerMolarWeight)
        throw std::logic_error("shear() called but but polymer molecular weight is disabled"); 
    
    return shear_rate_;
    }


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerConcentration_;
    Evaluation polymerMoleWeight_;
    Scalar polymerDeadPoreVolume_;
    Scalar polymerRockDensity_;
    Evaluation polymerAdsorption_;
    Evaluation polymerViscosityCorrection_;
    Evaluation waterViscosityCorrection_;
    
    Evaluation eta_rel_;
    Evaluation emvis_;
    Evaluation epvis_;
    Evaluation ewv_pol_;
    Evaluation shear_rate_;

};

template <class TypeTag>
class BlackOilPolymerIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:

    void polymerPropertiesUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned scvIdx OPM_UNUSED,
                                  unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& polymerMoleWeight() const
    { throw std::logic_error("polymerMoleWeight() called but polymer molecular weight is disabled"); }

    const Evaluation& polymerConcentration() const
    { throw std::runtime_error("polymerConcentration() called but polymers are disabled"); }

    const Evaluation& polymerDeadPoreVolume() const
    { throw std::runtime_error("polymerDeadPoreVolume() called but polymers are disabled"); }

    const Evaluation& polymerAdsorption() const
    { throw std::runtime_error("polymerAdsorption() called but polymers are disabled"); }

    const Evaluation& polymerRockDensity() const
    { throw std::runtime_error("polymerRockDensity() called but polymers are disabled"); }

    const Evaluation& polymerViscosityCorrection() const
    { throw std::runtime_error("polymerViscosityCorrection() called but polymers are disabled"); }

    const Evaluation& waterViscosityCorrection() const
    { throw std::runtime_error("waterViscosityCorrection() called but polymers are disabled"); }
    
    const Evaluation& polymerRelativeViscosity() const
    { throw std::runtime_error("polymerRelativeViscosity() called but polymers are disabled"); }
    
    const Evaluation& emvis() const
    { throw std::runtime_error("emvis() called but polymers are disabled"); }
    
    const Evaluation& epvis() const
    { throw std::runtime_error("epvis() called but polymers are disabled"); }
    
    const Evaluation& ewvpol() const
    { throw std::runtime_error("ewvpol() called but polymers are disabled"); }
    
    const Evaluation& shear() const
    { throw std::runtime_error("shear() called but polymers are disabled"); }
    
};


/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilPolymerExtensiveQuantities
 *
 * \brief Provides the polymer specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enablePolymerV = getPropValue<TypeTag, Properties::EnablePolymer>()>
class BlackOilPolymerExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    
    static constexpr bool enablePolymerMechanicalDegradation = getPropValue<TypeTag, Properties::EnablePolymerMechanicalDegradation>();

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr unsigned waterPhaseIdx =  FluidSystem::waterPhaseIdx;

    using Toolbox = Opm::MathToolbox<Evaluation>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;
    
    using StashedValues = typename PolymerModule::StashedValues;

public:
    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using permeabilities, i.e., *not* via transmissibilities.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateShearMultipliersPerm(const ElementContext& elemCtx OPM_UNUSED,
                                    unsigned scvfIdx OPM_UNUSED,
                                    unsigned timeIdx OPM_UNUSED)
    {
        throw std::runtime_error("The extension of the blackoil model for polymers is not yet "
                                 "implemented for problems specified using permeabilities.");
    }
    
    /*!
     * \brief Method which calculates in-situ shear rate based on flow velocity
     * 
     * Note that this implementation is separate from the ECLIPSE approach 
     * of using the SHRATE keyword and a corresponding shear rate factor.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
    // compiler errors if it is not instantiated!
    void calcInSituShearRate(const ElementContext& elemCtx,
                             unsigned scvfIdx,
                             unsigned timeIdx)
    {
        if(!enablePolymerMechanicalDegradation)
        {
            return;
        }
        
        // shear factors for water and polymer are not used for this option
        waterShearFactor_ = 1.0;
        polymerShearFactor_ = 1.0;

        const ExtensiveQuantities& extQuants = asImp_();
        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        unsigned interiorDofIdx = extQuants.interiorIndex();
        unsigned exteriorDofIdx = extQuants.exteriorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);
        
        // compute water velocity from flux
        Evaluation poroAvg = intQuantsIn.porosity()*0.5 + intQuantsEx.porosity()*0.5;
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvfIdx, timeIdx);
        const Evaluation& Sw = up.fluidState().saturation(waterPhaseIdx);
        unsigned cellIdx = elemCtx.globalSpaceIndex(scvfIdx, timeIdx);
        const auto& materialLawManager = elemCtx.problem().materialLawManager();
        const auto& scaledDrainageInfo =
                materialLawManager->oilWaterScaledEpsInfoDrainage(cellIdx);
        const Scalar& Swcr = scaledDrainageInfo.Swcr;
        
        // guard against zero porosity and no mobile water
        Evaluation denom = Opm::max(poroAvg * (Sw - Swcr), 1e-12); // NB!!!
        Evaluation waterVolumeVelocity = extQuants.volumeFlux(waterPhaseIdx) / denom;
        
        unsigned numInteriorFaces = elemCtx.numInteriorFaces(timeIdx);
        unsigned cellGlobalIdx = elemCtx.globalSpaceIndex(interiorDofIdx, timeIdx);
        
        StashedValues& stashedValues = PolymerModule::stashedValues_[cellGlobalIdx];
        if(scvfIdx == 0) // first face
        {
            stashedValues.velocity = 0.0;
            //stashedValues.distance = 0.0;
            stashedValues.average_cell_velocity = 0.0;
            stashedValues.cell_shear_rate = 0.0;
            stashedValues.debug_value = 0;
        }
        stashedValues.velocity[scvfIdx] = waterVolumeVelocity.value();

        // Compute shear rate from velocity:
        inSituShearRate_ = 0.0;
        
        // CRITICAL TODO: Rel. perm. is currently < 1 even when Sw=1...
        const Evaluation relWater = 1.0; //up.relativePermeability(waterPhaseIdx);
        
        // hack: currently need to (manually) edit this in TWO places....
        const bool useCarreauYasudaModel = true;
        const Scalar alpha_factor = useCarreauYasudaModel? 2.0 : 2.5;
        
        /*
        // THIS IS NOT USED....
        Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        if (trans > 0.0) {
             // Question: What if trans == 0 ?
            Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
            auto dist = elemCtx.pos(interiorDofIdx, timeIdx) -  elemCtx.pos(exteriorDofIdx, timeIdx);
            // compute permeability from transmissibility.
            const Scalar absPerm = trans / faceArea * dist.two_norm();
            const Scalar absWaterVolumeVelocity = Opm::abs(waterVolumeVelocity.value());
            inSituShearRate_ = absWaterVolumeVelocity*alpha_factor*Opm::sqrt(2.0*poroAvg*Sw / (relWater*absPerm));
            assert(Opm::isfinite(inSituShearRate_));
        }
        */
        
        // TO DOs (generalize to 2D & 3D): 
        // - How to weigh different velocities in the average?
        // - How to adjust for directional permeabilities??
        if (scvfIdx == numInteriorFaces - 1){ // last face for this interior cell
            Scalar sumVelSquared = 0.0;
            for(unsigned i=0; i < numInteriorFaces; ++i){
                //stashedValues.debug_value = static_cast<Scalar>(numInteriorFaces);//DEBUG ONLY
                Scalar velocity = stashedValues.velocity[i];
                sumVelSquared += velocity*velocity;
            }
            const Scalar avgCellVelocity = sqrt(sumVelSquared/numInteriorFaces); // m/sec
            
            //Note from source code: "This method is specific to the finite volume discretizations."
            const auto& K = elemCtx.problem().intrinsicPermeability(cellGlobalIdx); // 3x3 matrix
            const auto Kx = K[0][0];
            const auto Ky = K[1][1];
            const auto Kz = K[2][2];
            const Scalar absPermCell = Kx;
            Scalar cellShearRate = 0.0;
            Scalar the_rest = 0.0;
            if(absPermCell > 0 && relWater > 0){
                const Evaluation gamma = avgCellVelocity*alpha_factor*Opm::sqrt(2.0*poroAvg*Sw / (relWater*absPermCell));
                cellShearRate = gamma.value();
            }
            
            stashedValues.average_cell_velocity = avgCellVelocity;
            stashedValues.cell_shear_rate = cellShearRate;
            stashedValues.debug_value = relWater.value();//avgCellVelocity; // DEBUG ONLY
        }
    }

    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using transmissibilities, i.e., *not* via permeabilities.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
    // compiler errors if it is not instantiated!
    void updateShearMultipliers(const ElementContext& elemCtx,
                                unsigned scvfIdx,
                                unsigned timeIdx)
    {

        waterShearFactor_ = 1.0;
        polymerShearFactor_ = 1.0;

        if (!PolymerModule::hasPlyshlog()){
            return;
        }
        
        if(enablePolymerMechanicalDegradation){
            return;
        }

        const ExtensiveQuantities& extQuants = asImp_();
        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        unsigned interiorDofIdx = extQuants.interiorIndex();
        unsigned exteriorDofIdx = extQuants.exteriorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        // compute water velocity from flux
        Evaluation poroAvg = intQuantsIn.porosity()*0.5 + intQuantsEx.porosity()*0.5;
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvfIdx, timeIdx);
        const Evaluation& Sw = up.fluidState().saturation(waterPhaseIdx);
        unsigned cellIdx = elemCtx.globalSpaceIndex(scvfIdx, timeIdx);
        const auto& materialLawManager = elemCtx.problem().materialLawManager();
        const auto& scaledDrainageInfo =
                materialLawManager->oilWaterScaledEpsInfoDrainage(cellIdx);
        const Scalar& Swcr = scaledDrainageInfo.Swcr;

        // guard against zero porosity and no mobile water
        Evaluation denom = Opm::max(poroAvg * (Sw - Swcr), 1e-12);
        Evaluation waterVolumeVelocity = extQuants.volumeFlux(waterPhaseIdx) / denom;

        // if shrate is specified. Compute shrate based on the water velocity
        if (PolymerModule::hasShrate()) {
            const Evaluation& relWater = up.relativePermeability(waterPhaseIdx);
            Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
            if (trans > 0.0) {
                Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
                auto dist = elemCtx.pos(interiorDofIdx, timeIdx) -  elemCtx.pos(exteriorDofIdx, timeIdx);
                // compute permeability from transmissibility.
                Scalar absPerm = trans / faceArea * dist.two_norm();
                waterVolumeVelocity *=
                    PolymerModule::shrate(pvtnumRegionIdx)*Opm::sqrt(poroAvg*Sw / (relWater*absPerm));
                assert(Opm::isfinite(waterVolumeVelocity));
            }
        }

        // compute share factors for water and polymer
        waterShearFactor_ =
            PolymerModule::computeShearFactor(up.polymerConcentration(),
                                              pvtnumRegionIdx,
                                              waterVolumeVelocity);
        polymerShearFactor_ =
            PolymerModule::computeShearFactor(up.polymerConcentration(),
                                              pvtnumRegionIdx,
                                              waterVolumeVelocity*up.polymerViscosityCorrection());

    }

    const Evaluation& polymerShearFactor() const
    { return polymerShearFactor_; }

    const Evaluation& waterShearFactor() const
    { return waterShearFactor_; }
    
    const Evaluation& inSituShearRate() const
    { return inSituShearRate_; }


private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerShearFactor_;
    Evaluation waterShearFactor_;
    
    Evaluation inSituShearRate_;

};

template <class TypeTag>
class BlackOilPolymerExtensiveQuantities<TypeTag, false>
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    void updateShearMultipliers(const ElementContext& elemCtx OPM_UNUSED,
                                unsigned scvfIdx OPM_UNUSED,
                                unsigned timeIdx OPM_UNUSED)
    { }

    void updateShearMultipliersPerm(const ElementContext& elemCtx OPM_UNUSED,
                                    unsigned scvfIdx OPM_UNUSED,
                                    unsigned timeIdx OPM_UNUSED)
    { }
    
    void calcInSituShearRate(const ElementContext& elemCtx OPM_UNUSED,
                             unsigned scvfIdx OPM_UNUSED,
                             unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& polymerShearFactor() const
    { throw std::runtime_error("polymerShearFactor() called but polymers are disabled"); }

    const Evaluation& waterShearFactor() const
    { throw std::runtime_error("waterShearFactor() called but polymers are disabled"); }
    
    const Evaluation& inSituShearRate() const
    { throw std::runtime_error("inSituShearRate() called but polymers are disabled"); }
    
};


} // namespace Opm

#endif
