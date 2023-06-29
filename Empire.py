from __future__ import division
from pyomo.environ import *
from pyomo.common.tempfiles import TempfileManager
import csv
import sys
import cloudpickle
import time
import os
from datetime import datetime
from scipy.stats import norm

# import cartopy
# import cartopy.crs as ccrs
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D

def strfdelta(tdelta, fmt):
	d = {"days": tdelta.days}
	d["H"], rem = divmod(tdelta.seconds, 3600)
	d["H"] = str("{:02d}".format(d["H"]))
	d["M"], d["S"] = divmod(rem, 60)
	return fmt.format(**d)

# noinspection PyTypeChecker
def run_empire(name, tab_file_path, result_file_path, scenariogeneration, scenario_data_path,
			   solver, temp_dir, FirstHoursOfRegSeason, FirstHoursOfPeakSeason, lengthRegSeason,
			   lengthPeakSeason, Period, Operationalhour, Scenario, Season, HoursOfSeason,
			   NoOfNormalScenarios, NoOfHydrogenScenarios,
			   discountrate, WACC, LeapYearsInvestment, WRITE_LP,
			   PICKLE_INSTANCE, EMISSION_CAP, USE_TEMP_DIR, NoOfRegSeason, NoOfPeakSeason,
			   windfarmNodes = None, verboseResultWriting=False,
			   hydrogen=False, TIME_LIMIT=None,
			   h2storage=False, hydrogen_demand_percentage = 1.0, std_dev_percentage = 0):

	if USE_TEMP_DIR:
		TempfileManager.tempdir = temp_dir

	if not os.path.exists(result_file_path):
		os.makedirs(result_file_path)

	#Removing spaces from each element in offshoreNodesList, because the nodes have them removed too.
	if windfarmNodes is not None:
		for count,node in enumerate(windfarmNodes):
			windfarmNodes[count] = node.replace(' ', '')

	if (len(Scenario) is not NoOfNormalScenarios * NoOfHydrogenScenarios):
		print(f'Error: Number of scenarios ({len(Scenario)}) does not equal product of number of '
			  f'normal scenarios ({NoOfNormalScenarios}) and hydrogen scenarios ({NoOfHydrogenScenarios})!')
		exit()

	model = AbstractModel()

	###########
	##SOLVERS##
	###########

	if solver == "CPLEX":
		print("Solver: CPLEX")
	elif solver == "Xpress":
		print("Solver: Xpress")
	elif solver == "Gurobi":
		print("Solver: Gurobi")
	else:
		sys.exit("ERROR! Invalid solver! Options: CPLEX, Xpress, Gurobi")

	##########
	##MODULE##
	##########

	if WRITE_LP:
		print("Will write LP-file...")

	if PICKLE_INSTANCE:
		print("Will pickle instance...")

	if EMISSION_CAP:
		print("Absolute emission cap in each scenario...")
	else:
		print("No absolute emission cap...")

	if hydrogen is True:
		print("Optimizing with hydrogen component")

	########
	##SETS##
	########

	#Define the sets
	timeStart = datetime.now()
	print("Declaring sets...")

	#Supply technology sets
	model.Generator = Set(ordered=True) #g
	model.HydrogenGenerators = Set(ordered=True, within=model.Generator)
	model.Technology = Set(ordered=True) #t
	model.Storage =  Set() #b


	#Temporal sets
	model.Period = Set(ordered=True, initialize=Period) #i
	model.Operationalhour = Set(ordered=True, initialize=Operationalhour) #h
	model.Season = Set(ordered=True, initialize=Season) #s

	#Spatial sets
	model.Node = Set(ordered=True) #n
	model.DirectionalLink = Set(dimen=2, within=model.Node*model.Node, ordered=True) #a
	model.TransmissionType = Set(ordered=True)
	model.PipelineType = Set(ordered=True)

	if windfarmNodes is not None:
		#GD: Set of all offshore wind farm nodes. Need this set to restrict transmission through wind farms based on their invested capacity
		model.windfarmNodes = Set(ordered=True, within=model.Node, initialize=windfarmNodes)

	#Stochastic sets
	model.Scenario = Set(ordered=True, initialize=Scenario) #w

	#Subsets
	model.GeneratorsOfTechnology=Set(dimen=2) #(t,g) for all t in T, g in G_t
	model.GeneratorsOfNode = Set(dimen=2) #(n,g) for all n in N, g in G_n
	model.TransmissionTypeOfDirectionalLink = Set(dimen=3) #(n1,n2,t) for all (n1,n2) in L, t in T
	model.PipelineTypeOfDirectionalLink = Set(dimen=3) #(n1,n2,t) for all (n1,n2) in L, t in T
	model.ThermalGenerators = Set(within=model.Generator) #g_ramp
	model.RegHydroGenerator = Set(within=model.Generator) #g_reghyd
	model.HydroGenerator = Set(within=model.Generator) #g_hyd
	model.StoragesOfNode = Set(dimen=2) #(n,b) for all n in N, b in B_n
	model.DependentStorage = Set() #b_dagger
	model.HoursOfSeason = Set(dimen=2, ordered=True, initialize=HoursOfSeason) #(s,h) for all s in S, h in H_s
	model.FirstHoursOfRegSeason = Set(within=model.Operationalhour, ordered=True, initialize=FirstHoursOfRegSeason)
	model.FirstHoursOfPeakSeason = Set(within=model.Operationalhour, ordered=True, initialize=FirstHoursOfPeakSeason)


	print("Reading sets...")

	#Load the data

	data = DataPortal()
	data.load(filename=tab_file_path + "/" + 'Sets_Generator.tab',format="set", set=model.Generator)
	data.load(filename=tab_file_path + '/' + 'Sets_HydrogenGenerators.tab', format="set", set=model.HydrogenGenerators)
	data.load(filename=tab_file_path + "/" + 'Sets_ThermalGenerators.tab',format="set", set=model.ThermalGenerators)
	data.load(filename=tab_file_path + "/" + 'Sets_HydroGenerator.tab',format="set", set=model.HydroGenerator)
	data.load(filename=tab_file_path + "/" + 'Sets_HydroGeneratorWithReservoir.tab',format="set", set=model.RegHydroGenerator)
	data.load(filename=tab_file_path + "/" + 'Sets_Storage.tab',format="set", set=model.Storage)
	data.load(filename=tab_file_path + "/" + 'Sets_DependentStorage.tab',format="set", set=model.DependentStorage)
	data.load(filename=tab_file_path + "/" + 'Sets_Technology.tab',format="set", set=model.Technology)
	data.load(filename=tab_file_path + "/" + 'Sets_Node.tab',format="set", set=model.Node)
	data.load(filename=tab_file_path + "/" + 'Sets_DirectionalLines.tab',format="set", set=model.DirectionalLink)
	data.load(filename=tab_file_path + "/" + 'Sets_LineType.tab',format="set", set=model.TransmissionType)
	data.load(filename=tab_file_path + "/" + 'Sets_LineTypeOfDirectionalLines.tab',format="set", set=model.TransmissionTypeOfDirectionalLink)
	data.load(filename=tab_file_path + "/" + 'Sets_GeneratorsOfTechnology.tab',format="set", set=model.GeneratorsOfTechnology)
	data.load(filename=tab_file_path + "/" + 'Sets_GeneratorsOfNode.tab',format="set", set=model.GeneratorsOfNode)
	data.load(filename=tab_file_path + "/" + 'Sets_StorageOfNodes.tab',format="set", set=model.StoragesOfNode)
	data.load(filename=tab_file_path + "/" + 'Sets_PipelineType.tab',format="set", set=model.PipelineType)
	data.load(filename=tab_file_path + "/" + 'Sets_PipelineTypeOfLines.tab',format="set", set=model.PipelineTypeOfDirectionalLink)

	print("Constructing sub sets...")

	#Build arc subsets

	def NodesLinked_init(model, node):
		retval = []
		for (i,j) in model.DirectionalLink:
			if j == node:
				retval.append(i)
		return retval
	model.NodesLinked = Set(model.Node, initialize=NodesLinked_init)

	def BidirectionalArc_init(model):
		retval = []
		for (i,j) in model.DirectionalLink:
			if i != j and (not (j,i) in retval):
				retval.append((i,j))
		return retval
	model.BidirectionalArc = Set(dimen=2, initialize=BidirectionalArc_init, ordered=True) #l


	##############
	##PARAMETERS##
	##############

	#Define the parameters

	print("Declaring parameters...")

	#Scaling

	model.discountrate = Param(initialize=discountrate)
	model.WACC = Param(initialize=WACC)
	model.LeapYearsInvestment = Param(initialize=LeapYearsInvestment)
	model.operationalDiscountrate = Param(mutable=True)
	model.sceProbab = Param(model.Scenario, mutable=True)
	model.seasScale = Param(model.Season, initialize=1.0, mutable=True)
	model.lengthRegSeason = Param(initialize=lengthRegSeason)
	model.lengthPeakSeason = Param(initialize=lengthPeakSeason)

	#Cost

	model.genCapitalCost = Param(model.Generator, model.Period, default=0, mutable=True)
	model.transmissionTypeConverterFixedCost = Param(model.TransmissionType, model.Period, default=0, mutable=True)
	model.transmissionTypeConverterVariableCost = Param(model.TransmissionType, model.Period, default=0, mutable=True)
	model.transmissionTypeCableFixedCost = Param(model.TransmissionType, model.Period, default=0, mutable=True)
	model.transmissionTypeCableVariableCost = Param(model.TransmissionType, model.Period, default=0, mutable=True)
	model.storPWCapitalCost = Param(model.Storage, model.Period, default=0, mutable=True)
	model.storENCapitalCost = Param(model.Storage, model.Period, default=0, mutable=True)
	model.genFixedOMCost = Param(model.Generator, model.Period, default=0, mutable=True)
	model.transmissionTypeFixedOMCost = Param(model.TransmissionType, model.Period, default=0, mutable=True)
	model.storPWFixedOMCost = Param(model.Storage, model.Period, default=0, mutable=True)
	model.storENFixedOMCost = Param(model.Storage, model.Period, default=0, mutable=True)
	model.genInvCost = Param(model.Generator, model.Period, default=9000000, mutable=True)
	model.transmissionFixInvCost = Param(model.BidirectionalArc, model.Period, default=3000000, mutable=True)
	model.transmissionVarInvCost = Param(model.BidirectionalArc, model.Period, default=3000000, mutable=True)
	model.storPWInvCost = Param(model.Storage, model.Period, default=1000000, mutable=True)
	model.storENInvCost = Param(model.Storage, model.Period, default=800000, mutable=True)
	model.transmissionLength = Param(model.BidirectionalArc, default=0, mutable=True)
	model.genVariableOMCost = Param(model.Generator, default=0.0, mutable=True)
	model.genFuelCostRaw = Param(model.Generator, model.Period, default=0.0, mutable=True)
	model.genFuelCost = Param(model.Generator, model.Period, default=0.0, mutable=True)
	model.genMargCost = Param(model.Generator, model.Period, default=600, mutable=True)
	model.genCO2TypeFactor = Param(model.Generator, default=0.0, mutable=True)
	model.genCO2Captured = Param(model.Generator, default=0.0, mutable=True)
	model.nodeLostLoadCost = Param(model.Node, model.Period, default=22000.0)
	model.CO2price = Param(model.Period, default=0.0, mutable=True)
	# model.CCSCostTSFix = Param(initialize=1149873.72) #NB! Hard-coded
	# model.CCSCostTSVariable = Param(model.Period, default=0.0, mutable=True)
	# model.CCSRemFrac = Param(initialize=0.9)

	#Node dependent technology limitations

	model.genRefInitCap = Param(model.GeneratorsOfNode, default=0.0, mutable=True)
	model.genScaleInitCap = Param(model.Generator, model.Period, default=0.0, mutable=True)
	model.genInitCap = Param(model.GeneratorsOfNode, model.Period, default=0.0, mutable=True)
	model.transmissionInitCap = Param(model.BidirectionalArc, model.Period, default=0.0, mutable=True)
	model.storPWInitCap = Param(model.StoragesOfNode, model.Period, default=0.0, mutable=True)
	model.storENInitCap = Param(model.StoragesOfNode, model.Period, default=0.0, mutable=True)
	model.genMaxBuiltCap = Param(model.Node, model.Technology, model.Period, default=500000.0, mutable=True)
	model.transmissionMaxBuiltCap = Param(model.BidirectionalArc, model.Period, default=20000.0, mutable=True)
	model.storPWMaxBuiltCap = Param(model.StoragesOfNode, model.Period, default=500000.0, mutable=True)
	model.storENMaxBuiltCap = Param(model.StoragesOfNode, model.Period, default=500000.0, mutable=True)
	model.genMaxInstalledCapRaw = Param(model.Node, model.Technology, default=0.0, mutable=True)
	model.genMaxInstalledCap = Param(model.Node, model.Technology, model.Period, default=0.0, mutable=True)
	model.transmissionMaxInstalledCapRaw = Param(model.BidirectionalArc, model.Period, default=0.0)
	model.transmissionMaxInstalledCap = Param(model.BidirectionalArc, model.Period, default=0.0, mutable=True)
	model.storPWMaxInstalledCap = Param(model.StoragesOfNode, model.Period, default=0.0, mutable=True)
	model.storPWMaxInstalledCapRaw = Param(model.StoragesOfNode, default=0.0, mutable=True)
	model.storENMaxInstalledCap = Param(model.StoragesOfNode, model.Period, default=0.0, mutable=True)
	model.storENMaxInstalledCapRaw = Param(model.StoragesOfNode, default=0.0, mutable=True)

	#Type dependent technology limitations

	model.genLifetime = Param(model.Generator, default=0.0, mutable=True)
	model.transmissionLifetime = Param(model.BidirectionalArc, default=40.0, mutable=True)
	model.storageLifetime = Param(model.Storage, default=0.0, mutable=True)
	model.genEfficiency = Param(model.Generator, model.Period, default=1.0, mutable=True)
	model.lineEfficiency = Param(model.DirectionalLink, default=0.97, mutable=True)
	model.storageChargeEff = Param(model.Storage, default=1.0, mutable=True)
	model.storageDischargeEff = Param(model.Storage, default=1.0, mutable=True)
	model.storageBleedEff = Param(model.Storage, default=1.0, mutable=True)
	model.genRampUpCap = Param(model.ThermalGenerators, default=0.0, mutable=True)
	model.storageDiscToCharRatio = Param(model.Storage, default=1.0, mutable=True) #NB! Hard-coded
	model.storagePowToEnergy = Param(model.DependentStorage, default=1.0, mutable=True)


	#Stochastic input

	model.sloadRaw = Param(model.Node, model.Operationalhour, model.Scenario, model.Period, default=0.0, mutable=True)
	model.sloadAnnualDemand = Param(model.Node, model.Period, default=0.0, mutable=True)
	model.sload = Param(model.Node, model.Operationalhour, model.Period, model.Scenario, default=0.0, mutable=True)
	model.genCapAvailTypeRaw = Param(model.Generator, default=1.0, mutable=True)
	model.genCapAvailStochRaw = Param(model.GeneratorsOfNode, model.Operationalhour, model.Scenario, model.Period, default=0.0, mutable=True)
	model.genCapAvail = Param(model.GeneratorsOfNode, model.Operationalhour, model.Scenario, model.Period, default=0.0, mutable=True)
	model.maxRegHydroGenRaw = Param(model.Node, model.Period, model.HoursOfSeason, model.Scenario, default=1.0, mutable=True)
	model.maxRegHydroGen = Param(model.Node, model.Period, model.Season, model.Scenario, default=1.0, mutable=True)
	model.maxHydroNode = Param(model.Node, default=0.0, mutable=True)
	model.storOperationalInit = Param(model.Storage, default=0.0, mutable=True) #Percentage of installed energy capacity initially

	if EMISSION_CAP:
		model.CO2cap = Param(model.Period, default=5000.0, mutable=True)

	#SÆVAREID: Coordinates for map visualization
	model.Latitude = Param(model.Node, default=0.0, mutable=True)
	model.Longitude = Param(model.Node, default=0.0, mutable=True)

	#Load the parameters

	print("Reading parameters...")

	data.load(filename=tab_file_path + "/" + 'Generator_CapitalCosts.tab', param=model.genCapitalCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_FixedOMCosts.tab', param=model.genFixedOMCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_VariableOMCosts.tab', param=model.genVariableOMCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_FuelCosts.tab', param=model.genFuelCostRaw, format="table")
	# data.load(filename=tab_file_path + "/" + 'Generator_CCSCostTSVariable.tab', param=model.CCSCostTSVariable, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_Efficiency.tab', param=model.genEfficiency, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_RefInitialCap.tab', param=model.genRefInitCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_ScaleFactorInitialCap.tab', param=model.genScaleInitCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_InitialCapacity.tab', param=model.genInitCap, format="table") #node_generator_intial_capacity.xlsx
	data.load(filename=tab_file_path + "/" + 'Generator_MaxBuiltCapacity.tab', param=model.genMaxBuiltCap, format="table")#?
	data.load(filename=tab_file_path + "/" + 'Generator_MaxInstalledCapacity.tab', param=model.genMaxInstalledCapRaw, format="table")#maximum_capacity_constraint_040317_high
	data.load(filename=tab_file_path + "/" + 'Generator_CO2Content.tab', param=model.genCO2TypeFactor, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_CO2Captured.tab', param=model.genCO2Captured, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_RampRate.tab', param=model.genRampUpCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_GeneratorTypeAvailability.tab', param=model.genCapAvailTypeRaw, format="table")
	data.load(filename=tab_file_path + "/" + 'Generator_Lifetime.tab', param=model.genLifetime, format="table")

	data.load(filename=tab_file_path + "/" + 'Transmission_InitialCapacity.tab', param=model.transmissionInitCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_MaxBuiltCapacity.tab', param=model.transmissionMaxBuiltCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_MaxInstallCapacityRaw.tab', param=model.transmissionMaxInstalledCapRaw, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_Length.tab', param=model.transmissionLength, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_TypeConverterFixedCost.tab', param=model.transmissionTypeConverterFixedCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_TypeConverterVariableCost.tab', param=model.transmissionTypeConverterVariableCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_TypeCableFixedCost.tab', param=model.transmissionTypeCableFixedCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_TypeCableVariableCost.tab', param=model.transmissionTypeCableVariableCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_TypeFixedOMCost.tab', param=model.transmissionTypeFixedOMCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_lineEfficiency.tab', param=model.lineEfficiency, format="table")
	data.load(filename=tab_file_path + "/" + 'Transmission_Lifetime.tab', param=model.transmissionLifetime, format="table")

	data.load(filename=tab_file_path + "/" + 'Storage_StorageBleedEfficiency.tab', param=model.storageBleedEff, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_StorageChargeEff.tab', param=model.storageChargeEff, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_StorageDischargeEff.tab', param=model.storageDischargeEff, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_StoragePowToEnergy.tab', param=model.storagePowToEnergy, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_EnergyCapitalCost.tab', param=model.storENCapitalCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_EnergyFixedOMCost.tab', param=model.storENFixedOMCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_EnergyInitialCapacity.tab', param=model.storENInitCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_EnergyMaxBuiltCapacity.tab', param=model.storENMaxBuiltCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_EnergyMaxInstalledCapacity.tab', param=model.storENMaxInstalledCapRaw, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_StorageInitialEnergyLevel.tab', param=model.storOperationalInit, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_PowerCapitalCost.tab', param=model.storPWCapitalCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_PowerFixedOMCost.tab', param=model.storPWFixedOMCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_InitialPowerCapacity.tab', param=model.storPWInitCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_PowerMaxBuiltCapacity.tab', param=model.storPWMaxBuiltCap, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_PowerMaxInstalledCapacity.tab', param=model.storPWMaxInstalledCapRaw, format="table")
	data.load(filename=tab_file_path + "/" + 'Storage_Lifetime.tab', param=model.storageLifetime, format="table")

	data.load(filename=tab_file_path + "/" + 'Node_NodeLostLoadCost.tab', param=model.nodeLostLoadCost, format="table")
	data.load(filename=tab_file_path + "/" + 'Node_ElectricAnnualDemand.tab', param=model.sloadAnnualDemand, format="table")
	data.load(filename=tab_file_path + "/" + 'Node_HydroGenMaxAnnualProduction.tab', param=model.maxHydroNode, format="table")

	#SÆVAREID: Coordinates
	data.load(filename=tab_file_path + "/" + 'Node_Latitude.tab', param=model.Latitude, format="table")
	data.load(filename=tab_file_path + "/" + 'Node_Longitude.tab', param=model.Longitude, format="table")

	if scenariogeneration:
		scenariopath = tab_file_path
	else:
		scenariopath = scenario_data_path

	# data.load(filename=scenariopath + "/" + f'Stochastic_HydroGenMaxSeasonalProduction.tab', param=model.maxRegHydroGenRaw, format="table")
	# data.load(filename=scenariopath + "/" + f'Stochastic_StochasticAvailability.tab', param=model.genCapAvailStochRaw, format="table")
	# data.load(filename=scenariopath + "/" + f'Stochastic_ElectricLoadRaw.tab', param=model.sloadRaw, format="table")
	data.load(filename=scenariopath + "/" + f'Stochastic_HydroGenMaxSeasonalProduction_h2_{NoOfHydrogenScenarios}_scen_{NoOfNormalScenarios}.tab', param=model.maxRegHydroGenRaw, format="table")
	data.load(filename=scenariopath + "/" + f'Stochastic_StochasticAvailability_h2_{NoOfHydrogenScenarios}_scen_{NoOfNormalScenarios}.tab', param=model.genCapAvailStochRaw, format="table")
	data.load(filename=scenariopath + "/" + f'Stochastic_ElectricLoadRaw_h2_{NoOfHydrogenScenarios}_scen_{NoOfNormalScenarios}.tab', param=model.sloadRaw, format="table")

	# data.load(filename=tab_file_path + "/" + 'General_seasonScale.tab', param=model.seasScale, format="table")

	if EMISSION_CAP:
		data.load(filename=tab_file_path + "/" + 'General_CO2Cap.tab', param=model.CO2cap, format="table")
	# else:
	# 	data.load(filename=tab_file_path + "/" + 'General_CO2Price.tab', param=model.CO2price, format="table")
	data.load(filename=tab_file_path + "/" + 'General_CO2Price.tab', param=model.CO2price, format="table")

	print("Constructing parameter values...")

	def prepSceProbab_rule(model):
		#Build an equiprobable probability distribution for scenarios

		for sce in model.Scenario:
			model.sceProbab[sce] = value(1/len(model.Scenario))

	model.build_SceProbab = BuildAction(rule=prepSceProbab_rule)

	def prepSeasScale(model):
		for s in model.Season:
			if s in ["winter", "spring", "summer", "fall"]:
				model.seasScale[s] = (8760 - lengthPeakSeason * NoOfPeakSeason) / (NoOfRegSeason * lengthRegSeason)
			else:
				model.seasScale[s] = 1
	model.build_seasScale = BuildAction(rule=prepSeasScale)

	def prepInvCost_rule(model):
		#Build investment cost for generators, storages and transmission. Annual cost is calculated for the lifetime of the generator and discounted for a year.
		#Then cost is discounted for the investment period (or the remaining lifetime). CCS generators has additional fixed costs depending on emissions.

		#Generator
		for g in model.Generator:
			for i in model.Period:
				costperyear=(model.WACC / (1 - ((1+model.WACC) ** (1-model.genLifetime[g])))) * model.genCapitalCost[g,i] + model.genFixedOMCost[g,i]
				costperperiod = costperyear * 1000 * (1 - (1+model.discountrate) **-(min(value((len(model.Period)-i+1)*5), value(model.genLifetime[g]) )))/ (1 - (1 / (1 + model.discountrate)))
				# if ('CCS',g) in model.GeneratorsOfTechnology:
					# 	costperperiod+=model.CCSCostTSFix*model.CCSRemFrac*model.genCO2TypeFactor[g]*(3.6/model.genEfficiency[g,i])
				model.genInvCost[g,i]=costperperiod

		#Storage
		for b in model.Storage:
			for i in model.Period:
				costperyearPW=(model.WACC/(1-((1+model.WACC)**(1-model.storageLifetime[b]))))*model.storPWCapitalCost[b,i]+model.storPWFixedOMCost[b,i]
				costperperiodPW=costperyearPW*1000*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.storageLifetime[b]))))/(1-(1/(1+model.discountrate)))
				model.storPWInvCost[b,i]=costperperiodPW
				costperyearEN=(model.WACC/(1-((1+model.WACC)**(1-model.storageLifetime[b]))))*model.storENCapitalCost[b,i]+model.storENFixedOMCost[b,i]
				costperperiodEN=costperyearEN*1000*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.storageLifetime[b]))))/(1-(1/(1+model.discountrate)))
				model.storENInvCost[b,i]=costperperiodEN

		#Transmission
		for (n1,n2) in model.BidirectionalArc:
			for i in model.Period:
				for t in model.TransmissionType:
					if (n1,n2,t) in model.TransmissionTypeOfDirectionalLink:
						costperyear=(model.WACC/(1-((1+model.WACC)**(1-model.transmissionLifetime[n1,n2]))))*(model.transmissionLength[n1,n2]*model.transmissionTypeCableVariableCost[t,i] + model.transmissionTypeConverterVariableCost[t,i]) + model.transmissionLength[n1,n2]* model.transmissionTypeFixedOMCost[t,i]
						costperperiod=costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.transmissionLifetime[n1,n2]))))/(1-(1/(1+model.discountrate)))
						model.transmissionVarInvCost[n1,n2,i]=costperperiod

						costperyear=(model.WACC/(1-((1+model.WACC)**(1-model.transmissionLifetime[n1,n2]))))*(model.transmissionLength[n1,n2]*model.transmissionTypeCableFixedCost[t,i] + model.transmissionTypeConverterFixedCost[t,i])
						costperperiod=costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.transmissionLifetime[n1,n2]))))/(1-(1/(1+model.discountrate)))
						model.transmissionFixInvCost[n1,n2,i]=costperperiod
	model.build_InvCost = BuildAction(rule=prepInvCost_rule)

	def prepFuelCost_rule(model):
		for g in model.Generator:
			for i in model.Period:
				if "gas" in g.lower():
					model.genFuelCost[g,i] = model.genFuelCostRaw[g,i]
				else:
					model.genFuelCost[g,i] = model.genFuelCostRaw[g,i]
	model.build_genFuelCost = BuildAction(rule=prepFuelCost_rule)

	def prepOperationalCostGen_rule(model):
		#Build generator short term marginal costs

		for g in model.Generator:
			for i in model.Period:
				if not EMISSION_CAP:
					if ('CCS',g) in model.GeneratorsOfTechnology:
						#CCSRemFrac seems to be leftover code, and I think it should be removed. The 90% reduction is already accounted for in the input data for genCO2TypeFactor

						# costperenergyunit=(3.6/model.genEfficiency[g,i])*(model.genFuelCost[g,i]+(1-model.CCSRemFrac)*model.genCO2TypeFactor[g]*model.CO2price[i])+ \
						# 				  (3.6/model.genEfficiency[g,i])*(model.CCSRemFrac*model.genCO2TypeFactor[g]*model.CCSCostTSVariable[i])+ \
						# 				  model.genVariableOMCost[g]
						costperenergyunit=(3.6/model.genEfficiency[g,i])*(model.genFuelCost[g,i]+model.genCO2TypeFactor[g]*model.CO2price[i])+ \
										  model.genVariableOMCost[g]
					else:
						costperenergyunit=(3.6/model.genEfficiency[g,i])*(model.genFuelCost[g,i]+model.genCO2TypeFactor[g]*model.CO2price[i])+ \
										  model.genVariableOMCost[g]
				else:
					if ('CCS',g) in model.GeneratorsOfTechnology:
						# costperenergyunit=(3.6/model.genEfficiency[g,i])*(model.genFuelCost[g,i])+ \
						# 				  (3.6/model.genEfficiency[g,i])*(model.CCSRemFrac*model.genCO2TypeFactor[g]*model.CCSCostTSVariable[i])+ \
						# 				  model.genVariableOMCost[g]
						costperenergyunit=(3.6/model.genEfficiency[g,i])*(model.genFuelCost[g,i])+ \
										  model.genVariableOMCost[g]
					else:
						costperenergyunit=(3.6/model.genEfficiency[g,i])*(model.genFuelCost[g,i])+ \
										  model.genVariableOMCost[g]
				model.genMargCost[g,i]=costperenergyunit

	model.build_OperationalCostGen = BuildAction(rule=prepOperationalCostGen_rule)

	def prepInitialCapacityNodeGen_rule(model):
		#Build initial capacity for generator type in node

		for (n,g) in model.GeneratorsOfNode:
			for i in model.Period:
				if value(model.genInitCap[n,g,i]) == 0:
					model.genInitCap[n,g,i] = model.genRefInitCap[n,g]*(1-model.genScaleInitCap[g,i])

	model.build_InitialCapacityNodeGen = BuildAction(rule=prepInitialCapacityNodeGen_rule)

	def prepInitialCapacityTransmission_rule(model):
		#Build initial capacity for transmission lines to ensure initial capacity is the upper installation bound if infeasible

		for (n1,n2) in model.BidirectionalArc:
			for i in model.Period:
				if value(model.transmissionMaxInstalledCapRaw[n1,n2,i]) <= value(model.transmissionInitCap[n1,n2,i]):
					model.transmissionMaxInstalledCap[n1,n2,i] = model.transmissionInitCap[n1,n2,i]
				else:
					model.transmissionMaxInstalledCap[n1,n2,i] = model.transmissionMaxInstalledCapRaw[n1,n2,i]
	model.build_InitialCapacityTransmission = BuildAction(rule=prepInitialCapacityTransmission_rule)

	def prepOperationalDiscountrate_rule(model):
		#Build operational discount rate

		model.operationalDiscountrate = sum((1+model.discountrate)**(-j) for j in list(range(0,value(model.LeapYearsInvestment))))
	model.build_operationalDiscountrate = BuildAction(rule=prepOperationalDiscountrate_rule)

	def prepGenMaxInstalledCap_rule(model):
		#Build resource limit (installed limit) for all periods. Avoid infeasibility if installed limit lower than initially installed cap.

		for t in model.Technology:
			for n in model.Node:
				for i in model.Period:
					if value(model.genMaxInstalledCapRaw[n,t] <= sum(model.genInitCap[n,g,i] for g in model.Generator if (n,g) in model.GeneratorsOfNode and (t,g) in model.GeneratorsOfTechnology)):
						model.genMaxInstalledCap[n,t,i]=sum(model.genInitCap[n,g,i] for g in model.Generator if (n,g) in model.GeneratorsOfNode and (t,g) in model.GeneratorsOfTechnology)
					else:
						model.genMaxInstalledCap[n,t,i]=model.genMaxInstalledCapRaw[n,t]
	model.build_genMaxInstalledCap = BuildAction(rule=prepGenMaxInstalledCap_rule)

	def storENMaxInstalledCap_rule(model):
		#Build installed limit (resource limit) for storEN

		#Why is this here? Why not just use storENMaxInstalledCapRaw in the constraints?

		for (n,b) in model.StoragesOfNode:
			for i in model.Period:
				model.storENMaxInstalledCap[n,b,i]=model.storENMaxInstalledCapRaw[n,b]

	model.build_storENMaxInstalledCap = BuildAction(rule=storENMaxInstalledCap_rule)

	def storPWMaxInstalledCap_rule(model):
		#Build installed limit (resource limit) for storPW

		#Why is this here? Why not just use storPWMaxInstalledCapRaw in the constraints?

		for (n,b) in model.StoragesOfNode:
			for i in model.Period:
				model.storPWMaxInstalledCap[n,b,i]=model.storPWMaxInstalledCapRaw[n,b]

	model.build_storPWMaxInstalledCap = BuildAction(rule=storPWMaxInstalledCap_rule)

	def prepRegHydro_rule(model):
		#Build hydrolimits for all periods

		for n in model.Node:
			for s in model.Season:
				for i in model.Period:
					for sce in model.Scenario:
						model.maxRegHydroGen[n,i,s,sce]=sum(model.maxRegHydroGenRaw[n,i,s,h,sce] for h in model.Operationalhour if (s,h) in model.HoursOfSeason)

	model.build_maxRegHydroGen = BuildAction(rule=prepRegHydro_rule)

	def prepGenCapAvail_rule(model):
		#Build generator availability for all periods

		for (n,g) in model.GeneratorsOfNode:
			for h in model.Operationalhour:
				for s in model.Scenario:
					for i in model.Period:
						if value(model.genCapAvailTypeRaw[g]) == 0:
							if value(model.genCapAvailStochRaw[n,g,h,s,i]) >= 0.001:
								model.genCapAvail[n,g,h,s,i] = model.genCapAvailStochRaw[n,g,h,s,i]
							else:
								model.genCapAvail[n,g,h,s,i] = 0
						else:
							model.genCapAvail[n,g,h,s,i]=model.genCapAvailTypeRaw[g]

	model.build_genCapAvail = BuildAction(rule=prepGenCapAvail_rule)

	def prepSload_rule(model):
		#Build load profiles for all periods

		counter = 0
		f = open(result_file_path + '/AdjustedNegativeLoad_' + name + '.txt', 'w')
		for n in model.Node:
			for i in model.Period:
				nodeaverageload = 0
				for h in model.Operationalhour:
					if value(h) < value(model.FirstHoursOfRegSeason[-1] + model.lengthRegSeason):
						for sce in model.Scenario:
							nodeaverageload += model.sloadRaw[n, h, sce, i].value
				nodeaverageload = nodeaverageload / value(
					(model.FirstHoursOfRegSeason[-1] + model.lengthRegSeason - 1) * len(model.Scenario))

				hourlyadjustment = value(model.sloadAnnualDemand[n, i].value / 8760) - value(nodeaverageload)
				for h in model.Operationalhour:
					for sce in model.Scenario:
						if value(model.sloadRaw[n, h, sce, i].value + hourlyadjustment) > 0:
							model.sload[n, h, i, sce] = model.sloadRaw[n, h, sce, i].value + hourlyadjustment
							if value(model.sload[n,h,i,sce]) < 0:
								f.write('Adjusted electricity load: ' + str(value(model.sload[n,h,i,sce])) + ', 10 MW for hour ' + str(h) + ' and scenario ' + str(sce) + ' in ' + str(n) + "\n")
								model.sload[n,h,i,sce] = 10
								counter += 1
						else:
							f.write('Adjusted electricity load: ' + str(value(model.sloadRaw[n,h,sce,i].value + hourlyadjustment)) + ', 0 MW for hour ' + str(h) + ' and scenario ' + str(sce) + ' in ' + str(n) + "\n")
							model.sload[n,h,i,sce] = 0
							counter += 1
		f.write('Hours with too small raw electricity load: ' + str(counter))
		f.close()

	model.build_sload = BuildAction(rule=prepSload_rule)

	stopReading = startConstraints = datetime.now()
	print("Sets and parameters declared and read...")

	#############
	##VARIABLES##
	#############

	print("Declaring variables...")

	model.genInvCap = Var(model.GeneratorsOfNode, model.Period, domain=NonNegativeReals)
	model.transmissionFixedCostPaid = Var(model.BidirectionalArc, model.Period, domain=Binary)
	model.transmissionInvCap = Var(model.BidirectionalArc, model.Period, domain=NonNegativeReals)
	model.storPWInvCap = Var(model.StoragesOfNode, model.Period, domain=NonNegativeReals)
	model.storENInvCap = Var(model.StoragesOfNode, model.Period, domain=NonNegativeReals)
	model.genOperational = Var(model.GeneratorsOfNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
	model.storOperational = Var(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
	model.transmissionOperational = Var(model.DirectionalLink, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals) #flow
	model.storCharge = Var(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
	model.storDischarge = Var(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
	model.loadShed = Var(model.Node, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
	model.genInstalledCap = Var(model.GeneratorsOfNode, model.Period, domain=NonNegativeReals)
	model.transmissionInstalledCap = Var(model.BidirectionalArc, model.Period, domain=NonNegativeReals)
	model.storPWInstalledCap = Var(model.StoragesOfNode, model.Period, domain=NonNegativeReals)
	model.storENInstalledCap = Var(model.StoragesOfNode, model.Period, domain=NonNegativeReals)

	if hydrogen is True:
		#Hydrogen sets
		model.HydrogenProdNode = Set(ordered=True, within=model.Node)
		model.ReformerLocations = Set(ordered=True, within=model.HydrogenProdNode)
		model.ReformerPlants = Set(ordered=True)


		#Reading sets
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ProductionNodes.tab', format="set", set=model.HydrogenProdNode)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerLocations.tab', format="set", set=model.ReformerLocations)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerPlants.tab', format="set", set=model.ReformerPlants)


		def HydrogenLinks_init(model):
			retval= []
			for (n1,n2) in model.DirectionalLink:
				if n1 in model.HydrogenProdNode and n2 in model.HydrogenProdNode:
					retval.append((n1,n2))
			return retval
		model.AllowedHydrogenLinks = Set(dimen=2, initialize=HydrogenLinks_init, ordered=True)
		# model.AllowedHydrogenLinks = Set(dimen=2, within=model.Node * model.Node, ordered=True) # Depcreated; The links are now instead defined by the transmission links, but only between the production nodes
		# data.load(filename=tab_file_path + '/' + 'Hydrogen_Links.tab', format="set", set=model.AllowedHydrogenLinks) # Deprecated; The links are now instead defined by the transmission links, but only between the production nodes

		def HydrogenBidirectionPipelines_init(model):
			retval = []
			for (n1,n2) in model.BidirectionalArc:
				if n1 in model.HydrogenProdNode and n2 in model.HydrogenProdNode:
					retval.append((n1,n2))
			return retval
		model.HydrogenBidirectionPipelines = Set(dimen=2, initialize=HydrogenBidirectionPipelines_init, ordered=True)

		def HydrogenLinks_init(model, node):
			retval = []
			for (i,j) in model.AllowedHydrogenLinks:
				if j == node:
					retval.append(i)
			return retval
		model.HydrogenLinks = Set(model.Node, initialize=HydrogenLinks_init)

		# Hydrogen parameters
		model.hydrogenDemandRaw = Param(model.HydrogenProdNode, model.Period, default=0, mutable=True)
		model.hydrogenDemand = Param (model.HydrogenProdNode, model.Period, model.Scenario, default=0, mutable=True)

		model.elyzerPlantCapitalCostRaw = Param(model.Period, default=99999, mutable=True)
		model.elyzerPlantCapitalCost = Param(model.Period, default=99999, mutable=True)
		model.elyzerStackCapitalCostRaw = Param(model.Period, default=99999, mutable=True)
		model.elyzerStackCapitalCost = Param(model.Period, default=99999, mutable=True)
		model.elyzerFixedOMCostRaw = Param(model.Period, default=99999, mutable=True)
		model.elyzerFixedOMCost = Param(model.Period, default=99999, mutable=True)
		model.elyzerPowerConsumptionPerTon = Param(model.Period, default=99999, mutable=True)
		model.elyzerLifetime = Param(default=20, mutable=True)
		model.elyzerInvCost = Param(model.Period, default=99999, mutable=True)

		model.ReformerPlantsCapitalCostRaw = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantsCapitalCost = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantFixedOMCostRaw = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantFixedOMCost = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantVarOMCostRaw = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantVarOMCost = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantInvCost = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantEfficiency = Param(model.ReformerPlants, model.Period, default=0, mutable=True)
		model.ReformerPlantElectricityUse = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerPlantLifetime = Param(model.ReformerPlants, default=25, mutable=True)
		model.ReformerEmissionFactor = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerCO2CaptureFactor = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)
		model.ReformerMargCost = Param(model.ReformerPlants, model.Period, default=99999, mutable=True)

		model.hydrogenPipelineLifetime = Param(default=40)
		model.hydrogenPipelineCapCost = Param(model.PipelineType, default=99999, mutable=True)
		model.hydrogenPipelineOMCost = Param(model.PipelineType, default=99999, mutable=True)
		model.hydrogenPipelineCapacity = Param(model.PipelineType, default=0, mutable=True)
		model.hydrogenPipelineInvCost = Param(model.HydrogenBidirectionPipelines, model.PipelineType, default=999999, mutable=True)
		model.PipelineLength = Param(model.HydrogenBidirectionPipelines, mutable=True, default=9999)
		model.hydrogenPipelineCompressorElectricityUsage = Param(model.PipelineType, default=99999, mutable=True)
		model.hydrogenPipelinePowerDemandPerTon = Param(model.HydrogenBidirectionPipelines, model.PipelineType, default=99999, mutable=True)


		if h2storage is False:
			#Cost of storing the produced hydrogen intraseasonally. Have to have this because we have implicit free storage without.
			model.averageHydrogenSeasonalStorageCost = Param(default=0.35, mutable=True) #Source: Levelized cost of storage from Table 5 in Picturing the value of underground gas storage to the European hydrogen system by Gas Infrastructure Europe (GIE)
		else:
			model.hydrogenMaxStorageCapacity = Param(model.HydrogenProdNode, default=0, mutable=True)
			model.hydrogenStorageCapitalCost = Param(model.Period, default=99999, mutable=True)
			model.hydrogenStorageFixedOMCost = Param(model.Period, default=99999, mutable=True)
			model.hydrogenStorageInvCost = Param(model.Period, default=99999, mutable=True)
			model.hydrogenStorageInitOperational = Param(default=0.5)
			model.hydrogenStorageLifetime = Param(default=30)

		model.hydrogenLHV_kg = Param(default=33.3/1000, mutable=False) #LHV of hydrogen is 33.3 kWh/kg = 0.0333 MWh / kg = 33.3 MWh/ton

		#Hydrogen variables
		#Operational
		model.meetHydrogenDemand = Var(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
		model.hydrogenProducedElectro = Var(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
		model.hydrogenProducedReformer_ton = Var(model.ReformerLocations, model.ReformerPlants, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
		model.hydrogenProducedReformer_MWh = Var(model.ReformerLocations, model.ReformerPlants, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
		model.hydrogenSentPipeline = Var(model.PipelineTypeOfDirectionalLink, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
		model.powerForHydrogen = Var(model.HydrogenProdNode, model.Period, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals) #Two period indexes because one describes the year it was bought (the first index), the other describes when it is used (second index)
		model.h2_demand_shed = Var(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)

		if h2storage is True:
			model.hydrogenStorageOperational = Var(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
			model.hydrogenChargeStorage = Var(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
			model.hydrogenDischargeStorage = Var(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)

		model.hydrogenForPower = Var(model.HydrogenGenerators ,model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals,initialize=0.0)

		#Strategic
		model.elyzerCapBuilt = Var(model.HydrogenProdNode, model.Period, model.Period, domain=NonNegativeReals)
		model.elyzerTotalCap = Var(model.HydrogenProdNode, model.Period, model.Period, domain=NonNegativeReals)
		model.ReformerCapBuilt = Var(model.ReformerLocations, model.ReformerPlants, model.Period, domain=NonNegativeReals) #Capacity  of MW H2 production built in period i
		model.ReformerTotalCap = Var(model.ReformerLocations, model.ReformerPlants, model.Period, domain=NonNegativeReals) #Total capacity of MW H2 production
		# model.hydrogenPipelineBuilt = Var(model.HydrogenBidirectionPipelines, model.PipelineType, model.Period, domain=NonNegativeReals)
		model.hydrogenPipelineBuilt = Var(model.HydrogenBidirectionPipelines, model.PipelineType, model.Period, domain=NonNegativeIntegers)
		model.totalHydrogenPipelineCapacity = Var(model.HydrogenBidirectionPipelines, model.PipelineType, model.Period, domain=NonNegativeReals)
		if h2storage is True:
			model.hydrogenStorageBuilt = Var(model.HydrogenProdNode, model.Period, domain=NonNegativeReals)
			model.hydrogenTotalStorage = Var(model.HydrogenProdNode, model.Period, domain=NonNegativeReals)

		#Reading parameters
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ElectrolyzerPlantCapitalCost.tab', format="table", param=model.elyzerPlantCapitalCostRaw)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ElectrolyzerStackCapitalCost.tab', format="table", param=model.elyzerStackCapitalCostRaw)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ElectrolyzerFixedOMCost.tab', format="table", param=model.elyzerFixedOMCostRaw)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ElectrolyzerMWhPerTon.tab', format="table", param=model.elyzerPowerConsumptionPerTon)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ElectrolyzerLifetime.tab', format="table", param=model.elyzerLifetime)

		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerCapitalCost.tab', format='table', param=model.ReformerPlantsCapitalCostRaw)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerFixedOMCost.tab', format='table', param=model.ReformerPlantFixedOMCostRaw)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerVariableOMCost.tab', format='table', param=model.ReformerPlantVarOMCostRaw)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerEfficiency.tab', format='table', param=model.ReformerPlantEfficiency)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerElectricityUse.tab', format='table', param=model.ReformerPlantElectricityUse)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerLifetime.tab', format='table', param=model.ReformerPlantLifetime)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerEmissionFactor.tab', format='table', param=model.ReformerEmissionFactor)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_ReformerCO2CaptureFactor.tab', format='table', param=model.ReformerCO2CaptureFactor)

		data.load(filename=tab_file_path + '/' + 'Hydrogen_PipelineCapitalCost.tab', format="table", param=model.hydrogenPipelineCapCost)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_PipelineOMCostPerKM.tab', format="table", param=model.hydrogenPipelineOMCost)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_PipelineCapacity.tab', format="table", param=model.hydrogenPipelineCapacity)
		data.load(filename=tab_file_path + '/' + 'Hydrogen_PipelineCompressorPowerUsage.tab', format="table", param=model.hydrogenPipelineCompressorElectricityUsage)
		# data.load(filename=tab_file_path + '/' + 'Hydrogen_Distances.tab', format="table", param=model.PipelineLength) # Depecreated; Distances are now copied from the transmission distances


		data.load(filename=tab_file_path + '/' + 'Hydrogen_Demand.tab', format="table", param=model.hydrogenDemandRaw) # Depecreated; Distances are now copied from the transmission distances

		if h2storage is True:
			data.load(filename=tab_file_path + '/' + 'Hydrogen_StorageCapitalCost.tab', format="table", param=model.hydrogenStorageCapitalCost)
			data.load(filename=tab_file_path + '/' + 'Hydrogen_StorageFixedOMCost.tab', format="table", param=model.hydrogenStorageFixedOMCost)
			data.load(filename=tab_file_path + '/' + 'Hydrogen_StorageMaxCapacity.tab', format="table", param=model.hydrogenMaxStorageCapacity)

		def prepPipelineLength_rule(model):
			for (n1,n2) in model.HydrogenBidirectionPipelines:
				if (n1,n2) in model.BidirectionalArc:
					model.PipelineLength[n1,n2] = model.transmissionLength[n1,n2]
				elif (n2,n1) in model.BidirectionalArc:
					model.PipelineLength[n1,n2] = model.transmissionLength[n2,n1]
				else:
					print('Error constructing hydrogen pipeline length for bidirectional pipeline ' + n1 + ' and ' + n2)
					exit()
		model.build_PipelineLength = BuildAction(rule=prepPipelineLength_rule)

		def prepElectrolyzerCost(model):
			for i in model.Period:
				model.elyzerPlantCapitalCost[i] = model.elyzerPlantCapitalCostRaw[i]
				model.elyzerStackCapitalCost[i] = model.elyzerStackCapitalCostRaw[i]
				model.elyzerFixedOMCost[i] = model.elyzerFixedOMCostRaw[i]
		model.build_electrolyzerCost = BuildAction(rule=prepElectrolyzerCost)

		def prepRefomerCost(model):
			for p in model.ReformerPlants:
				for i in model.Period:
					model.ReformerPlantsCapitalCost[p,i] = model.ReformerPlantsCapitalCostRaw[p,i]
					model.ReformerPlantFixedOMCost[p,i] = model.ReformerPlantFixedOMCostRaw[p,i]
					model.ReformerPlantVarOMCost[p,i] = model.ReformerPlantVarOMCostRaw[p,i]
		model.build_reformerCost = BuildAction(rule=prepRefomerCost)


		def prepElectrolyzerInvCost_rule(model):
			for i in model.Period:
				costperyear = (model.WACC / (1 - ((1 + model.WACC) ** (1 - model.elyzerLifetime)))) * model.elyzerPlantCapitalCost[i] + model.elyzerFixedOMCost[i] + ((model.WACC/(1-((1+model.WACC)**(1-8)))) + (model.WACC/(1-((1+model.WACC)**(1-16))))) * model.elyzerStackCapitalCost[i]
				costperperiod = costperyear * (1 - (1 + model.discountrate) ** -(min(value((len(model.Period) - i + 1) * 5), value(model.elyzerLifetime)))) / (1 - (1 / (1 + model.discountrate)))
				model.elyzerInvCost[i] = costperperiod
		model.build_elyzerInvCost = BuildAction(rule=prepElectrolyzerInvCost_rule)

		def prepReformerPlantInvCost_rule(model):
			for p in model.ReformerPlants:
				for i in model.Period:
					costperyear = (model.WACC/(1-((1+model.WACC)**(1-model.ReformerPlantLifetime[p]))))*model.ReformerPlantsCapitalCost[p,i]+model.ReformerPlantFixedOMCost[p,i]
					costperperiod = costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.ReformerPlantLifetime[p]))))/(1-(1/(1+model.discountrate)))
					model.ReformerPlantInvCost[p,i] = costperperiod
		model.build_ReformerPlantInvCost = BuildAction(rule=prepReformerPlantInvCost_rule)

		def prepReformerMargCost_rule(model):
			for p in model.ReformerPlants:
				for i in model.Period:
					natural_gas_price = model.genFuelCost['GasOCGT',i]*3600/1000 # Need to convert from GJ to MWh.
					NG_cost_per_kg_H2 = (natural_gas_price/model.ReformerPlantEfficiency[p,i])*model.hydrogenLHV_kg
					model.ReformerMargCost[p,i] = NG_cost_per_kg_H2 * 1000 + model.ReformerPlantVarOMCost[p,i]
		model.build_ReformerMargCost = BuildAction(rule=prepReformerMargCost_rule)

		def reformer_operational_cost_rule(model, i):
			return sum(model.operationalDiscountrate*model.seasScale[s]*model.sceProbab[w]*model.ReformerMargCost[p,i]*model.hydrogenProducedReformer_ton[n,p,h,i,w] for p in model.ReformerPlants for n in model.ReformerLocations for (s,h) in model.HoursOfSeason for w in model.Scenario)
		model.reformerOperationalCost = Expression(model.Period, rule=reformer_operational_cost_rule)

		def reformer_emissions_rule(model,i,w): #Calculates kg of CO2 emissions per kg H2 produced with Reformer
			return sum(model.seasScale[s]*model.hydrogenProducedReformer_ton[n,p,h,i,w]*model.ReformerEmissionFactor[p,i] for (s,h) in model.HoursOfSeason for n in model.ReformerLocations for p in model.ReformerPlants)
		model.reformerEmissions = Expression(model.Period, model.Scenario, rule=reformer_emissions_rule)

		def reformer_emissions_cost_rule(model,i):
			return sum(model.operationalDiscountrate*model.sceProbab[w]*model.reformerEmissions[i,w]*model.CO2price[i] for w in model.Scenario)
		model.reformerEmissionsCost = Expression(model.Period, rule=reformer_emissions_cost_rule)

		def CO2_captured_generators_rule(model, n, h, i, w):
			return sum(model.genCO2Captured[g] * model.genOperational[n,g,h,i,w] * 3.6 / model.genEfficiency[g,i] for g in model.Generator if (n,g) in model.GeneratorsOfNode)
		model.co2_captured_generators = Expression(model.Node, model.Operationalhour, model.Period, model.Scenario, rule=CO2_captured_generators_rule)

		def CO2_captured_reformers_rule(model, n, h, i, w):
			return sum(model.ReformerCO2CaptureFactor[r,i] * model.hydrogenProducedReformer_ton[n,r,h,i,w] for r in model.ReformerPlants)
		model.co2_captured_reformers = Expression(model.ReformerLocations, model.Operationalhour, model.Period, model.Scenario, rule=CO2_captured_reformers_rule)

		def prepPipelineInvcost_rule(model):
			for i in model.Period:
				for (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
					if (n1,n2) in model.HydrogenBidirectionPipelines:
						# for t in model.PipelineType:
						# 	if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
						costperyear= (model.WACC/(1-((1+model.WACC)**(1-model.hydrogenPipelineLifetime))))*model.PipelineLength[n1,n2]*(model.hydrogenPipelineCapCost[t]) + model.PipelineLength[n1,n2]*model.hydrogenPipelineOMCost[t]
						costperperiod =costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.hydrogenPipelineLifetime))))/(1-(1/(1+model.discountrate)))
						model.hydrogenPipelineInvCost[n1,n2,t] = costperperiod
					# elif (n2,n1) in model.BidirectionalArc:
					# 	costperyear= (model.WACC/(1-((1+model.WACC)**(1-model.hydrogenPipelineLifetime))))*model.PipelineLength[n2,n1]*(model.hydrogenPipelineCapCost[t]) + model.PipelineLength[n2,n1]*model.hydrogenPipelineOMCost[t]
					# 	costperperiod =costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.hydrogenPipelineLifetime))))/(1-(1/(1+model.discountrate)))
					# 	model.hydrogenPipelineInvCost[n2,n1,t] = costperperiod
		model.build_pipelineInvCost = BuildAction(rule=prepPipelineInvcost_rule)

		if h2storage is True:
			def prepHydrogenStorageInvcost_rule(model):
				for i in model.Period:
					costperyear =(model.WACC/(1-((1+model.WACC)**(1-model.hydrogenStorageLifetime))))*model.hydrogenStorageCapitalCost[i]+model.hydrogenStorageFixedOMCost[i]
					costperperiod = costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.hydrogenStorageLifetime))))/(1-(1/(1+model.discountrate)))
					model.hydrogenStorageInvCost[i] = costperperiod
			model.build_hydrogenStorageInvCost = BuildAction(rule=prepHydrogenStorageInvcost_rule)

		def prepHydrogenDemand_rule(model):
			for n in model.HydrogenProdNode:
				for i in model.Period:
					mean = value(model.hydrogenDemandRaw[n,i])
					std_dev = std_dev_percentage * mean
					for w_count, w in enumerate(model.Scenario, start=1):
						if len(Scenario) > NoOfNormalScenarios:
							# model.hydrogenDemand[n,i,w] = NoOfNormalScenarios/len(Scenario) * ceil(w_count/NoOfHydrogenScenarios) * model.hydrogenDemandRaw[n,i]
							# model.hydrogenDemand[n,i,w] = NoOfNormalScenarios/len(Scenario) * ceil(w_count/NoOfNormalScenarios) * model.hydrogenDemandRaw[n,i]
							if mean > 0 and std_dev > 0:
								h2_scen = ceil(w_count/NoOfNormalScenarios)
								demand = norm.ppf((h2_scen / NoOfHydrogenScenarios -  (h2_scen - (h2_scen - 1)) / (2 * NoOfHydrogenScenarios)), mean, std_dev)
							else:
								demand = 0
							if demand > 0:
								model.hydrogenDemand[n,i,w] = hydrogen_demand_percentage * demand
							else:
								model.hydrogenDemand[n,i,w] = 0
						else:
							model.hydrogenDemand[n,i,w] = hydrogen_demand_percentage * model.hydrogenDemandRaw[n,i]
		model.build_hydrogenDemand = BuildAction(rule=prepHydrogenDemand_rule)

		def prepHydrogenCompressorElectricityUsage_rule(model):
			for (n1,n2) in model.HydrogenBidirectionPipelines:
				for t in model.PipelineType:
					if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
						model.hydrogenPipelinePowerDemandPerTon[n1,n2,t] = model.PipelineLength[n1,n2] * model.hydrogenPipelineCompressorElectricityUsage[t]
		model.build_hydrogenPipelineCompressorPowerDemand = BuildAction(rule=prepHydrogenCompressorElectricityUsage_rule)

		# CO2 part
		model.OnshoreNodes = Set(initialize=model.HydrogenProdNode)
		model.CO2DirectionalLinks = Set(dimen=2, initialize=model.AllowedHydrogenLinks)
		model.CO2BidirectionalPipelines = Set(dimen=2, initialize=model.HydrogenBidirectionPipelines)
		model.CO2SequestrationNodes = Set(within=model.Node, initialize = ["Norway"])

		model.CO2StorageSiteCapitalCost = Param(model.CO2SequestrationNodes, default=999999999, mutable=True)
		model.CO2StorageSiteInvCost = Param(model.CO2SequestrationNodes, model.Period, default=999999999, mutable=True)
		model.StorageSiteFixedOMCost = Param(model.CO2SequestrationNodes, default=999999999, mutable=True)
		model.CO2PipelineLifetime = Param(default=40, mutable=False)
		model.CO2PipelineCapCost = Param(model.PipelineType, default=99999, mutable=True)
		model.CO2PipelineOMCost = Param(model.PipelineType, default=99999, mutable=True)
		model.CO2PipelineInvCost = Param(model.CO2BidirectionalPipelines, model.PipelineType, model.Period, default=99999, mutable=True)
		model.CO2PipelineCapacity = Param(model.PipelineType, default=99999, mutable=True)
		model.CO2PipelineElectricityUsage = Param(model.PipelineType, default=99999, mutable=True)
		model.CO2PipelinePowerDemandPerTon = Param(model.CO2BidirectionalPipelines, model.PipelineType, default=99999, mutable=True)

		# model.CO2LiquefierCapitalCost = Param(default=99999, mutable=True)
		# model.CO2LiquefierFixedOMCost = Param(default=99999, mutable=True)
		# model.CO2LiquefierInvCost = Param(model.Period, default=99999, mutable=True)
		# model.CO2LiquefierLifetime = Param(default=20, mutable=True)
		# model.CO2LiquefierElectricityUse = Param(default=1, mutable=True)

		# model.CO2LiquidStorageCapitalCost = Param(default=99999, mutable=True)
		# model.CO2LiquidStorageFixedOMCost = Param(default=99999, mutable=True)
		# model.CO2LiquidStorageInvCost = Param(model.Period, default=99999, mutable=True)
		# model.CO2LiquidStorageLifetime = Param(default=20, mutable=True)
		# model.CO2LiquidStorageInitOperational = Param(default=0.5, mutable=True)

		# model.CO2LiquidShipCapitalCost = Param(default=99999, mutable=True)
		# model.CO2LiquidShipFixedOMCost = Param(default=99999, mutable=True)
		# model.CO2LiquidShipInvCost = Param(model.Period, default=99999, mutable=True)
		# model.CO2LiquidShipVarCost = Param(default=99999, mutable=True)
		# model.CO2LiquidShipVarUnitCost = Param(model.CO2DirectionalLinks,default=99999, mutable=True)
		# model.CO2LiquidShipCapacity = Param(default=99999, mutable=True)
		# model.CO2LiquidShipLifetime = Param(default=20, mutable=True)
		# model.CO2LiquidShipLoadingDischargeTime = Param(default=100, mutable=True)
		# model.CO2LiquidShipSpeed = Param(default=27.78, mutable=False) # 15 nm/hr ~ 27.78 km/h shipping speed
		# model.ShippingTimes = Param(model.CO2DirectionalLinks, default=1000, mutable=True)

		data.load(filename=tab_file_path + '/' + 'CO2_CO2SequestrationNodes.tab', format="set", set=model.CO2SequestrationNodes)
		data.load(filename=tab_file_path + '/' + 'CO2_StorageSiteCapitalCost.tab', format="table", param=model.CO2StorageSiteCapitalCost)
		data.load(filename=tab_file_path + '/' + 'CO2_StorageSiteFixedOMCost.tab', format="table", param=model.StorageSiteFixedOMCost)
		data.load(filename=tab_file_path + '/' + 'CO2_PipelineCapacity.tab', format="table", param=model.CO2PipelineCapacity)
		data.load(filename=tab_file_path + '/' + 'CO2_PipelineCapitalCost.tab', format="table", param=model.CO2PipelineCapCost)
		data.load(filename=tab_file_path + '/' + 'CO2_PipelineFixedOM.tab', format="table", param=model.CO2PipelineOMCost)
		data.load(filename=tab_file_path + '/' + 'CO2_PipelineElectricityUsage.tab', format="table", param=model.CO2PipelineElectricityUsage)

		# data.load(filename=tab_file_path + '/' + 'CO2_LiquefierCapitalCost.tab', format="table", param=model.CO2LiquefierCapitalCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquefierFixedOMCost.tab', format="table", param=model.CO2LiquefierFixedOMCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquefierElectricityUse.tab', format="table", param=model.CO2LiquefierElectricityUse)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidStorageCapitalCost.tab', format="table", param=model.CO2LiquidStorageCapitalCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidStorageFixedOM.tab', format="table", param=model.CO2LiquidStorageFixedOMCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidShipCapitalCost.tab', format="table", param=model.CO2LiquidShipCapitalCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidShipFixedOM.tab', format="table", param=model.CO2LiquidShipFixedOMCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidShipVariableCost.tab', format="table", param=model.CO2LiquidShipVarCost)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidShipCapacity.tab', format="table", param=model.CO2LiquidShipCapacity)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidShipLifetime.tab', format="table", param=model.CO2LiquidShipLifetime)
		# data.load(filename=tab_file_path + '/' + 'CO2_LiquidShipLoadingDischargeTime.tab', format="table", param=model.CO2LiquidShipLoadingDischargeTime)


		model.CO2sentPipeline = Var(model.PipelineTypeOfDirectionalLink, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)
		model.CO2sequestered = Var(model.CO2SequestrationNodes, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals)

		# model.CO2Shipped = Var(model.CO2DirectionalLinks, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
		# model.CO2Liquefied = Var(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
		# model.CO2Regasified = Var(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquidStorageCharge = Var(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquidStorageDischarge = Var(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquidStorageOperational = Var(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, domain=NonNegativeReals, initialize=0)

		model.CO2PipelineBuilt = Var(model.CO2BidirectionalPipelines, model.PipelineType, model.Period, domain=NonNegativeIntegers)
		model.totalCO2PipelineCapacity = Var(model.CO2BidirectionalPipelines, model.PipelineType, model.Period, domain=NonNegativeReals)
		model.CO2SiteCapacityDeveloped = Var(model.CO2SequestrationNodes, model.Period, domain=NonNegativeReals)

		# model.CO2ShipBought = Var(model.CO2DirectionalLinks, model.Period, domain=NonNegativeIntegers, initialize=0)
		# model.totalCO2ShipCapacity = Var(model.CO2DirectionalLinks, model.Period, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquidStorageBuilt = Var(model.OnshoreNodes, model.Period, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquidTotalStorage = Var(model.OnshoreNodes, model.Period, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquefierCapacityBuilt = Var(model.OnshoreNodes, model.Period, domain=NonNegativeReals, initialize=0)
		# model.CO2LiquefierTotalCapacity = Var(model.OnshoreNodes, model.Period, domain=NonNegativeReals, initialize=0)

		def prepCO2InvCosts_rule(model):
			for i in model.Period:
				for (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
					if (n1,n2) in model.CO2BidirectionalPipelines:
						costperyear = (model.WACC / (1 - ((1 + model.WACC) ** (1-model.CO2PipelineLifetime)))) * model.PipelineLength[n1,n2] * model.CO2PipelineCapCost[t] + model.PipelineLength[n1,n2] * model.CO2PipelineOMCost[t]
						costperperiod = costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.CO2PipelineLifetime))))/(1-(1/(1+model.discountrate)))
						model.CO2PipelineInvCost[n1,n2,t,i] = costperperiod

				for n in model.CO2SequestrationNodes:
					#Assume infinite lifetime of storage site -> annual cost approaches WACC * capital cost instead of the previous formula for equivalent annual cost
					costperyear = model.WACC * model.CO2StorageSiteCapitalCost[n] + model.StorageSiteFixedOMCost[n]
					costperperiod = costperyear * (1 - (1 + model.discountrate)**-(value((len(model.Period)-i+1)*5)))/(1-(1/(1+model.discountrate)))
					model.CO2StorageSiteInvCost[n,i] = costperperiod

				# #Investment cost for liquid ship
				# costperyear = (model.WACC / (1 - ((1 + model.WACC) ** (1-model.CO2LiquidShipLifetime)))) * model.CO2LiquidShipCapitalCost + model.CO2LiquidShipFixedOMCost
				# costperperiod = costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.CO2LiquidShipLifetime))))/(1-(1/(1+model.discountrate)))
				# model.CO2LiquidShipInvCost[i] = costperperiod
				#
				# #Investment cost for liquid storage
				# (model.WACC / (1 - ((1 + model.WACC) ** (1-model.CO2LiquidStorageLifetime)))) * model.CO2LiquidStorageCapitalCost + model.CO2LiquidStorageFixedOMCost
				# costperperiod = costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.CO2LiquidStorageCapitalCost))))/(1-(1/(1+model.discountrate)))
				# model.CO2LiquidStorageInvCost[i] = costperperiod
				#
				# #Investment cost for liquefier
				# (model.WACC / (1 - ((1 + model.WACC) ** (1-model.CO2LiquefierLifetime)))) * model.CO2LiquefierCapitalCost + model.CO2LiquefierFixedOMCost
				# costperperiod = costperyear*(1-(1+model.discountrate)**-(min(value((len(model.Period)-i+1)*5), value(model.CO2LiquefierLifetime))))/(1-(1/(1+model.discountrate)))
				# model.CO2LiquefierInvCost[i] = costperperiod

		model.build_CO2InvCosts = BuildAction(rule=prepCO2InvCosts_rule)

		def prep_CO2PipelinePowerUse(model):
			for (n1,n2) in model.CO2BidirectionalPipelines:
					for t in model.PipelineType:
						if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
							# model.CO2PipelinePowerDemandPerTon[n1,n2,t] = model.PipelineLength[n1,n2] * model.CO2PipelineElectricityUsage[t]
							model.CO2PipelinePowerDemandPerTon[n1,n2,t] = model.CO2PipelineElectricityUsage[t]
		model.build_co2_pipeline_power_use = BuildAction(rule=prep_CO2PipelinePowerUse)

		def prep_CO2_storage_cost(model, i):
			return sum(model.CO2StorageSiteInvCost[n,i] * model.CO2SiteCapacityDeveloped[n,i] for n in model.CO2SequestrationNodes)
		model.co2_storage_site_development_cost = Expression(model.Period, rule=prep_CO2_storage_cost)

		def prep_CO2_ship_variable_cost(model):
			for (n1,n2) in model.CO2DirectionalLinks:
				if (n1,n2) in model.CO2BidirectionalPipelines:
					pipeline_length = model.PipelineLength[n1,n2]
				else:
					pipeline_length = model.PipelineLength[n2,n1]
				model.ShippingTimes[n1,n2] = ceil(pipeline_length / model.CO2LiquidShipSpeed + model.CO2LiquidShipLoadingDischargeTime)
				# Linearizing variable cost between n1 and n2
				model.CO2LiquidShipVarUnitCost[n1,n2] = model.CO2LiquidShipVarCost * model.ShippingTimes[n1,n2] / model.CO2LiquidShipCapacity
		# model.build_co2_ship_var_cost = BuildAction(rule=prep_CO2_ship_variable_cost)

		def prep_CO2_ship_operational_expense(model, i):
			operational_expense = sum(model.operationalDiscountrate * model.seasScale[s] * model.sceProbab[w] * model.CO2LiquidShipVarUnitCost[n1,n2] * model.CO2Shipped[n1,n2,h,i,w] for (n1,n2) in model.CO2DirectionalLinks for w in model.Scenario for (s,h) in model.HoursOfSeason)
			return operational_expense
		# model.co2_ship_operational_expense = Expression(model.Period, rule=prep_CO2_ship_operational_expense)


	###############
	##EXPRESSIONS##
	###############

	def multiplier_rule(model,period):
		coeff=1
		if period>1:
			coeff=pow(1.0+model.discountrate,(-5*(int(period)-1)))
		return coeff
	model.discount_multiplier = Expression(model.Period, rule=multiplier_rule)

	def shed_component_rule(model,i):
		hydrogen_shed_cost = 1000*33.3
		load_lost_cost = 0
		load_lost_cost += sum(model.operationalDiscountrate * model.seasScale[s] * model.sceProbab[w] * model.nodeLostLoadCost[n,i]*model.loadShed[n,h,i,w] for n in model.Node for w in model.Scenario for (s,h) in model.HoursOfSeason)
		load_lost_cost += sum(model.operationalDiscountrate *model.seasScale[s] * model.sceProbab[w] * hydrogen_shed_cost * model.h2_demand_shed[n,h,i,w] for n in model.HydrogenProdNode for w in model.Scenario for (s,h) in model.HoursOfSeason)
		return load_lost_cost
	model.shedcomponent = Expression(model.Period, rule=shed_component_rule)

	def operational_cost_rule(model,i):
		return sum(model.operationalDiscountrate*model.seasScale[s]*model.sceProbab[w]*model.genMargCost[g,i]*model.genOperational[n,g,h,i,w] for (n,g) in model.GeneratorsOfNode for (s,h) in model.HoursOfSeason for w in model.Scenario)
	model.operationalcost = Expression(model.Period, rule=operational_cost_rule)

	if hydrogen is True and h2storage is False:
		def hydrogen_operational_storage_cost_rule(model,i):
				electrolyzerOperationalCost = sum(model.operationalDiscountrate*model.seasScale[s]*model.sceProbab[w]*(model.hydrogenProducedElectro[n,h,i,w])*model.averageHydrogenSeasonalStorageCost for n in model.HydrogenProdNode for (s,h) in model.HoursOfSeason for w in model.Scenario)
				reformerOperationalCost = sum(model.operationalDiscountrate*model.seasScale[s]*model.sceProbab[w]*sum(model.hydrogenProducedReformer_ton[n,p,h,i,w] for p in model.ReformerPlants)*model.averageHydrogenSeasonalStorageCost for n in model.ReformerLocations for (s,h) in model.HoursOfSeason for w in model.Scenario)
				return electrolyzerOperationalCost + reformerOperationalCost
		model.hydrogen_operational_storage_cost = Expression(model.Period, rule=hydrogen_operational_storage_cost_rule)

	#############
	##OBJECTIVE##
	#############

	def Obj_rule(model):
		if hydrogen is True:
			returnSum = sum(model.discount_multiplier[i]*(sum(model.genInvCost[g,i]* model.genInvCap[n,g,i] for (n,g) in model.GeneratorsOfNode ) + \
														  sum(model.transmissionVarInvCost[n1,n2,i]*model.transmissionInvCap[n1,n2,i] + model.transmissionFixInvCost[n1,n2,i] * model.transmissionFixedCostPaid[n1,n2,i] for (n1,n2) in model.BidirectionalArc ) + \
														  sum((model.storPWInvCost[b,i]*model.storPWInvCap[n,b,i]+model.storENInvCost[b,i]*model.storENInvCap[n,b,i]) for (n,b) in model.StoragesOfNode ) + \
														  model.shedcomponent[i] + \
														  model.operationalcost[i] + \
														  sum(model.elyzerInvCost[i] * model.elyzerCapBuilt[n,i,i] for n in model.HydrogenProdNode) + \
														  sum(model.hydrogenPipelineInvCost[n1,n2,t] * model.hydrogenPipelineBuilt[n1,n2,t,i] for (n1,n2) in model.HydrogenBidirectionPipelines for t in model.PipelineType if (n1,n2,t) in model.PipelineTypeOfDirectionalLink) + \
														  sum(model.ReformerPlantInvCost[p,i] * model.ReformerCapBuilt[n,p,i] for n in model.ReformerLocations for p in model.ReformerPlants) + \
														  model.reformerOperationalCost[i] + model.reformerEmissionsCost[i] +
														  sum(model.CO2PipelineInvCost[n1,n2,t,i] * model.CO2PipelineBuilt[n1,n2,t,i] for (n1,n2) in model.CO2BidirectionalPipelines for t in model.PipelineType if (n1,n2,t) in model.PipelineTypeOfDirectionalLink) + \
														  model.co2_storage_site_development_cost[i])# + \
														  # sum(model.CO2LiquidShipInvCost[i] * model.CO2ShipBought[n1,n2,i] for (n1,n2) in model.CO2DirectionalLinks) + \
														  # model.co2_ship_operational_expense[i] + \
														  # sum(model.CO2LiquefierInvCost[i] * model.CO2LiquefierCapacityBuilt[n,i] for n in model.OnshoreNodes) + \
														  # sum(model.CO2LiquidStorageInvCost[i] * model.CO2LiquidStorageBuilt[n,i] for n in model.OnshoreNodes))
							for i in model.Period)
			if h2storage is True:
				returnSum += sum(model.discount_multiplier[i]*sum(model.hydrogenStorageBuilt[n,i] * model.hydrogenStorageInvCost[i] for n in model.HydrogenProdNode) for i in model.Period)
			else:
				returnSum += sum(model.discount_multiplier[i] * model.hydrogen_operational_storage_cost[i] for i in model.Period)
			return returnSum
		else:
			return sum(model.discount_multiplier[i] * (sum(model.genInvCost[g, i] * model.genInvCap[n, g, i] for (n, g) in model.GeneratorsOfNode) + \
													   sum(model.transmissionVarInvCost[n1, n2, i] * model.transmissionInvCap[n1, n2, i] for (n1, n2) in model.BidirectionalArc) + \
													   sum((model.storPWInvCost[b, i] * model.storPWInvCap[n, b, i] + model.storENInvCost[b, i] * model.storENInvCap[n, b, i]) for (n, b) in model.StoragesOfNode) + \
													   model.shedcomponent[i] + model.operationalcost[i]) for i in model.Period)
	model.Obj = Objective(rule=Obj_rule, sense=minimize)

	###############
	##CONSTRAINTS##
	###############

	def FlowBalance_rule(model, n, h, i, w):
		if hydrogen is False or n not in model.HydrogenProdNode:
			retVal = sum(model.genOperational[n,g,h,i,w] for g in model.Generator if (n,g) in model.GeneratorsOfNode) \
				   + sum((model.storageDischargeEff[b]*model.storDischarge[n,b,h,i,w]-model.storCharge[n,b,h,i,w]) for b in model.Storage if (n,b) in model.StoragesOfNode) \
				   + sum((model.lineEfficiency[link,n]*model.transmissionOperational[link,n,h,i,w] - model.transmissionOperational[n,link,h,i,w]) for link in model.NodesLinked[n]) \
				   - model.sload[n,h,i,w]
			# retVal += model.loadShed[n,h,i,w]
			return  retVal == 0

		else:
			returnSum = (sum(model.genOperational[n, g, h, i, w] for g in model.Generator if (n, g) in model.GeneratorsOfNode) \
						+ sum((model.storageDischargeEff[b] * model.storDischarge[n, b, h, i, w] - model.storCharge[n, b, h, i, w]) for b in model.Storage if (n, b) in model.StoragesOfNode) \
						+ sum((model.lineEfficiency[link, n] * model.transmissionOperational[link, n, h, i, w] - model.transmissionOperational[n, link, h, i, w]) for link in model.NodesLinked[n]) \
						- model.sload[n, h, i, w] \
						- sum(model.powerForHydrogen[n,j,h,i,w] for j in model.Period if j<=i))
			# returnSum += model.loadShed[n, h, i, w]
			for n2 in model.HydrogenLinks[n]:
				if (n,n2) in model.HydrogenBidirectionPipelines:
					returnSum -= 0.5 * sum(model.hydrogenPipelinePowerDemandPerTon[n,n2,t] * (model.hydrogenSentPipeline[n,n2,t,h,i,w] + model.hydrogenSentPipeline[n2,n,t,h,i,w]) for t in model.PipelineType if (n,n2,t) in model.PipelineTypeOfDirectionalLink)
					returnSum -= 0.5 * sum(model.CO2PipelinePowerDemandPerTon[n,n2,t] * (model.CO2sentPipeline[n,n2,t,h,i,w] + model.CO2sentPipeline[n2,n,t,h,i,w]) for t in model.PipelineType if (n,n2,t) in model.PipelineTypeOfDirectionalLink)
				elif (n2,n) in model.HydrogenBidirectionPipelines:
					returnSum -= 0.5 * sum(model.hydrogenPipelinePowerDemandPerTon[n2,n,t] * (model.hydrogenSentPipeline[n,n2,t,h,i,w] + model.hydrogenSentPipeline[n2,n,t,h,i,w]) for t in model.PipelineType if (n,n2,t) in model.PipelineTypeOfDirectionalLink)
					returnSum -= 0.5 * sum(model.CO2PipelinePowerDemandPerTon[n2,n,t] * (model.CO2sentPipeline[n,n2,t,h,i,w] + model.CO2sentPipeline[n2,n,t,h,i,w]) for t in model.PipelineType if (n,n2,t) in model.PipelineTypeOfDirectionalLink)
			if n in model.ReformerLocations:
				returnSum -= sum(model.ReformerPlantElectricityUse[p,i] * model.hydrogenProducedReformer_ton[n,p,h,i,w] for p in model.ReformerPlants)
			# if n in model.OnshoreNodes:
			# 	returnSum -= model.CO2LiquefierElectricityUse * model.CO2Liquefied[n,h,i,w]
			return returnSum == 0
			#Hydrogen pipeline compressor power usage is split 50/50 between sending node and receiving node
	model.FlowBalance = Constraint(model.Node, model.Operationalhour, model.Period, model.Scenario, rule=FlowBalance_rule)

	#################################################################

	def genMaxProd_rule(model, n, g, h, i, w):
		return model.genOperational[n,g,h,i,w] - model.genCapAvail[n,g,h,w,i]*model.genInstalledCap[n,g,i] <= 0
	model.maxGenProduction = Constraint(model.GeneratorsOfNode, model.Operationalhour, model.Period, model.Scenario, rule=genMaxProd_rule)

	#################################################################

	def ramping_rule(model, n, g, h, i, w):
		if h in model.FirstHoursOfRegSeason or h in model.FirstHoursOfPeakSeason:
			return Constraint.Skip
		else:
			if g in model.ThermalGenerators:
				return model.genOperational[n,g,h,i,w]-model.genOperational[n,g,(h-1),i,w] - model.genRampUpCap[g]*model.genInstalledCap[n,g,i] <= 0   #
			else:
				return Constraint.Skip
	model.ramping = Constraint(model.GeneratorsOfNode, model.Operationalhour, model.Period, model.Scenario, rule=ramping_rule)

	#################################################################

	def storage_energy_balance_rule(model, n, b, h, i, w):
		if h in model.FirstHoursOfRegSeason or h in model.FirstHoursOfPeakSeason:
			return model.storOperationalInit[b]*model.storENInstalledCap[n,b,i] + model.storageChargeEff[b]*model.storCharge[n,b,h,i,w]-model.storDischarge[n,b,h,i,w]-model.storOperational[n,b,h,i,w] == 0   #
		else:
			return model.storageBleedEff[b]*model.storOperational[n,b,(h-1),i,w] + model.storageChargeEff[b]*model.storCharge[n,b,h,i,w]-model.storDischarge[n,b,h,i,w]-model.storOperational[n,b,h,i,w] == 0   #
	model.storage_energy_balance = Constraint(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, rule=storage_energy_balance_rule)

	#################################################################

	def storage_seasonal_net_zero_balance_rule(model, n, b, h, i, w):
		if h in model.FirstHoursOfRegSeason:
			return model.storOperational[n,b,h+value(model.lengthRegSeason)-1,i,w] - model.storOperationalInit[b]*model.storENInstalledCap[n,b,i] == 0  #
		elif h in model.FirstHoursOfPeakSeason:
			return model.storOperational[n,b,h+value(model.lengthPeakSeason)-1,i,w] - model.storOperationalInit[b]*model.storENInstalledCap[n,b,i] == 0  #
		else:
			return Constraint.Skip
	model.storage_seasonal_net_zero_balance = Constraint(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, rule=storage_seasonal_net_zero_balance_rule)

	#################################################################

	def storage_operational_cap_rule(model, n, b, h, i, w):
		return model.storOperational[n,b,h,i,w] - model.storENInstalledCap[n,b,i]  <= 0   #
	model.storage_operational_cap = Constraint(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, rule=storage_operational_cap_rule)

	#################################################################

	def storage_power_discharg_cap_rule(model, n, b, h, i, w):
		return model.storDischarge[n,b,h,i,w] - model.storageDiscToCharRatio[b]*model.storPWInstalledCap[n,b,i] <= 0   #
	model.storage_power_discharg_cap = Constraint(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, rule=storage_power_discharg_cap_rule)

	#################################################################

	def storage_power_charg_cap_rule(model, n, b, h, i, w):
		return model.storCharge[n,b,h,i,w] - model.storPWInstalledCap[n,b,i] <= 0   #
	model.storage_power_charg_cap = Constraint(model.StoragesOfNode, model.Operationalhour, model.Period, model.Scenario, rule=storage_power_charg_cap_rule)

	#################################################################

	def hydro_gen_limit_rule(model, n, g, s, i, w):
		if g in model.RegHydroGenerator:
			return sum(model.genOperational[n,g,h,i,w] for h in model.Operationalhour if (s,h) in model.HoursOfSeason) - model.maxRegHydroGen[n,i,s,w] <= 0
		else:
			return Constraint.Skip  #
	model.hydro_gen_limit = Constraint(model.GeneratorsOfNode, model.Season, model.Period, model.Scenario, rule=hydro_gen_limit_rule)

	#################################################################

	def hydro_node_limit_rule(model, n, i):
		return sum(model.genOperational[n,g,h,i,w]*model.seasScale[s]*model.sceProbab[w] for g in model.HydroGenerator if (n,g) in model.GeneratorsOfNode for (s,h) in model.HoursOfSeason for w in model.Scenario) - model.maxHydroNode[n] <= 0   #
	model.hydro_node_limit = Constraint(model.Node, model.Period, rule=hydro_node_limit_rule)


	#################################################################

	def transmission_cap_rule(model, n1, n2, h, i, w):
		if (n1,n2) in model.BidirectionalArc:
			return model.transmissionOperational[(n1,n2),h,i,w]  - model.transmissionInstalledCap[(n1,n2),i] <= 0
		elif (n2,n1) in model.BidirectionalArc:
			return model.transmissionOperational[(n1,n2),h,i,w]  - model.transmissionInstalledCap[(n2,n1),i] <= 0
	model.transmission_cap = Constraint(model.DirectionalLink, model.Operationalhour, model.Period, model.Scenario, rule=transmission_cap_rule)

	def transmission_pay_fixed_cost(model, n1, n2, i):
		for t in model.TransmissionType:
			if (n1,n2,t) in model.TransmissionTypeOfDirectionalLink:
				if "hvdc" in t.lower():
					return model.transmissionInvCap[n1,n2,i] <= 40000 * model.transmissionFixedCostPaid[n1,n2,i]
				else:
					#If it is not HVDC, then there is no fixed cost, and we just force the binary variable to 0.
					return model.transmissionFixedCostPaid[n1,n2,i] == 0
	model.tranmission_pay_fixed_cost = Constraint(model.BidirectionalArc, model.Period, rule=transmission_pay_fixed_cost)

	if windfarmNodes is not None:
		#This constraints restricts the transmission through offshore wind farms, so that the total transmission capacity cannot be bigger than the invested generation capacity
		# def wind_farm_tranmission_cap_rule(model, n, i):
		# 	sumCap = 0
		# 	for n2 in model.NodesLinked[n]:
		# 		if (n,n2) in model.BidirectionalArc:
		# 			sumCap += model.transmissionInstalledCap[(n,n2),i]
		# 		else:
		# 			sumCap += model.transmissionInstalledCap[(n2,n),i]
		# 	return sumCap <= sum(model.genInstalledCap[n,g,i] for g in model.Generator if (n,g) in model.GeneratorsOfNode)
		# model.wind_farm_transmission_cap = Constraint(model.windfarmNodes, model.Period, rule=wind_farm_tranmission_cap_rule)
		def wind_farm_tranmission_cap_rule(model, n1, n2, i):
			if n1 in model.windfarmNodes or n2 in model.windfarmNodes:
				if (n1,n2) in model.BidirectionalArc:
					if n1 in model.windfarmNodes:
						return model.transmissionInstalledCap[(n1,n2),i] <= sum(model.genInstalledCap[n1,g,i] for g in model.Generator if (n1,g) in model.GeneratorsOfNode)
					else:
						return model.transmissionInstalledCap[(n1,n2),i] <= sum(model.genInstalledCap[n2,g,i] for g in model.Generator if (n2,g) in model.GeneratorsOfNode)
				elif (n2,n1) in model.BidirectionalArc:
					if n1 in model.windfarmNodes:
						return model.transmissionInstalledCap[(n2,n1),i] <= sum(model.genInstalledCap[n1,g,i] for g in model.Generator if (n1,g) in model.GeneratorsOfNode)
					else:
						return model.transmissionInstalledCap[(n2,n1),i] <= sum(model.genInstalledCap[n2,g,i] for g in model.Generator if (n2,g) in model.GeneratorsOfNode)
				else:
					return Constraint.Skip
			else:
				return Constraint.Skip
		model.wind_farm_transmission_cap = Constraint(model.Node, model.Node, model.Period, rule=wind_farm_tranmission_cap_rule)



	#################################################################

	if EMISSION_CAP:
		def emission_cap_rule(model, i, w):
			return sum(model.seasScale[s]*model.genCO2TypeFactor[g]*(3.6/model.genEfficiency[g,i])*model.genOperational[n,g,h,i,w] for (n,g) in model.GeneratorsOfNode for (s,h) in model.HoursOfSeason)/1000000 \
				   - model.CO2cap[i] <= 0   #
		model.emission_cap = Constraint(model.Period, model.Scenario, rule=emission_cap_rule)

	#################################################################

	def lifetime_rule_gen(model, n, g, i):
		startPeriod=1
		if value(1+i-(model.genLifetime[g]/model.LeapYearsInvestment))>startPeriod:
			startPeriod=value(1+i-model.genLifetime[g]/model.LeapYearsInvestment)
		return sum(model.genInvCap[n,g,j]  for j in model.Period if j>=startPeriod and j<=i )- model.genInstalledCap[n,g,i] + model.genInitCap[n,g,i]== 0   #
	model.installedCapDefinitionGen = Constraint(model.GeneratorsOfNode, model.Period, rule=lifetime_rule_gen)

	#################################################################

	def lifetime_rule_storEN(model, n, b, i):
		startPeriod=1
		if value(1+i-model.storageLifetime[b]*(1/model.LeapYearsInvestment))>startPeriod:
			startPeriod=value(1+i-model.storageLifetime[b]/model.LeapYearsInvestment)
		return sum(model.storENInvCap[n,b,j]  for j in model.Period if j>=startPeriod and j<=i )- model.storENInstalledCap[n,b,i] + model.storENInitCap[n,b,i]== 0   #
	model.installedCapDefinitionStorEN = Constraint(model.StoragesOfNode, model.Period, rule=lifetime_rule_storEN)

	#################################################################

	def lifetime_rule_storPOW(model, n, b, i):
		startPeriod=1
		if value(1+i-model.storageLifetime[b]*(1/model.LeapYearsInvestment))>startPeriod:
			startPeriod=value(1+i-model.storageLifetime[b]/model.LeapYearsInvestment)
		return sum(model.storPWInvCap[n,b,j]  for j in model.Period if j>=startPeriod and j<=i )- model.storPWInstalledCap[n,b,i] + model.storPWInitCap[n,b,i]== 0   #
	model.installedCapDefinitionStorPOW = Constraint(model.StoragesOfNode, model.Period, rule=lifetime_rule_storPOW)

	#################################################################

	def lifetime_rule_trans(model, n1, n2, i):
		startPeriod=1
		if value(1+i-model.transmissionLifetime[n1,n2]*(1/model.LeapYearsInvestment))>startPeriod:
			startPeriod=value(1+i-model.transmissionLifetime[n1,n2]/model.LeapYearsInvestment)
		return sum(model.transmissionInvCap[n1,n2,j]  for j in model.Period if j>=startPeriod and j<=i )- model.transmissionInstalledCap[n1,n2,i] + model.transmissionInitCap[n1,n2,i] == 0   #
	model.installedCapDefinitionTrans = Constraint(model.BidirectionalArc, model.Period, rule=lifetime_rule_trans)

	#################################################################

	def investment_gen_cap_rule(model, t, n, i):
		return sum(model.genInvCap[n,g,i] for g in model.Generator if (n,g) in model.GeneratorsOfNode and (t,g) in model.GeneratorsOfTechnology) - model.genMaxBuiltCap[n,t,i] <= 0
	model.investment_gen_cap = Constraint(model.Technology, model.Node, model.Period, rule=investment_gen_cap_rule)

	#################################################################

	def investment_trans_cap_rule(model, n1, n2, i):
		return model.transmissionInvCap[n1,n2,i] - model.transmissionMaxBuiltCap[n1,n2,i] <= 0
	model.investment_trans_cap = Constraint(model.BidirectionalArc, model.Period, rule=investment_trans_cap_rule)

	#################################################################

	def investment_storage_power_cap_rule(model, n, b, i):
		return model.storPWInvCap[n,b,i] - model.storPWMaxBuiltCap[n,b,i] <= 0
	model.investment_storage_power_cap = Constraint(model.StoragesOfNode, model.Period, rule=investment_storage_power_cap_rule)

	#################################################################

	def investment_storage_energy_cap_rule(model, n, b, i):
		return model.storENInvCap[n,b,i] - model.storENMaxBuiltCap[n,b,i] <= 0
	model.investment_storage_energy_cap = Constraint(model.StoragesOfNode, model.Period, rule=investment_storage_energy_cap_rule)

	################################################################

	def installed_gen_cap_rule(model, t, n, i):
		return sum(model.genInstalledCap[n,g,i] for g in model.Generator if (n,g) in model.GeneratorsOfNode and (t,g) in model.GeneratorsOfTechnology) - model.genMaxInstalledCap[n,t,i] <= 0
	model.installed_gen_cap = Constraint(model.Technology, model.Node, model.Period, rule=installed_gen_cap_rule)

	#################################################################

	def installed_trans_cap_rule(model, n1, n2, i):
		return model.transmissionInstalledCap[n1, n2, i] - model.transmissionMaxInstalledCap[n1, n2, i] <= 0
	model.installed_trans_cap = Constraint(model.BidirectionalArc, model.Period, rule=installed_trans_cap_rule)

	#################################################################

	def installed_storage_power_cap_rule(model, n, b, i):
		return model.storPWInstalledCap[n,b,i] - model.storPWMaxInstalledCap[n,b,i] <= 0
	model.installed_storage_power_cap = Constraint(model.StoragesOfNode, model.Period, rule=installed_storage_power_cap_rule)

	#################################################################

	def installed_storage_energy_cap_rule(model, n, b, i):
		return model.storENInstalledCap[n,b,i] - model.storENMaxInstalledCap[n,b,i] <= 0
	model.installed_storage_energy_cap = Constraint(model.StoragesOfNode, model.Period, rule=installed_storage_energy_cap_rule)

	#################################################################

	def power_energy_relate_rule(model, n, b, i):
		if b in model.DependentStorage:
			return model.storPWInstalledCap[n,b,i] - model.storagePowToEnergy[b]*model.storENInstalledCap[n,b,i] == 0   #
		else:
			return Constraint.Skip
	model.power_energy_relate = Constraint(model.StoragesOfNode, model.Period, rule=power_energy_relate_rule)

	#################################################################

	# To any reader of this code, the next two constraints are very ugly, and there is likely a better implementation that achieves the same. They were put together as quick fixes, and will be fixed if I remember, have time and can be bothered (in that order of priority). The last is most likely to fail.
	def powerFromHydrogenRule(model, n, g, h, i, w):
		if hydrogen is True:
			if g in model.HydrogenGenerators:
				if (n,g) in model.GeneratorsOfNode:
					return model.genOperational[n,g,h,i,w] == model.genEfficiency[g,i] * model.hydrogenForPower[g,n,h,i,w] * model.hydrogenLHV_kg * 1000
				else:
					return model.genOperational[n,g,h,i,w] == 0
			else:
				return Constraint.Skip
		else:
			if g in model.HydrogenGenerators:
				return model.genOperational[n,g,h,i,w] == 0
			else:
				return Constraint.Skip
	model.powerFromHydrogen = Constraint(model.GeneratorsOfNode, model.Operationalhour, model.Period, model.Scenario, rule=powerFromHydrogenRule)

	# Commented out this constraint, because we initialize model.hydrogenForPower to 0 in the variable definition.
	# Early testing suggests that this is enough (all generators other than hydrogen generators have their associated
	# hydrogenForPower set to 0. Will not fully delete the constraint (yet) in case other results suggest we need to
	# explicitly limit hydrogenForPower for non-hydrogen generators.

	# def hydrogenToGenerator_rule(model,n,g,h,i,w):
	# 	if n in model.HydrogenProdNode and g not in model.HydrogenGenerators:
	# 		return model.hydrogenForPower[g,n,h,i,w]==0
	# 	else:
	# 		return Constraint.Skip
	# model.hydrogenToGenerator = Constraint(model.GeneratorsOfNode,model.Operationalhour,model.Period,model.Scenario, rule=hydrogenToGenerator_rule)

	# Hydrogen constraints
	if hydrogen is True:
		def lifetime_rule_pipeline(model,n1,n2,t,i):
			if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
				startPeriod = 1
				if value(1+i-model.hydrogenPipelineLifetime/model.LeapYearsInvestment)>startPeriod:
					startPeriod=value(1+i-model.hydrogenPipelineLifetime/model.LeapYearsInvestment)
				return sum(model.hydrogenPipelineCapacity[t] * model.hydrogenPipelineBuilt[n1,n2,t,j] for j in model.Period if j>=startPeriod and j<=i) - model.totalHydrogenPipelineCapacity[n1,n2,t,i] == 0
			else:
				return Constraint.Skip
		model.installedCapDefinitionPipe = Constraint(model.HydrogenBidirectionPipelines, model.PipelineType, model.Period, rule=lifetime_rule_pipeline)

		def lifetime_rule_elyzer(model,n,j,i): # j = year the electrolyzer is built, i = current year when the capacity is checked (i.e., has the electrolyzer worn out?)
			startPeriod = 1
			if value(1+i-model.elyzerLifetime/model.LeapYearsInvestment)>startPeriod:
				startPeriod=value(1+i-model.elyzerLifetime/model.LeapYearsInvestment)
			return sum(model.elyzerCapBuilt[n,j,k] for k in model.Period if k>=startPeriod and k<=i) - model.elyzerTotalCap[n,j,i] == 0
		model.installedCapDefinitionElyzer = Constraint(model.HydrogenProdNode,model.Period, model.Period, rule=lifetime_rule_elyzer)

		def link_electrolyzer_year_built_and_period_rule(model,n,i,j):
			if i != j:
				return model.elyzerCapBuilt[n,j,i] == 0
			else:
				return Constraint.Skip
		model.link_electrolyzer_year_built_and_period = Constraint(model.HydrogenProdNode, model.Period, model.Period, rule=link_electrolyzer_year_built_and_period_rule)

		def lifetime_rule_reformer(model,n,p,i):
			startPeriod = 1
			if value(1+i-model.ReformerPlantLifetime[p]/model.LeapYearsInvestment)>startPeriod:
				startPeriod = value(1+i-model.ReformerPlantLifetime[p]/model.LeapYearsInvestment)
			return sum(model.ReformerCapBuilt[n,p,j] for j in model.Period if j>=startPeriod and j<=i) - model.ReformerTotalCap[n,p,i] == 0
		model.installedCapDefinitionReformer = Constraint(model.ReformerLocations, model.ReformerPlants, model.Period, rule=lifetime_rule_reformer)

		def pipeline_cap_rule(model,n1,n2,t,h,i,w):
			if (n1,n2) in model.HydrogenBidirectionPipelines:
				if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
					return model.hydrogenSentPipeline[(n1,n2),t,h,i,w] - model.totalHydrogenPipelineCapacity[(n1,n2),t,i] <= 0
				else:
					return Constraint.Skip
			elif (n2,n1) in model.HydrogenBidirectionPipelines:
				if (n2,n1,t) in model.PipelineTypeOfDirectionalLink:
					return model.hydrogenSentPipeline[(n1,n2),t,h,i,w] - model.totalHydrogenPipelineCapacity[(n2,n1),t,i] <= 0
				else:
					return Constraint.Skip
			else:
				print('Problem creating max pipeline capacity constraint for nodes ' + n1 +' and ' + n2)
				exit()
		model.pipeline_cap = Constraint(model.AllowedHydrogenLinks, model.PipelineType, model.Operationalhour, model.Period, model.Scenario, rule=pipeline_cap_rule)

		def hydrogen_flow_balance_rule(model,n,h,i,w):
			balance = - model.meetHydrogenDemand[n,h,i,w] \
					  - sum(model.hydrogenForPower[g,n,h,i,w] for g in model.HydrogenGenerators) \
					  + sum(model.hydrogenSentPipeline[(n2,n), t, h, i, w] - model.hydrogenSentPipeline[(n,n2), t, h, i, w] for n2 in model.HydrogenLinks[n] for t in model.PipelineType if (n,n2,t) in model.PipelineTypeOfDirectionalLink) \
					  # + model.h2_demand_shed[n,h,i,w]
			if n in model.HydrogenProdNode:
				balance += model.hydrogenProducedElectro[n,h,i,w]
			if n in model.ReformerLocations:
				balance += sum(model.hydrogenProducedReformer_ton[n,p,h,i,w] for p in model.ReformerPlants)
			if h2storage is True:
				balance += model.hydrogenDischargeStorage[n,h,i,w] - model.hydrogenChargeStorage[n,h,i,w]
			return balance == 0
		model.hydrogen_flow_balance = Constraint(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_flow_balance_rule)

		def meet_hydrogen_demand_rule(model,n,h,i,w):
			return model.meetHydrogenDemand[n,h,i,w] == model.hydrogenDemand[n,i,w]
		model.meet_hydrogen_demand = Constraint(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, rule=meet_hydrogen_demand_rule)

		def hydrogen_production_rule(model,n,h,i,w):
			return model.hydrogenProducedElectro[n,h,i,w] == sum(model.powerForHydrogen[n,j,h,i,w] / model.elyzerPowerConsumptionPerTon[j] for j in model.Period if j<=i)
		model.hydrogen_production = Constraint(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_production_rule)

		def hydrogen_production_electrolyzer_capacity_rule(model,n,j,h,i,w):
			return model.powerForHydrogen[n,j,h,i,w] <= model.elyzerTotalCap[n,j,i] #j = Year the electrolyzer was bought
		model.hydrogen_production_electrolyzer_capacity = Constraint(model.HydrogenProdNode, model.Period, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_production_electrolyzer_capacity_rule)

		def hydrogen_production_reformer_capacity_rule(model,n,p,h,i,w):
			return model.hydrogenProducedReformer_MWh[n,p,h,i,w] <= model.ReformerTotalCap[n,p,i]
		model.hydrogen_production_reformer_capacity = Constraint(model.ReformerLocations, model.ReformerPlants, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_production_reformer_capacity_rule)

		def hydrogen_link_reformer_kg_MWh_rule(model,n,p,h,i,w):
			return model.hydrogenProducedReformer_ton[n,p,h,i,w] == model.hydrogenProducedReformer_MWh[n,p,h,i,w]/(model.hydrogenLHV_kg * 1000)
		model.hydrogen_link_reformer_kg_MWh = Constraint(model.ReformerLocations, model.ReformerPlants, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_link_reformer_kg_MWh_rule)

		if h2storage is True:
			def hydrogen_storage_balance_rule(model,n,h,i,w):
				if h in model.FirstHoursOfRegSeason or h in model.FirstHoursOfPeakSeason:
					return model.hydrogenStorageInitOperational * model.hydrogenTotalStorage[n,i]+model.hydrogenChargeStorage[n,h,i,w]-model.hydrogenDischargeStorage[n,h,i,w] - model.hydrogenStorageOperational[n,h,i,w] == 0
				else:
					return model.hydrogenStorageOperational[n,h-1,i,w] + model.hydrogenChargeStorage[n,h,i,w] - model.hydrogenDischargeStorage[n,h,i,w] - model.hydrogenStorageOperational[n,h,i,w] == 0
			model.hydrogen_storage_balance = Constraint(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_storage_balance_rule)

			def hydrogen_storage_max_capacity_rule(model,n,i):
				return model.hydrogenTotalStorage[n,i] <= model.hydrogenMaxStorageCapacity[n]
			model.hydrogen_storage_max_capacity = Constraint(model.HydrogenProdNode, model.Period, rule=hydrogen_storage_max_capacity_rule)

			def hydrogen_storage_operational_capacity_rule(model,n,h,i,w):
				return model.hydrogenStorageOperational[n,h,i,w] <= model.hydrogenTotalStorage[n,i]
			model.hydrogen_storage_operational_capacity = Constraint(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_storage_operational_capacity_rule)

			def hydrogen_balance_storage_rule(model,n,h,i,w):
				if h in model.FirstHoursOfRegSeason:
					return model.hydrogenStorageOperational[n,h+value(model.lengthRegSeason)-1,i,w] - model.hydrogenStorageInitOperational * model.hydrogenTotalStorage[n,i] == 0
				elif h in model.FirstHoursOfPeakSeason:
					return model.hydrogenStorageOperational[n,h+value(model.lengthPeakSeason)-1,i,w] - model.hydrogenStorageInitOperational * model.hydrogenTotalStorage[n,i] == 0
				else:
					return Constraint.Skip
			model.hydrogen_balance_storage = Constraint(model.HydrogenProdNode, model.Operationalhour, model.Period, model.Scenario, rule=hydrogen_balance_storage_rule)

			def hydrogen_storage_lifetime_rule(model,n,i):
				startPeriod=1
				if value(1+i-model.hydrogenStorageLifetime/model.LeapYearsInvestment)>startPeriod:
					startPeriod = value(1+i-model.hydrogenStorageLifetime/model.LeapYearsInvestment)
				return sum(model.hydrogenStorageBuilt[n,j] for j in model.Period if j>=startPeriod and j <=i) - model.hydrogenTotalStorage[n,i] == 0
			model.hydrogen_storage_lifetime = Constraint(model.HydrogenProdNode, model.Period, rule=hydrogen_storage_lifetime_rule)

		#CO2 constraints

		def co2_pipeline_lifetime_rule(model,n1,n2,t,i):
			if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
				startPeriod = 1
				if value(1+i-model.CO2PipelineLifetime/model.LeapYearsInvestment)>startPeriod:
					startPeriod=value(1+i-model.CO2PipelineLifetime/model.LeapYearsInvestment)
				return sum(model.CO2PipelineCapacity[t] * model.CO2PipelineBuilt[n1,n2,t,j] for j in model.Period if j>=startPeriod and j<=i) - model.totalCO2PipelineCapacity[n1,n2,t,i] == 0
			else:
				return Constraint.Skip
		model.co2_pipeline_lifetime = Constraint(model.CO2BidirectionalPipelines, model.PipelineType, model.Period, rule=co2_pipeline_lifetime_rule)

		def co2_pipeline_cap_rule(model, n1, n2, t, h, i, w):
			if (n1,n2) in model.CO2BidirectionalPipelines:
				if (n1,n2,t) in model.PipelineTypeOfDirectionalLink:
					return model.CO2sentPipeline[(n1,n2),t,h,i,w] - model.totalCO2PipelineCapacity[(n1,n2),t,i] <= 0
				else:
					return Constraint.Skip
			elif (n2,n1) in model.CO2BidirectionalPipelines:
				if (n2,n1,t) in model.PipelineTypeOfDirectionalLink:
					return model.CO2sentPipeline[(n1,n2),t,h,i,w] - model.totalCO2PipelineCapacity[(n2,n1),t,i] <= 0
				else:
					return Constraint.Skip
		model.co2_pipeline_cap = Constraint(model.CO2DirectionalLinks, model.PipelineType, model.Operationalhour, model.Period, model.Scenario, rule=co2_pipeline_cap_rule)

		def co2_flow_balance_rule(model, n, h, i, w):
			balance = 0
			balance += model.co2_captured_generators[n,h,i,w]
			balance += model.co2_captured_reformers[n,h,i,w]
			balance += sum(model.CO2sentPipeline[n2,n,t,h,i,w] - model.CO2sentPipeline[n,n2,t,h,i,w] for n2 in model.OnshoreNodes if (n,n2) in model.CO2DirectionalLinks for t in model.PipelineType if (n,n2,t) in model.PipelineTypeOfDirectionalLink)
			# balance += model.CO2Regasified[n,h,i,w] - model.CO2Liquefied[n,h,i,w]
			if n in model.CO2SequestrationNodes:
				balance -= model.CO2sequestered[n,h,i,w]
			return balance == 0
		model.co2_flow_balance = Constraint(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, rule=co2_flow_balance_rule)

		def co2_sequestering_max_capacity_rule(model, n, h, i, w):
			return model.CO2sequestered[n,h,i,w] <= sum(model.CO2SiteCapacityDeveloped[n,j] for j in model.Period if j<=i)
		model.co2_sequestering_max_capacity = Constraint(model.CO2SequestrationNodes, model.Operationalhour, model.Period, model.Scenario, rule=co2_sequestering_max_capacity_rule)

		def co2_ship_lifetime_rule(model,n1,n2,i):
			startPeriod = 1
			if value(1 + i - model.CO2LiquidShipLifetime / model.LeapYearsInvestment) > startPeriod:
				startPeriod = value(1 + i - model.CO2LiquidShipLifetime / model.LeapYearsInvestment)
			return sum(model.CO2LiquidShipCapacity * model.CO2ShipBought[n1,n2,j] for j in model.Period if j>=startPeriod and j<=i) - model.totalCO2ShipCapacity[n1,n2,i] == 0
		# model.co2_ship_lifetime = Constraint(model.CO2DirectionalLinks, model.Period, rule=co2_ship_lifetime_rule)

		def co2_ship_capacity_rule(model, n1, n2, h, i, w):
			if h < min(model.FirstHoursOfPeakSeason): # We are in normal seasons
				season_starting_hour = 1 + floor((h-1)/lengthRegSeason) * lengthRegSeason
				season_ending_hour = season_starting_hour + lengthRegSeason - 1
				if value(model.ShippingTimes[n1,n2]) >= lengthRegSeason: # If the shipping time is longer than the representative period, then we can only ship one full ship in the investigated time
					return sum(model.CO2Shipped[n1,n2,h_cand,i,w] for h_cand in model.Operationalhour if h_cand >= season_starting_hour and h_cand <= season_ending_hour) - model.totalCO2ShipCapacity[n1,n2,i] <= 0
				else: # We can ship several times in the investigated season
					if h - season_starting_hour + 1 > value(model.ShippingTimes[n1,n2]): # If the hour is high enough, then we won't loop back to the end of the representative season
						startHour = h - value(model.ShippingTimes[n1,n2]) + 1
						return sum(model.CO2Shipped[n1,n2,h_candidate,i,w] for h_candidate in model.Operationalhour if  h_candidate >= startHour and h_candidate <= h) - model.totalCO2ShipCapacity[n1,n2,i] <= 0
					else: # We need to consider the ending hours of this representative season
						overshooting_hours = value(model.ShippingTimes[n1,n2]) + season_starting_hour - h
						end_of_season_index = overshooting_hours % lengthRegSeason
						return sum(model.CO2Shipped[n1,n2,h_candidate,i,w] for h_candidate in model.Operationalhour if  h_candidate >= end_of_season_index and h_candidate <= season_ending_hour) + sum(model.CO2Shipped[n1,n2,h_candidate,i,w] for h_candidate in model.Operationalhour if  h_candidate >= season_starting_hour and h_candidate <= h) - model.totalCO2ShipCapacity[n1,n2,i] <= 0
			else: # We are in peak seasons
				season_starting_hour = 1 + floor((h - 1 - NoOfRegSeason * lengthRegSeason)/lengthPeakSeason) * lengthPeakSeason + NoOfRegSeason * lengthRegSeason
				season_ending_hour = season_starting_hour + lengthPeakSeason - 1
				if value(model.ShippingTimes[n1,n2]) >= lengthPeakSeason:
					return sum(model.CO2Shipped[n1,n2,h_cand,i,w] for h_cand in model.Operationalhour if h_cand >= season_starting_hour and h_cand <= season_ending_hour) - model.totalCO2ShipCapacity[n1,n2,i] <= 0
				else:
					if h - season_starting_hour + 1 > value(model.ShippingTimes[n1,n2]):
						startHour = h - value(model.ShippingTimes[n1,n2]) + 1
						return sum(model.CO2Shipped[n1,n2,h_candidate,i,w] for h_candidate in model.Operationalhour if  h_candidate >= startHour and h_candidate <= h) - model.totalCO2ShipCapacity[n1,n2,i] <= 0
					else:
						overshooting_hours = value(model.ShippingTimes[n1,n2]) + season_starting_hour - h
						end_of_season_index = overshooting_hours % lengthPeakSeason
						return sum(model.CO2Shipped[n1,n2,h_candidate,i,w] for h_candidate in model.Operationalhour if  h_candidate >= end_of_season_index and h_candidate <= season_ending_hour) + sum(model.CO2Shipped[n1,n2,h_candidate,i,w] for h_candidate in model.Operationalhour if  h_candidate >= season_starting_hour and h_candidate <= h) - model.totalCO2ShipCapacity[n1,n2,i] <= 0
		# model.co2_ship_capacity = Constraint(model.CO2DirectionalLinks, model.Operationalhour, model.Period, model.Scenario, rule=co2_ship_capacity_rule)

		def co2_liquefier_lifetime_rule(model, n, i):
			startPeriod = 1
			if value(1 + i - model.CO2LiquefierLifetime / model.LeapYearsInvestment) > startPeriod:
				startPeriod = value(1 + i - model.CO2LiquefierLifetime / model.LeapYearsInvestment)
			return sum(model.CO2LiquefierCapacityBuilt[n,j] for j in model.Period if j >= startPeriod and j <= i) - model.CO2LiquefierTotalCapacity[n,i] == 0
		# model.co2_liquefier_lifetime = Constraint(model.OnshoreNodes, model.Period, rule=co2_liquefier_lifetime_rule)

		def co2_liquefier_capacity_rule(model, n, h, i, w):
			return model.CO2Liquefied[n,h,i,w] <= model.CO2LiquefierTotalCapacity[n,i]
		# model.co2_liquefier_capacity = Constraint(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, rule=co2_liquefier_capacity_rule)

		def liquid_co2_storage_balance_rule(model, n, h, i, w):
			if h in model.FirstHoursOfRegSeason or h in model.FirstHoursOfPeakSeason:
				return model.CO2LiquidStorageInitOperational * model.CO2LiquidTotalStorage[n,i] + model.CO2LiquidStorageCharge[n,h,i,w] - model.CO2LiquidStorageDischarge[n,h,i,w] == model.CO2LiquidStorageOperational[n,h,i,w]
			else:
				return model.CO2LiquidStorageOperational[n,h-1,i,w] + model.CO2LiquidStorageCharge[n,h,i,w] - model.CO2LiquidStorageDischarge[n,h,i,w] == model.CO2LiquidStorageOperational[n,h,i,w]
		# model.liquid_co2_storage_balance = Constraint(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, rule=liquid_co2_storage_balance_rule)

		def liquid_co2_storage_operational_capacity_rule(model, n, h, i, w):
			return model.CO2LiquidStorageOperational[n,h,i,w] <= model.CO2LiquidTotalStorage[n,i]
		# model.liquid_co2_storage_operational_capacity = Constraint(model.OnshoreNodes, model.Operationalhour ,model.Period, model.Scenario, rule=liquid_co2_storage_operational_capacity_rule)

		def liquid_co2_storage_net_zero_rule(model,n,h,i,w):
			if h in model.FirstHoursOfRegSeason:
				return model.CO2LiquidStorageOperational[n, h+lengthRegSeason-1,i,w] - model.CO2LiquidStorageInitOperational * model.CO2LiquidTotalStorage[n,i] == 0
			elif h in model.FirstHoursOfPeakSeason:
				return model.CO2LiquidStorageOperational[n, h+lengthPeakSeason-1,i,w] - model.CO2LiquidStorageInitOperational * model.CO2LiquidTotalStorage[n,i] == 0
			else:
				return Constraint.Skip
		# model.liquid_co2_storage_net_zero = Constraint(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, rule=liquid_co2_storage_net_zero_rule)

		def liquid_co2_storage_lifetime_rule(model,n,i):
			startPeriod = 1
			if value(1+i-model.CO2LiquidStorageLifetime/model.LeapYearsInvestment)>startPeriod:
				startPeriod = value(1+i-model.CO2LiquidStorageLifetime/model.LeapYearsInvestment)
			return sum(model.CO2LiquidStorageBuilt[n,j] for j in model.Period if j>=startPeriod and j <=i) - model.CO2LiquidTotalStorage[n,i] == 0
		# model.liquid_co2_storage_lifetime = Constraint(model.OnshoreNodes, model.Period, rule=liquid_co2_storage_lifetime_rule)

		def liquid_co2_flow_balance_rule(model, n, h, i, w):
			balance = 0
			balance += model.CO2Liquefied[n,h,i,w] - model.CO2Regasified[n,h,i,w]
			balance += model.CO2LiquidStorageDischarge[n,h,i,w] - model.CO2LiquidStorageCharge[n,h,i,w]
			for n2 in model.OnshoreNodes:
				if (n,n2) in model.CO2DirectionalLinks:
					half_shipping_time = ceil(value(model.ShippingTimes[n,n2]) / 2)
					if h < min(model.FirstHoursOfPeakSeason): # We are in normal seasons
						season_starting_hour = 1 + floor((h-1)/lengthRegSeason) * lengthRegSeason
						starting_hour_receiving_ship = h - half_shipping_time
						while starting_hour_receiving_ship < season_starting_hour:
							starting_hour_receiving_ship += lengthRegSeason
					else: # We are in peak seasons
						season_starting_hour = 1 + floor((h - 1 - NoOfRegSeason * lengthRegSeason)/lengthPeakSeason) * lengthPeakSeason + NoOfRegSeason * lengthRegSeason
						starting_hour_receiving_ship = h-half_shipping_time
						while starting_hour_receiving_ship < season_starting_hour:
							starting_hour_receiving_ship += lengthPeakSeason
					balance -= model.CO2Shipped[n,n2,h,i,w] - model.CO2Shipped[n2,n,starting_hour_receiving_ship,i,w]
			return balance == 0
		# model.liquid_co2_flow_balance = Constraint(model.OnshoreNodes, model.Operationalhour, model.Period, model.Scenario, rule=liquid_co2_flow_balance_rule)


	stopConstraints = startBuild = datetime.now()
	#######
	##RUN##
	#######

	print("Objective and constraints read...")

	print("{hour}:{minute}:{second}: Building instance...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))

	start = time.time()

	instance = model.create_instance(data) #, report_timing=True)
	# instance.dual = Suffix(direction=Suffix.IMPORT) #Make sure the dual value is collected into solver results (if solver supplies dual information)
	# instance.seasScale.pprint()

	inv_per = []
	for i in instance.Period:
		my_string = str(value(2015+int(i)*5))+"-"+str(value(2020+int(i)*5))
		inv_per.append(my_string)

	print("{hour}:{minute}:{second}: Writing load data to data_electric_load.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(tab_file_path + "/" + 'data_electric_load.csv', 'w', newline='')
	writer = csv.writer(f)
	my_header = ["Node","Period","Scenario","Season","Hour",'Electric load [MW]']
	writer.writerow(my_header)
	for n in instance.Node:
		for i in instance.Period:
			for w in instance.Scenario:
				for (s,h) in instance.HoursOfSeason:
					my_string = [n,inv_per[int(i-1)],w,s,h,
								 value(instance.sload[n,h,i,w])]
					writer.writerow(my_string)
	f.close()

	end = time.time()
	print("{hour}:{minute}:{second}: Building instance took:".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")) + str(end - start))

	endBuild = startOptimization = datetime.now()

	#import pdb; pdb.set_trace()

	print("\n----------------------Problem Statistics---------------------")
	print("Nodes: "+ str(len(instance.Node)))
	print("Lines: "+str(len(instance.BidirectionalArc)))
	print("")
	print("GeneratorTypes: "+str(len(instance.Generator)))
	print("TotalGenerators: "+str(len(instance.GeneratorsOfNode)))
	print("StorageTypes: "+str(len(instance.Storage)))
	print("TotalStorages: "+str(len(instance.StoragesOfNode)))
	print("")
	print("InvestmentYears: "+str(len(instance.Period)))
	print("Scenarios: "+str(len(instance.Scenario)))
	print("TotalOperationalHoursPerScenario: "+str(len(instance.Operationalhour)))
	print("TotalOperationalHoursPerInvYear: "+str(len(instance.Operationalhour)*len(instance.Scenario)))
	print("Seasons: "+str(len(instance.Season)))
	print("RegularSeasons: "+str(len(instance.FirstHoursOfRegSeason)))
	print("LengthRegSeason: "+str(value(instance.lengthRegSeason)))
	print("PeakSeasons: "+str(len(instance.FirstHoursOfPeakSeason)))
	print("LengthPeakSeason: "+str(value(instance.lengthPeakSeason)))
	print("")
	print("Discount rate: "+str(value(instance.discountrate)))
	print(f"Operational discount scale: {value(instance.operationalDiscountrate):.3f}")
	print("Optimizing with hydrogen: " + str(hydrogen))
	print("--------------------------------------------------------------\n")

	if WRITE_LP:
		print("Writing LP-file...")
		start = time.time()
		lpstring = 'LP_' + name + '.lp'
		if USE_TEMP_DIR:
			lpstring = './LP_'+ name + '.lp'
		instance.write(lpstring, io_options={'symbolic_solver_labels': True})
		end = time.time()
		print("Writing LP-file took:")
		print(str(end - start))

	print("{hour}:{minute}:{second}: Solving...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))

	if solver == "CPLEX":
		opt = SolverFactory("cplex", Verbose=True)
		opt.options["lpmethod"] = 4
		opt.options["barrier crossover"] = -1
		if TIME_LIMIT is not None and TIME_LIMIT > 0:
			opt.options['timelimit'] = TIME_LIMIT
		#instance.display('outputs_cplex.txt')
	if solver == "Xpress":
		opt = SolverFactory("xpress") #Verbose=True
		opt.options["defaultAlg"] = 4
		opt.options["crossover"] = 0
		opt.options["lpLog"] = 1
		opt.options["Trace"] = 1
		if TIME_LIMIT is not None and TIME_LIMIT > 0:
			opt.options['maxtime'] = TIME_LIMIT
		#instance.display('outputs_xpress.txt')
	if solver == "Gurobi":
		opt = SolverFactory('gurobi', Verbose=True)
		opt.options["Crossover"]=-1
		# opt.options["CrossoverBasis"]=0

		# opt.options["MIRCuts"]=2
		# opt.options["FlowCoverCuts"]=1
		# opt.options["Cuts"]=0


		# opt.options["FlowPathCuts"]=2
		# opt.options["Presolve"]=2
		opt.options['NumericFocus']=3
		opt.options['PreSparsify']=1

		opt.options['BarConvTol']=1e-20
		opt.options['MIPGap']=1e-2
		# opt.options['MIPGap']=2e-1
		# opt.options['Threads']=32
		# opt.options['slog'] = 1
		# opt.options['FeasibilityTol']=10**(-9)
		if TIME_LIMIT is not None and TIME_LIMIT > 0:
			opt.options['TimeLimit'] = TIME_LIMIT
		opt.options["Method"]=2
		opt.options["NodeMethod"]=2

	# try:
	results = opt.solve(instance, tee=True, logfile=result_file_path + '/logfile_' + name + '.log')#, keepfiles=True, symbolic_solver_labels=True)
	# except:
	# 	print('{hour}:{minute}:{second}: ERROR: Could not load results. Likely cause: time limit reached'.format(
	# 		hour=datetime.now().strftime("%H"), minute = datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	# 	exit()

	if PICKLE_INSTANCE:
		start = time.time()
		picklestring = 'instance' + name + '.pkl'
		if USE_TEMP_DIR:
			picklestring = temp_dir + '/instance' + name + '.pkl'
		with open(picklestring, mode='wb') as file:
			cloudpickle.dump(instance, file)
		end = time.time()
		print("Pickling instance took [sec]:")
		print(str(end - start))

	endOptimization = StartReporting = datetime.now()

	#instance.display('outputs_gurobi.txt')

	#import pdb; pdb.set_trace()

	###########
	##RESULTS##
	###########

	def calculatePowerEmissionIntensity(n,h,i,w,m=None):
			#print(f'Evaluating {n}')
			emissions = 1000 * value(sum(instance.genOperational[n,g,h,i,w]*instance.genCO2TypeFactor[g]*(3.6/instance.genEfficiency[g,i]) for g in instance.Generator if (n,g) in instance.GeneratorsOfNode))
			total_power = value(sum(instance.genOperational[n,g,h,i,w] for g in instance.Generator if (n,g) in instance.GeneratorsOfNode))
			for n2 in instance.NodesLinked[n]:
				if value(instance.lineEfficiency[n2,n]*instance.transmissionOperational[n2,n,h,i,w]) > 200 and value(instance.lineEfficiency[n,n2]*instance.transmissionOperational[n,n2,h,i,w]) < 1:
					if n2==m:
						print(f'Warning: We have recursion loop between {n} and {m} in calculatePowerEmissionIntensity!')
					emissions += calculatePowerEmissionIntensity(n2,h,i,w,n) * value(instance.lineEfficiency[n2,n]*instance.transmissionOperational[n2,n,h,i,w])
					total_power += value(instance.lineEfficiency[n2,n]*instance.transmissionOperational[n2,n,h,i,w])
				else:
					emissions += 0
					total_power += 0
			if total_power > 0:
				emission_factor = emissions/total_power
			else:
				emission_factor = 0
				# print(f'Warning: Total power in {n} in hour {h} in {inv_per[int(i-1)]} in {w} is 0!')
			# print(f'Node {n}, period {i}, hour {h}, {w}:\tEm.fac.:{emission_factor:.3f} kg/MWh')
			return emission_factor

	print(("{hour}:{minute}:{second}: Writing results in " + result_file_path + '/\n').format(
		hour=datetime.now().strftime("%H"), minute = datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))

	f = open(result_file_path + "/" + 'results_objective.csv', 'w', newline='')
	writer = csv.writer(f)
	writer.writerow(["Objective function value:" + str(value(instance.Obj))])
	writer.writerow(["Scientific notation:", str(value(instance.Obj))])
	writer.writerow(["Solver status:",results.solver.status])
	f.close()

	print("{hour}:{minute}:{second}: Writing transmission investment decisions to results_output_transmission.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_output_transmission.csv', 'w', newline='')
	writer = csv.writer(f)
	writer.writerow(["BetweenNode","AndNode","Period","transmissionBuilt?","transmissionInvCap_MW","transmissionInvCapMax_MW","transmissionInstalledCap_MW","transmissionInstalledCapMax_MW","DiscountedInvestmentCost_Euro","transmissionExpectedAnnualVolume_GWh","ExpectedAnnualLosses_GWh"])
	for (n1,n2) in instance.BidirectionalArc:
		for i in instance.Period:
			writer.writerow([n1,n2,inv_per[int(i-1)],
							 value(instance.transmissionFixedCostPaid[n1,n2,i]),
							 value(instance.transmissionInvCap[n1,n2,i]), value(instance.transmissionMaxBuiltCap[n1,n2,i]),
							 value(instance.transmissionInstalledCap[n1,n2,i]), value(instance.transmissionMaxInstalledCap[n1,n2,i]),
							 value(instance.discount_multiplier[i]*instance.transmissionInvCap[n1,n2,i]*instance.transmissionVarInvCost[n1,n2,i]),
							 value(sum(instance.sceProbab[w]*instance.seasScale[s]*(instance.transmissionOperational[n1,n2,h,i,w]+instance.transmissionOperational[n2,n1,h,i,w])/1000 for (s,h) in instance.HoursOfSeason for w in instance.Scenario)),
							 value(sum(instance.sceProbab[w]*instance.seasScale[s]*((1 - instance.lineEfficiency[n1,n2])*instance.transmissionOperational[n1,n2,h,i,w] + (1 - instance.lineEfficiency[n2,n1])*instance.transmissionOperational[n2,n1,h,i,w])/1000 for (s,h) in instance.HoursOfSeason for w in instance.Scenario))])
	f.close()

	print("{hour}:{minute}:{second}: Writing generator investment decisions to results_output_gen.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))

	f = open(result_file_path + "/" + 'results_output_gen.csv', 'w', newline='')
	writer = csv.writer(f)
	my_string = ["Node","GeneratorType","Period","genInvCap_MW","genInstalledCap_MW","genExpectedCapacityFactor","DiscountedInvestmentCost_Euro","genExpectedAnnualProduction_GWh"]
	writer.writerow(my_string)
	for (n,g) in instance.GeneratorsOfNode:
		for i in instance.Period:
			my_string=[n,g,inv_per[int(i-1)],value(instance.genInvCap[n,g,i]),value(instance.genInstalledCap[n,g,i]),
					   value(sum(instance.sceProbab[w]*instance.seasScale[s]*instance.genOperational[n,g,h,i,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario)/(instance.genInstalledCap[n,g,i]*8760) if value(instance.genInstalledCap[n,g,i]) != 0 and value(instance.genInstalledCap[n,g,i]) > 3 else 0),
					   value(instance.discount_multiplier[i]*instance.genInvCap[n,g,i]*instance.genInvCost[g,i]),
					   value(sum(instance.seasScale[s]*instance.sceProbab[w]*instance.genOperational[n,g,h,i,w]/1000 for (s,h) in instance.HoursOfSeason for w in instance.Scenario) if value(instance.genInstalledCap[n,g,i]) > 3 else 0)]
			writer.writerow(my_string)
	f.close()

	print("{hour}:{minute}:{second}: Writing storage investment decisions to results_output_stor.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_output_stor.csv', 'w', newline='')
	writer = csv.writer(f)
	writer.writerow(["Node","StorageType","Period","storPWInvCap_MW","storPWInstalledCap_MW","storENInvCap_MWh","storENInstalledCap_MWh","DiscountedInvestmentCostPWEN_EuroPerMWMWh","ExpectedAnnualDischargeVolume_GWh","ExpectedAnnualLossesChargeDischarge_GWh"])
	for (n,b) in instance.StoragesOfNode:
		for i in instance.Period:
			writer.writerow([n,b,inv_per[int(i-1)],value(instance.storPWInvCap[n,b,i]),value(instance.storPWInstalledCap[n,b,i]),
							 value(instance.storENInvCap[n,b,i]),value(instance.storENInstalledCap[n,b,i]),
							 value(instance.discount_multiplier[i]*(instance.storPWInvCap[n,b,i]*instance.storPWInvCost[b,i] + instance.storENInvCap[n,b,i]*instance.storENInvCost[b,i])),
							 value(sum(instance.sceProbab[w]*instance.seasScale[s]*instance.storDischarge[n,b,h,i,w]/1000 for (s,h) in instance.HoursOfSeason for w in instance.Scenario)),
							 value(sum(instance.sceProbab[w]*instance.seasScale[s]*((1 - instance.storageDischargeEff[b])*instance.storDischarge[n,b,h,i,w] + (1 - instance.storageChargeEff[b])*instance.storCharge[n,b,h,i,w])/1000 for (s,h) in instance.HoursOfSeason for w in instance.Scenario))])
	f.close()

	print("{hour}:{minute}:{second}: Writing transmission operational decisions to results_output_transmission_operational.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_output_transmission_operational.csv', 'w', newline='')
	writer = csv.writer(f)
	writer.writerow(["FromNode","ToNode","Period","Season","Scenario","Hour","TransmissionReceived_MW","Losses_MW"])
	for (n1,n2) in instance.DirectionalLink:
		for i in instance.Period:
			for (s,h) in instance.HoursOfSeason:
				for w in instance.Scenario:
					transmissionSent = value(instance.transmissionOperational[n1,n2,h,i,w])
					writer.writerow([n1,n2,inv_per[int(i-1)],s,w,h,
									 value(instance.lineEfficiency[n1,n2])*transmissionSent,
									 value((1 - instance.lineEfficiency[n1,n2]))*transmissionSent])
	f.close()

	print(
		"{hour}:{minute}:{second}: Writing power balances to results_power_balance.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"),
			second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_power_balance.csv', 'w', newline='')
	writer = csv.writer(f)
	header = ["Node", "Period", "Season", "Hour", "Scenario", "Available power [MWh]", "Power generation [MWh]", "Power curtailed [MWh]", "Power transmission in [MWh]","Power storage discharge [MWh]", "Power transmission out [MWh]", "Power storage charge [MWh]", "Power load [MWh]", "Power shed [MWh]"]
	if hydrogen is True:
		header.append("Power for hydrogen [MWh]")
	writer.writerow(header)
	for n in instance.Node:
		for i in instance.Period:
			for (s,h) in instance.HoursOfSeason:
				for w in instance.Scenario:
					row = [n,inv_per[int(i-1)],s,h,w]
					row.append(value(sum(instance.genCapAvail[n,g,h,w,i]*instance.genInstalledCap[n,g,i] for g in instance.Generator if (n,g) in instance.GeneratorsOfNode)))
					row.append(value(sum(instance.genOperational[n,g,h,i,w] for g in instance.Generator if (n,g) in instance.GeneratorsOfNode)))
					row.append(value(sum((instance.genCapAvail[n,g,h,w,i]*instance.genInstalledCap[n,g,i] - instance.genOperational[n,g,h,i,w]) for g in instance.Generator if (n,g) in instance.GeneratorsOfNode)))
					row.append(value(sum(instance.lineEfficiency[link,n]*instance.transmissionOperational[link,n,h,i,w] for link in instance.NodesLinked[n])))
					row.append(value(sum(instance.storageDischargeEff[b] * instance.storDischarge[n, b, h, i, w] for b in instance.Storage if (n, b) in instance.StoragesOfNode)))
					row.append(value(sum(instance.transmissionOperational[n,link,h,i,w] for link in instance.NodesLinked[n])))
					row.append(value(sum(instance.storCharge[n,b,h,i,w] for b in instance.Storage if (n,b) in instance.StoragesOfNode)))
					row.append(value(instance.sload[n,h,i,w]))
					row.append(value(instance.loadShed[n,h,i,w]))
					if hydrogen is True and n in instance.HydrogenProdNode:
						row.append(value(sum(instance.powerForHydrogen[n,j,h,i,w] for j in instance.Period if j<=i)))
					writer.writerow(row)
	f.close()

	print("{hour}:{minute}:{second}: Writing curtailed power to results_output_curtailed_prod.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_output_curtailed_prod.csv', 'w', newline='')
	writer = csv.writer(f)
	writer.writerow(["Node","RESGeneratorType","Period","ExpectedAnnualCurtailment_GWh", "Expected total available power_GWh", "Expected annual curtailment ratio of total capacity_%"])
	for t in instance.Technology:
		if t == 'Hydro_ror' or t == 'Wind_onshr' or t == 'Wind_offshr_grounded' or t == 'Wind_offshr_floating' or t == 'Solar':
			for (n,g) in instance.GeneratorsOfNode:
				if (t,g) in instance.GeneratorsOfTechnology:
					for i in instance.Period:
						curtailedPower = value(sum(instance.sceProbab[w]*instance.seasScale[s]*(instance.genCapAvail[n,g,h,w,i]*instance.genInstalledCap[n,g,i] - instance.genOperational[n,g,h,i,w])/1000 for w in instance.Scenario for (s,h) in instance.HoursOfSeason))
						totalPowerProduction = value(sum(instance.sceProbab[w]*instance.seasScale[s]*(instance.genCapAvail[n,g,h,w,i]*instance.genInstalledCap[n,g,i])/1000 for w in instance.Scenario for (s,h) in instance.HoursOfSeason))
						row = [n,g,inv_per[int(i-1)], curtailedPower, totalPowerProduction]
						if totalPowerProduction > 0:
							row.append(curtailedPower/totalPowerProduction*100)
						else:
							row.append(0)
						writer.writerow(row)
	f.close()

	#Commenting out this plotting because it is not currently interesting.

	# print("{hour}:{minute}:{second}: Writing plotting file to results_output_EuropePlot.csv...".format(
	#     hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	# f = open(result_file_path + "/" + 'results_output_EuropePlot.csv', 'w', newline='')
	# writer = csv.writer(f)
	# writer.writerow(["Period","genInstalledCap_MW"])
	# my_string=[""]
	# for g in instance.Generator:
	#     my_string.append(g)
	# writer.writerow(my_string)
	# my_string=["Initial"]
	# for g in instance.Generator:
	#     my_string.append((value(sum(instance.genInitCap[n,g,1] for n in instance.Node if (n,g) in instance.GeneratorsOfNode))))
	# writer.writerow(my_string)
	# for i in instance.Period:
	#     my_string=[inv_per[int(i-1)]]
	#     for g in instance.Generator:
	#         my_string.append(value(sum(instance.genInstalledCap[n,g,i] for n in instance.Node if (n,g) in instance.GeneratorsOfNode)))
	#     writer.writerow(my_string)
	# writer.writerow([""])
	# writer.writerow(["Period","genExpectedAnnualProduction_GWh"])
	# my_string=[""]
	# for g in instance.Generator:
	#     my_string.append(g)
	# writer.writerow(my_string)
	# for i in instance.Period:
	#     my_string=[inv_per[int(i-1)]]
	#     for g in instance.Generator:
	#         my_string.append(value(sum(instance.sceProbab[w]*instance.seasScale[s]*instance.genOperational[n,g,h,i,w]/1000 for n in instance.Node if (n,g) in instance.GeneratorsOfNode for (s,h) in instance.HoursOfSeason for w in instance.Scenario)))
	#     writer.writerow(my_string)
	# writer.writerow([""])
	# writer.writerow(["Period","storPWInstalledCap_MW"])
	# my_string=[""]
	# for b in instance.Storage:
	#     my_string.append(b)
	# writer.writerow(my_string)
	# for i in instance.Period:
	#     my_string=[inv_per[int(i-1)]]
	#     for b in instance.Storage:
	#         my_string.append(value(sum(instance.storPWInstalledCap[n,b,i] for n in instance.Node if (n,b) in instance.StoragesOfNode)))
	#     writer.writerow(my_string)
	# writer.writerow([""])
	# writer.writerow(["Period","storENInstalledCap_MW"])
	# my_string=[""]
	# for b in instance.Storage:
	#     my_string.append(b)
	# writer.writerow(my_string)
	# for i in instance.Period:
	#     my_string=[inv_per[int(i-1)]]
	#     for b in instance.Storage:
	#         my_string.append(value(sum(instance.storENInstalledCap[n,b,i] for n in instance.Node if (n,b) in instance.StoragesOfNode)))
	#     writer.writerow(my_string)
	# writer.writerow([""])
	# writer.writerow(["Period","storExpectedAnnualDischarge_GWh"])
	# my_string=[""]
	# for b in instance.Storage:
	#     my_string.append(b)
	# writer.writerow(my_string)
	# for i in instance.Period:
	#     my_string=[inv_per[int(i-1)]]
	#     for b in instance.Storage:
	#         my_string.append(value(sum(instance.sceProbab[w]*instance.seasScale[s]*instance.storDischarge[n,b,h,i,w]/1000 for n in instance.Node if (n,b) in instance.StoragesOfNode for (s,h) in instance.HoursOfSeason for w in instance.Scenario)))
	#     writer.writerow(my_string)
	# f.close()

	print("{hour}:{minute}:{second}: Writing summary file to results_output_EuropeSummary.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_output_EuropeSummary.csv', 'w', newline='')
	if solver == 'Xpress':
		fError = open(result_file_path + "/" + "errorLog.log",'w')
	writer = csv.writer(f)
	header = ["Period","Scenario","AnnualCO2emission_Ton","CO2Price_EuroPerTon","CO2Cap_Ton","AnnualGeneration_GWh","AvgCO2factor_TonPerMWh","AvgPowerPrice_Euro","TotAnnualCurtailedRES_GWh","TotAnnualLossesChargeDischarge_GWh","AnnualLossesTransmission_GWh"]
	if hydrogen is True:
		header.append("AvgH2MarginalCost_EuroPerKg")
		# header.append('AverageEmissionsH2_kgCO2PerKgH2')
	writer.writerow(header)
	for i in instance.Period:
		for w in instance.Scenario:
			power_co2_factor = value(sum(instance.seasScale[s]*instance.genOperational[n,g,h,i,w]*instance.genCO2TypeFactor[g]*(3.6/instance.genEfficiency[g,i]) for (n,g) in instance.GeneratorsOfNode for (s,h) in instance.HoursOfSeason)/sum(instance.seasScale[s]*instance.genOperational[n,g,h,i,w] for (n,g) in instance.GeneratorsOfNode for (s,h) in instance.HoursOfSeason))
			my_string=[inv_per[int(i-1)],w,
					   value(sum(instance.seasScale[s]*instance.genOperational[n,g,h,i,w]*instance.genCO2TypeFactor[g]*(3.6/instance.genEfficiency[g,i]) for (n,g) in instance.GeneratorsOfNode for (s,h) in instance.HoursOfSeason))]
			if EMISSION_CAP:
				try:
					my_string.append(-value(instance.dual[instance.emission_cap[i,w]]/(instance.discount_multiplier[i]*instance.operationalDiscountrate*instance.sceProbab[w]*1e6)))
					my_string.append(value(instance.CO2cap[i]*1e6))
				except:
					my_string.append(0)
					my_string.append(value(instance.CO2cap[i]*1e6))
			else:
				my_string.append(value(instance.CO2price[i]))
				my_string.append("INF")
			my_string.extend([value(sum(instance.seasScale[s]*instance.genOperational[n,g,h,i,w]/1000 for (n,g) in instance.GeneratorsOfNode for (s,h) in instance.HoursOfSeason)),
							  power_co2_factor,
							  # value(sum(instance.dual[instance.FlowBalance[n,h,i,w]]/(instance.discount_multiplier[i]*instance.operationalDiscountrate*instance.seasScale[s]*instance.sceProbab[w]) for n in instance.Node for (s,h) in instance.HoursOfSeason)/value(len(instance.HoursOfSeason)*len(instance.Node))),
							  0, # No duals with MIPs
							  value(sum(instance.seasScale[s]*(instance.genCapAvail[n,g,h,w,i]*instance.genInstalledCap[n,g,i] - instance.genOperational[n,g,h,i,w])/1000 for (n,g) in instance.GeneratorsOfNode if g == 'Hydrorun-of-the-river' or g == 'Windonshore' or g == 'Windoffshore' or g == 'Solar' for (s,h) in instance.HoursOfSeason)),
							  value(sum(instance.seasScale[s]*((1 - instance.storageDischargeEff[b])*instance.storDischarge[n,b,h,i,w] + (1 - instance.storageChargeEff[b])*instance.storCharge[n,b,h,i,w])/1000 for (n,b) in instance.StoragesOfNode for (s,h) in instance.HoursOfSeason)),
							  value(sum(instance.seasScale[s]*((1 - instance.lineEfficiency[n1,n2])*instance.transmissionOperational[n1,n2,h,i,w] + (1 - instance.lineEfficiency[n2,n1])*instance.transmissionOperational[n2,n1,h,i,w])/1000 for (n1,n2) in instance.BidirectionalArc for (s,h) in instance.HoursOfSeason))])
			if hydrogen is True:
				try:
					my_string.extend([value(sum(instance.dual[instance.hydrogen_flow_balance[n,h,i,w]]/(instance.discount_multiplier[i]*instance.operationalDiscountrate*instance.seasScale[s]*instance.sceProbab[w]) for n in instance.HydrogenProdNode for (s,h) in instance.HoursOfSeason)/value(len(HoursOfSeason)*len(instance.HydrogenProdNode)))])
				except:
					# print('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_output_EuropeSummary.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+')')
					# fError.write('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_output_EuropeSummary.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+'). Dual is likely 0 (poorly handled in Xpress)\n')
					my_string.extend([0])
				# green_h2_production = value(sum(instance.seasScale[s]*instance.powerForHydrogen[n,j,h,i,w] for j in instance.Period if j<=i for n in instance.HydrogenProdNode for (s,h) in instance.HoursOfSeason))
				# green_h2_emissions = 1000 * green_h2_production * power_co2_factor
				# blue_h2_production = value(sum(instance.seasScale[s] * instance.hydrogenProducedReformer_kg[n,p,h,i,w] for p in instance.ReformerPlants for n in instance.ReformerLocations for (s,h) in instance.HoursOfSeason))
				# blue_h2_direct_emissions = value(sum(instance.seasScale[s] * instance.ReformerEmissionFactor[p,i] * instance.hydrogenProducedReformer_kg[n,p,h,i,w] for n in instance.ReformerLocations for p in instance.ReformerPlants for (s,h) in instance.HoursOfSeason))
				# blue_h2_power_emissions = 1000 * power_co2_factor * value(sum(instance.seasScale[s] * instance.ReformerPlantElectricityUse[p,i] * instance.hydrogenProducedReformer_kg[n,p,h,i,w] for n in instance.ReformerLocations for p in instance.ReformerPlants for (s,h) in instance.HoursOfSeason))
				# my_string.extend([(green_h2_emissions + blue_h2_direct_emissions + blue_h2_power_emissions)/(green_h2_production + blue_h2_production)])
			writer.writerow(my_string)
	if solver == 'Xpress':
		fError.write('\n')
		fError.close()
	writer.writerow([""])
	writer.writerow(["GeneratorType","Period","genInvCap_MW","genInstalledCap_MW","TotDiscountedInvestmentCost_Euro","genExpectedAnnualProduction_GWh"])
	for g in instance.Generator:
		for i in instance.Period:
			expected_production = 0
			for n in instance.Node:
				if (n,g) in instance.GeneratorsOfNode:
					expected_production += value(sum(instance.seasScale[s] * instance.sceProbab[w] * instance.genOperational[n,g,h,i,w] / 1000 for (s,h) in instance.HoursOfSeason for w in instance.Scenario))
			writer.writerow([g,inv_per[int(i-1)],value(sum(instance.genInvCap[n,g,i] for n in instance.Node if (n,g) in instance.GeneratorsOfNode)),
							 value(sum(instance.genInstalledCap[n,g,i] for n in instance.Node if (n,g) in instance.GeneratorsOfNode)),
							 value(sum(instance.discount_multiplier[i]*instance.genInvCap[n,g,i]*instance.genInvCost[g,i] for n in instance.Node if (n,g) in instance.GeneratorsOfNode)),
							 expected_production])
	writer.writerow([""])
	writer.writerow(["StorageType","Period","storPWInvCap_MW","storPWInstalledCap_MW","storENInvCap_MWh","storENInstalledCap_MWh","TotDiscountedInvestmentCostPWEN_Euro","ExpectedAnnualDischargeVolume_GWh"])
	for b in instance.Storage:
		for i in instance.Period:
			expected_discharge = 0
			for n in instance.Node:
				if (n,b) in instance.StoragesOfNode:
					expected_discharge += value(sum(instance.seasScale[s] * instance.sceProbab[w] * instance.storDischarge[n,b,h,i,w] / 1000 for (s,h) in instance.HoursOfSeason or w in instance.Scenario))
			writer.writerow([b,inv_per[int(i-1)],value(sum(instance.storPWInvCap[n,b,i] for n in instance.Node if (n,b) in instance.StoragesOfNode)),
							 value(sum(instance.storPWInstalledCap[n,b,i] for n in instance.Node if (n,b) in instance.StoragesOfNode)),
							 value(sum(instance.storENInvCap[n,b,i] for n in instance.Node if (n,b) in instance.StoragesOfNode)),
							 value(sum(instance.storENInstalledCap[n,b,i] for n in instance.Node if (n,b) in instance.StoragesOfNode)),
							 value(sum(instance.discount_multiplier[i]*(instance.storPWInvCap[n,b,i]*instance.storPWInvCost[b,i] + instance.storENInvCap[n,b,i]*instance.storENInvCost[b,i]) for n in instance.Node if (n,b) in instance.StoragesOfNode)),
							 expected_discharge])
	f.close()

	print("{hour}:{minute}:{second}: Writing operational results to results_output_Operational.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
	f = open(result_file_path + "/" + 'results_output_Operational.csv', 'w', newline='')
	writer = csv.writer(f)
	my_header = ["Node","Period","Scenario","Season","Hour"]
	for g in instance.Generator:
		my_string = str(g)+"_MW"
		my_header.append(my_string)
	my_header.extend(["storCharge_MW","storDischarge_MW","storEnergyLevel_MWh","LossesChargeDischargeBleed_MW","FlowOut_MW","FlowIn_MW","LossesFlowIn_MW","LoadShed_MW","Price_EURperMWh","AvgCO2_kgCO2perMWh_PRODUCTION","AvgCO2_kgCO2perMWh_TOTAL"])
	writer.writerow(my_header)
	for n in instance.Node:
			for i in instance.Period:
				for w in instance.Scenario:
					for (s,h) in instance.HoursOfSeason:
						my_string=[n,inv_per[int(i-1)],w,s,h]
						for g in instance.Generator:
							if (n,g) in instance.GeneratorsOfNode:
								my_string.append(value(instance.genOperational[n,g,h,i,w]))
							else:
								my_string.append(0)
						my_string.extend([value(sum(-instance.storCharge[n,b,h,i,w] for b in instance.Storage if (n,b) in instance.StoragesOfNode)),
										  value(sum(instance.storDischarge[n,b,h,i,w] for b in instance.Storage if (n,b) in instance.StoragesOfNode)),
										  value(sum(instance.storOperational[n,b,h,i,w] for b in instance.Storage if (n,b) in instance.StoragesOfNode)),
										  value(sum(-(1 - instance.storageDischargeEff[b])*instance.storDischarge[n,b,h,i,w] - (1 - instance.storageChargeEff[b])*instance.storCharge[n,b,h,i,w] - (1 - instance.storageBleedEff[b])*instance.storOperational[n,b,h,i,w] for b in instance.Storage if (n,b) in instance.StoragesOfNode)),
										  value(sum(-instance.transmissionOperational[n,link,h,i,w] for link in instance.NodesLinked[n])),
										  value(sum(instance.transmissionOperational[link,n,h,i,w] for link in instance.NodesLinked[n])),
										  value(sum(-(1 - instance.lineEfficiency[link,n])*instance.transmissionOperational[link,n,h,i,w] for link in instance.NodesLinked[n])),
										  value(instance.loadShed[n,h,i,w]),
										  # value(instance.dual[instance.FlowBalance[n,h,i,w]]/(instance.discount_multiplier[i]*instance.operationalDiscountrate*instance.seasScale[s]*instance.sceProbab[w]))])
										  0]) # No duals with MIPs
						if value(sum(instance.genOperational[n,g,h,i,w] for g in instance.Generator if (n,g) in instance.GeneratorsOfNode)) > 0:
							my_string.extend([value(1000*sum(instance.genOperational[n,g,h,i,w]*instance.genCO2TypeFactor[g]*(3.6/instance.genEfficiency[g,i]) for g in instance.Generator if (n,g) in instance.GeneratorsOfNode)/sum(instance.genOperational[n,g,h,i,w] for g in instance.Generator if (n,g) in instance.GeneratorsOfNode))])
						my_string.extend([calculatePowerEmissionIntensity(n,h,i,w)])
						if verboseResultWriting is True:
							print("{hour}:{minute}:{second}: ".format(hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"),
																	  second=datetime.now().strftime("%S")) + str(my_string))
						writer.writerow(my_string)
	f.close()

	if hydrogen is True:
		print("{hour}:{minute}:{second}: Writing hydrogen investment results to results_hydrogen_production_investments.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + "/" + 'results_hydrogen_production_investments.csv', 'w', newline='')
		if solver == 'Xpress':
			fError = open(result_file_path + "/" + "errorLog.log",'a')
		writer = csv.writer(f)
		my_header = ["Node","Period","New electrolyzer capacity [MW]", "Total electrolyzer capacity [MW]", "New electrolyzer capacity [ton/h]", "Total electrolyzer capacity [ton/h]",
					 "Expected annual power usage [MWh]","Expected annual electrolyzer hydrogen production [ton]",
					 'Expected electrolyzer capacity factor', 'New Reformer capacity [ton/h]', 'Total Reformer capacity [ton/h]',
					 'Expected annual reformer hydrogen production [ton]',"Expected marginal price of hydrogen [EUR/ton]"]
		writer.writerow(my_header)
		for n in instance.HydrogenProdNode:
			for i in instance.Period:
				try:
					# hydrogenPrice = value(sum(instance.sceProbab[w]*instance.dual[instance.hydrogen_flow_balance[n,h,i,w]] for w in instance.Scenario))
					hydrogenPrice = 0
				except:
					fError.write('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_hydrogen_investments.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+'). Dual is likely 0 (poorly handled in Xpress). Setting dual = 0.\n')
					hydrogenPrice = 0
				if n in instance.ReformerLocations:
					ReformerCapBuilt = value(sum(instance.ReformerCapBuilt[n,p,i] for p in instance.ReformerPlants)/(instance.hydrogenLHV_kg * 1000))
					reformerCapTotal = value(sum(instance.ReformerTotalCap[n,p,i] for p in instance.ReformerPlants)/(instance.hydrogenLHV_kg * 1000))
					reformerExpectedProduction = value(sum(instance.seasScale[s] * instance.sceProbab[w] * instance.hydrogenProducedReformer_ton[n,p,h,i,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario for p in instance.ReformerPlants))
				else:
					ReformerCapBuilt = 0
					reformerCapTotal = 0
					reformerExpectedProduction = 0
				electrolyzerCapacity = value(sum(instance.elyzerTotalCap[n,j,i] / instance.elyzerPowerConsumptionPerTon[j] for j in instance.Period if j<=i))
				expectedElectrolyzerProduction = value(sum(instance.sceProbab[w] * instance.seasScale[s] * instance.hydrogenProducedElectro[n,h,i,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario))
				electrolyzerCapFactor = (expectedElectrolyzerProduction/(electrolyzerCapacity*8760) if electrolyzerCapacity > 10 else 0)
				writer.writerow([n,inv_per[int(i-1)],
								 value(instance.elyzerCapBuilt[n,i,i]),
								 value(sum(instance.elyzerTotalCap[n,j,i] for j in instance.Period if j <= i)),
								 value(sum(instance.elyzerCapBuilt[n,j,i] / instance.elyzerPowerConsumptionPerTon[j] for j in instance.Period if j<=i)),
								 electrolyzerCapacity,
								 value(sum(instance.seasScale[s] * instance.sceProbab[w] * instance.powerForHydrogen[n,j,h,i,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario for j in instance.Period if j<=i)),
								 expectedElectrolyzerProduction,
								 electrolyzerCapFactor,
								 ReformerCapBuilt,
								 reformerCapTotal,
								 reformerExpectedProduction,
								 hydrogenPrice])
		f.close()
		if solver == 'Xpress':
			fError.write('\n')
			fError.close()
		f.close()


		# print("{hour}:{minute}:{second}: Writing detailed Reformer investment results to results_hydrogen_electrolyzer_detailed_check.csv...".format(
		# 	hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + "/" + 'results_hydrogen_electrolyzer_detailed_check.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ['Node','Buying Period','Operating period','New capacity','Total capacity']
		# writer.writerow(my_header)
		# for n in instance.HydrogenProdNode:
		# 	for j in instance.Period:
		# 		for i in instance.Period:
		# 			my_string= [n,j,i,
		# 						value(instance.elyzerCapBuilt[n,j,i]),
		# 						value(instance.elyzerTotalCap[n,j,i])]
		# 			writer.writerow(my_string)
		# f.close()


		print("{hour}:{minute}:{second}: Writing detailed reformer investment results to results_hydrogen_reformer_detailed_investments.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + "/" + 'results_hydrogen_reformer_detailed_investments.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ['Node','Reformer plant type','Period','New capacity [MW]','Total capacity [MW]','New capacity [ton/h]','Total capacity [ton/h]',
					 'Expected production [ton H2/year]', 'Expected capacity factor [%]', 'Expected emissions [tons CO2/year]', 'Expected electricity consumption [GWh]']
		writer.writerow(my_header)
		for n in instance.ReformerLocations:
			for p in instance.ReformerPlants:
				for i in instance.Period:
					reformerCap = value(instance.ReformerTotalCap[n,p,i])
					reformerProduction = value(sum(instance.sceProbab[w] * instance.seasScale[s] * instance.hydrogenProducedReformer_ton[n,p,h,i,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario))
					capFactor = (reformerProduction / (8760*(reformerCap/value(instance.hydrogenLHV_kg * 1000))) if reformerCap > 1 else 0)
					my_string = [n,p,inv_per[int(i)-1],
								 value(instance.ReformerCapBuilt[n,p,i]),
								 reformerCap,
								 value(instance.ReformerCapBuilt[n,p,i]/(instance.hydrogenLHV_kg * 1000)),
								 reformerCap/(value(instance.hydrogenLHV_kg) * 1000),
								 reformerProduction,
								 capFactor,
								 reformerProduction * value(instance.ReformerEmissionFactor[p,i]),
								 reformerProduction * value(instance.ReformerPlantElectricityUse[p,i]/1000)]
					writer.writerow(my_string)
		f.close()


		if h2storage is True:
			print("{hour}:{minute}:{second}: Writing hydrogen storage investment results to results_hydrogen_storage_investments.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
			f = open(result_file_path + "/" + 'results_hydrogen_storage_investments.csv', 'w', newline='')
			writer = csv.writer(f)
			my_header = ['Node','Period','New storage capacity [ton]','Total storage capacity [ton]', 'Discounted cost of new capacity [EUR]','Discounted total cost [EUR]']
			writer.writerow(my_header)
			for n in instance.HydrogenProdNode:
				for i in instance.Period:
					my_string = [n,inv_per[int(i)-1],
								 value(instance.hydrogenStorageBuilt[n,i]),
								 value(instance.hydrogenTotalStorage[n,i]),
								 value(instance.hydrogenStorageBuilt[n,i] * instance.hydrogenStorageInvCost[i]),
								 value(sum(instance.hydrogenStorageBuilt[n,j] * instance.hydrogenStorageInvCost[j] for j in instance.Period if j<=i))]
					writer.writerow(my_string)
			f.close()

			print("{hour}:{minute}:{second}: Writing hydrogen storage operational results to results_hydrogen_storage_operational.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
			f = open(result_file_path + "/" + 'results_hydrogen_storage_operational.csv', 'w', newline='')
			writer = csv.writer(f)
			my_header = ['Node','Period','Scenario', 'Season',' Hour','Initial storage [ton]','Charge [ton]','Discharge [ton]','Final stored [ton]']
			writer.writerow(my_header)
			for n in instance.HydrogenProdNode:
				for i in instance.Period:
					for w in instance.Scenario:
						for (s,h) in instance.HoursOfSeason:
							my_string= [n,inv_per[i-1], w, s, h]
							if h in instance.FirstHoursOfRegSeason or h in instance.FirstHoursOfPeakSeason:
								my_string.extend([value(instance.hydrogenStorageInitOperational * instance.hydrogenTotalStorage[n,i])])
							else:
								my_string.extend([value(instance.hydrogenStorageOperational[n,h-1,i,w])])
							my_string.extend([value(instance.hydrogenChargeStorage[n,h,i,w]),
											  value(instance.hydrogenDischargeStorage[n,h,i,w]),
											  value(instance.hydrogenStorageOperational[n,h,i,w])])
							writer.writerow(my_string)
			f.close()


		print("{hour}:{minute}:{second}: Writing hydrogen production results to results_hydrogen_production.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_hydrogen_production.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["Node", "Period", "Scenario", "Season", "Hour", "Power for hydrogen [MWh]", "Electrolyzer production [ton]", "Electrolyzer production SCALED [ton]", 'Reformer production [ton]', 'Reformer production SCALED [ton]', 'Emissions per ton [ton CO2/ton H2]']
		writer.writerow(my_header)
		for n in instance.HydrogenProdNode:
			for i in instance.Period:
				for w in instance.Scenario:
					for (s,h) in instance.HoursOfSeason:
						my_string = [n, inv_per[int(i-1)], w, s, h,
									 value(sum(instance.powerForHydrogen[n,j,h,i,w] for j in instance.Period if j<=i)),
									 value(instance.hydrogenProducedElectro[n,h,i,w]),
									 value(instance.seasScale[s] * instance.hydrogenProducedElectro[n,h,i,w])]
						# power_emissions_kg_per_MWh = calculatePowerEmissionIntensity(n,h,i,w) # Commented out this calculation because it's not needed.
						if n in instance.ReformerLocations:
							blue_h2_production_ton = value(sum(instance.hydrogenProducedReformer_ton[n,p,h,i,w] for p in instance.ReformerPlants))
							blue_h2_direct_emissions_ton = value(sum(instance.ReformerEmissionFactor[p,i] * instance.hydrogenProducedReformer_ton[n,p,h,i,w] for p in instance.ReformerPlants))
							# blue_h2_emissions_from_power_ton = power_emissions_kg_per_MWh/1000 * value(sum(instance.ReformerPlantElectricityUse[p,i] * instance.hydrogenProducedReformer_ton[n,p,h,i,w] for p in instance.ReformerPlants))
							blue_h2_emissions_from_power_ton = 0 # Emissions from use of power is 0, because marginal emissions from a sector (here, there power sector) which is already capped on emissions is 0.
						else:
							blue_h2_production_ton = 0
							blue_h2_direct_emissions_ton = 0
							blue_h2_emissions_from_power_ton = 0
						my_string.extend([blue_h2_production_ton,
										  value(instance.seasScale[s] * blue_h2_production_ton)])

						green_h2_production_ton = value(instance.hydrogenProducedElectro[n,h,i,w])
						# green_h2_emissions_ton = power_emissions_kg_per_MWh/1000 * value(sum(instance.powerForHydrogen[n,j,h,i,w] for j in instance.Period if j<=i))
						green_h2_emissions_ton = 0 # Emissions from green H2 is 0, because marginal emissions from a sector (here, there power sector) which is already capped on emissions is 0.
						total_h2_production = blue_h2_production_ton + green_h2_production_ton
						if total_h2_production < .5:
							total_h2_emissions = 0
							my_string.extend([0])
						else:
							total_h2_emissions = blue_h2_direct_emissions_ton + blue_h2_emissions_from_power_ton + green_h2_emissions_ton
							my_string.extend([total_h2_emissions / total_h2_production])

						writer.writerow(my_string)
		f.close()

		print("{hour}:{minute}:{second}: Writing hydrogen sales results to results_hydrogen_use.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_hydrogen_use.csv', 'w', newline='')
		if solver == 'Xpress':
			fError = open(result_file_path + "/" + "errorLog.log",'a')
		writer = csv.writer(f)
		my_header = ["Node", "Period", "Scenario", "Season", "Hour", "Hydrogen stored [ton]", "Hydrogen withdrawn from storage [ton]", "Total demand in period [ton]", "Average hourly demand [ton]", "Demand shed [ton]", "Demand met [ton]", "Demand met SCALED [ton]", "Hydrogen burned for power [ton]", "Hydrogen price [EUR]", 'Hydrogen exported [ton]', 'Hydrogen imported [ton]']
		writer.writerow(my_header)
		if h2storage is True:
			for n in instance.HydrogenProdNode:
				for i in instance.Period:
					for w in instance.Scenario:
						for (s, h) in instance.HoursOfSeason:
							try:
								dualVar = value(instance.dual[instance.hydrogen_flow_balance[n,h,i,w]])
							except:
								# print('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_hydrogen_use.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+')')
								# fError.write('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_hydrogen_use.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+'). Dual is likely 0 (poorly handled in Xpress)\n')
								dualVar = 0
							my_string = [n, inv_per[int(i - 1)], w, s, h,
										 value(instance.hydrogenChargeStorage[n,h,i,w]),
										 value(instance.hydrogenDischargeStorage[n,h,i,w]),
										 value(instance.hydrogenDemand[n,i,w])*8760,
										 value(instance.hydrogenDemand[n,i,w]),
										 value(instance.h2_demand_shed[n,h,i,w]),
										 value(instance.meetHydrogenDemand[n,h,i,w]),
										 value(instance.seasScale[s] * instance.meetHydrogenDemand[n,h,i,w]),
										 value(sum(instance.hydrogenForPower[g,n,h,i,w] for g in instance.HydrogenGenerators)),
										 dualVar,
										 value(sum(instance.hydrogenSentPipeline[n,n2,t,h,i,w] for n2 in instance.HydrogenLinks[n] for t in instance.PipelineType if (n,n2,t) in instance.PipelineTypeOfDirectionalLink)),
										 value(sum(instance.hydrogenSentPipeline[n2,n,t,h,i,w] for n2 in instance.HydrogenLinks[n] for t in instance.PipelineType if (n2,n,t) in instance.PipelineTypeOfDirectionalLink))]
							writer.writerow(my_string)
			if solver == 'Xpress':
				fError.write('\n')
				fError.close()
		else:
			for n in instance.HydrogenProdNode:
				for i in instance.Period:
					for w in instance.Scenario:
						for (s, h) in instance.HoursOfSeason:
							try:
								dualVar = value(instance.dual[instance.hydrogen_flow_balance[n,h,i,w]])
							except:
								# print('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_hydrogen_use.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+')')
								# fError.write('Something went wrong when accessing dual for hydrogen_flow_balance with key for results_hydrogen_use.csv. Key (node,period,scenario)=('+n+','+str(i)+','+w+'). Dual is likely 0 (poorly handled in Xpress)\n')
								dualVar = 0
							my_string = [n, inv_per[int(i - 1)], w, s, h,
									 0,
									 0,
									 value(instance.hydrogenDemand[n,i,w]),
									 value(sum(instance.hydrogenForPower[g,n,h,i,w] for g in instance.HydrogenGenerators)),
									 dualVar,
									 value(sum(instance.hydrogenSentPipeline[n,n2,t,h,i,w] for n2 in instance.HydrogenLinks[n] for t in instance.PipelineType if (n,n2,t) in instance.PipelineTypeOfDirectionalLink))]
							writer.writerow(my_string)
			if solver == 'Xpress':
				fError.write('\n')
				fError.close()
		f.close()

		print("{hour}:{minute}:{second}: Writing hydrogen pipeline investment results to results_hydrogen_pipeline_investments.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_hydrogen_pipeline_investments.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["Between node", "And node", "Pipeline type", "Period", "Pipelines built [#]", "Pipeline total capacity [ton]",
					 "Discounted cost of (newly) built pipeline [EUR]", "Expected hydrogen transmission [tons]"]
		writer.writerow(my_header)
		for (n1,n2) in instance.HydrogenBidirectionPipelines:
			for i in instance.Period:
				for t in instance.PipelineType:
					if (n1,n2,t) in instance.PipelineTypeOfDirectionalLink:
						my_string = [n1, n2, t, inv_per[int(i-1)],
									 value(instance.hydrogenPipelineBuilt[n1,n2,t,i]),
									 value(instance.totalHydrogenPipelineCapacity[n1,n2,t,i]),
									 value(instance.discount_multiplier[i] * (instance.hydrogenPipelineBuilt[n1,n2,t,i] * instance.hydrogenPipelineInvCost[n1,n2,t])),
									 value(sum(instance.sceProbab[w]*instance.seasScale[s]*(instance.hydrogenSentPipeline[n1,n2,t,h,i,w] + instance.hydrogenSentPipeline[n2,n1,t,h,i,w]) for (s,h) in instance.HoursOfSeason for w in instance.Scenario))]
						writer.writerow(my_string)
		f.close()

		print(
			"{hour}:{minute}:{second}: Writing hydrogen pipeline operational results to results_hydrogen_pipeline_operational.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_hydrogen_pipeline_operational.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["From node", "To node", "PipelineType", "Period", "Season", "Scenario", "Hour", "Hydrogen sent [ton]", "Power consumed in each node for transport (MWh)"]
		writer.writerow(my_header)
		for (n1,n2) in instance.AllowedHydrogenLinks:
			if (n1,n2) in instance.HydrogenBidirectionPipelines:
				for t in instance.PipelineType:
					if (n1,n2,t) in instance.PipelineTypeOfDirectionalLink:
						for i in instance.Period:
							for (s,h) in instance.HoursOfSeason:
								for w in instance.Scenario:
									my_string = [n1,n2,t,inv_per[int(i-1)],s,w,h,
												 value(instance.hydrogenSentPipeline[n1,n2,t,h,i,w]),
												 value(0.5*(instance.hydrogenSentPipeline[n1,n2,t,h,i,w] * instance.hydrogenPipelinePowerDemandPerTon[n1,n2,t]))]
									writer.writerow(my_string)
			else:
				for i in instance.Period:
					for t in instance.PipelineType:
						if (n1,n2,t) in instance.PipelineTypeOfDirectionalLink:
							for (s,h) in instance.HoursOfSeason:
								for w in instance.Scenario:
									my_string = [n1,n2,t,inv_per[int(i-1)],s,w,h,
												 value(instance.hydrogenSentPipeline[n1,n2,t,h,i,w]),
												 value(0.5*(instance.hydrogenSentPipeline[n1,n2,t,h,i,w] * instance.hydrogenPipelinePowerDemandPerTon[n2,n1,t]))]
									writer.writerow(my_string)
		f.close()

		print("{hour}:{minute}:{second}: Writing CO2 pipeline investment results to results_CO2_pipeline_investments.csv...".format(
			hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_CO2_pipeline_investments.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["Between node", "And node", "Pipeline type", "Period", "Pipelines built [#]", "Pipeline total capacity [ton]",
					 "Discounted cost of (newly) built pipeline [EUR]", "Expected CO2 transmission [tons]"]
		writer.writerow(my_header)
		for (n1,n2) in instance.CO2BidirectionalPipelines:
			for i in instance.Period:
				for t in instance.PipelineType:
					if (n1,n2,t) in instance.PipelineTypeOfDirectionalLink:
						my_string = [n1, n2, t, inv_per[int(i-1)],
									 value(instance.CO2PipelineBuilt[n1,n2,t,i]),
									 value(instance.totalCO2PipelineCapacity[n1,n2,t,i]),
									 value(instance.discount_multiplier[i] * (instance.CO2PipelineBuilt[n1,n2,t,i] * instance.CO2PipelineInvCost[n1,n2,t,i])),
									 value(sum(instance.sceProbab[w]*instance.seasScale[s]*(instance.CO2sentPipeline[n1,n2,t,h,i,w] + instance.CO2sentPipeline[n2,n1,t,h,i,w]) for (s,h) in instance.HoursOfSeason for w in instance.Scenario))]
						writer.writerow(my_string)
		f.close()

		print("{hour}:{minute}:{second}: Writing CO2 pipeline operational results to results_CO2_pipeline_operational.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_CO2_pipeline_operational.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["From node", "To node", "PipelineType", "Period", "Season", "Scenario", "Hour", "CO2 sent [ton]", "Power consumed in each node for transport (MWh)"]
		writer.writerow(my_header)
		for (n1,n2) in instance.CO2DirectionalLinks:
			if (n1,n2) in instance.CO2BidirectionalPipelines:
				for t in instance.PipelineType:
					if (n1,n2,t) in instance.PipelineTypeOfDirectionalLink:
						for i in instance.Period:
							for (s,h) in instance.HoursOfSeason:
								for w in instance.Scenario:
									my_string = [n1,n2,t,inv_per[int(i-1)],s,w,h,
												 value(instance.CO2sentPipeline[n1,n2,t,h,i,w]),
												 value(0.5*(instance.CO2sentPipeline[n1,n2,t,h,i,w] * instance.CO2PipelinePowerDemandPerTon[n1,n2,t]))]
									writer.writerow(my_string)
			else:
				for i in instance.Period:
					for t in instance.PipelineType:
						if (n1,n2,t) in instance.PipelineTypeOfDirectionalLink:
							for (s,h) in instance.HoursOfSeason:
								for w in instance.Scenario:
									my_string = [n1,n2,t,inv_per[int(i-1)],s,w,h,
												 value(instance.CO2sentPipeline[n1,n2,t,h,i,w]),
												 value(0.5*(instance.CO2sentPipeline[n1,n2,t,h,i,w] * instance.CO2PipelinePowerDemandPerTon[n2,n1,t]))]
									writer.writerow(my_string)
		f.close()

		print("{hour}:{minute}:{second}: Writing CO2 sequestration investment results to results_CO2_sequestration_investments.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_CO2_sequestration_investments.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["Node","Period", "CO2 sequestration capacity built [ton/hr]", "Total CO2 sequestration capacity[ton/hr]", "Sequestered in period (scaled) [Mton]", "Total cumulative amount of CO2 sequestered (scaled)[Mton]"]
		writer.writerow(my_header)
		for n in instance.CO2SequestrationNodes:
			for i in instance.Period:
				writer_string = [n,inv_per[int(i-1)],
								 value(instance.CO2SiteCapacityDeveloped[n,i]),
								 value(sum(instance.CO2SiteCapacityDeveloped[n,j] for j in instance.Period if j<=i)),
								 value(sum(instance.sceProbab[w]*instance.seasScale[s]*instance.CO2sequestered[n,h,i,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario))/1e6,
								 value(sum(instance.sceProbab[w]*instance.seasScale[s]*instance.CO2sequestered[n,h,j,w] for (s,h) in instance.HoursOfSeason for w in instance.Scenario for j in instance.Period if j<=i))/1e6]
				writer.writerow(writer_string)
		f.close()

		print("{hour}:{minute}:{second}: Writing CO2 sequestration results to results_CO2_sequestration_operational.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_CO2_sequestration_operational.csv', 'w', newline='')
		writer = csv.writer(f)
		my_header = ["Node","Period", "Season", "Scenario", "Hour", "CO2 sequestered [ton]"]
		writer.writerow(my_header)
		for n in instance.CO2SequestrationNodes:
			for i in instance.Period:
				for (s,h) in instance.HoursOfSeason:
					for w in instance.Scenario:
						writer_string = [n,inv_per[int(i-1)],s,w,h,
							   value(instance.CO2sequestered[n,h,i,w])]
						writer.writerow(writer_string)
		f.close()

		print("{hour}:{minute}:{second}: Writing CO2 flow balance to results_CO2_flow_balance.csv...".format(
				hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		f = open(result_file_path + '/' + 'results_CO2_flow_balance.csv', 'w', newline='')
		writer = csv.writer(f)
		# my_header = ["Node","Period", "Season", "Scenario", "Hour", "CO2 captured from power generators [ton]", "CO2 captured from natural gas reformers [ton]", "CO2 exported by pipeline [ton]", "CO2 imported by pipeline [ton]", "CO2 liquefied [ton]", "CO2 regasified [ton]", "CO2 exported by ship [ton]", "CO2 imported by ship [ton]", "Liquid storage charge [ton]", "Liquid storage discharge [ton]", "CO2 sequested [ton]"]
		my_header = ["Node","Period", "Season", "Scenario", "Hour", "CO2 captured from power generators [ton]", "CO2 captured from natural gas reformers [ton]", "CO2 exported by pipeline [ton]", "CO2 imported by pipeline [ton]", "CO2 sequested [ton]"]
		writer.writerow(my_header)
		for n in instance.OnshoreNodes:
			for i in instance.Period:
				for (s,h) in instance.HoursOfSeason:
					for w in instance.Scenario:
						writer_string = [n,inv_per[int(i-1)],s,w,h]
						writer_string.extend([value(instance.co2_captured_generators[n,h,i,w]), value(instance.co2_captured_reformers[n,h,i,w])])
						writer_string.extend([value(sum(instance.CO2sentPipeline[n,n2,t,h,i,w] for n2 in instance.OnshoreNodes if (n,n2) in instance.CO2DirectionalLinks for t in instance.PipelineType if (n,n2,t) in instance.PipelineTypeOfDirectionalLink))])
						writer_string.extend([value(sum(instance.CO2sentPipeline[n2,n,t,h,i,w] for n2 in instance.OnshoreNodes if (n,n2) in instance.CO2DirectionalLinks for t in instance.PipelineType if (n,n2,t) in instance.PipelineTypeOfDirectionalLink))])
						# writer_string.extend([value(instance.CO2Liquefied[n,h,i,w])])
						# writer_string.extend([value(instance.CO2Regasified[n,h,i,w])])
						# ship_arrival = 0
						# for n2 in instance.OnshoreNodes:
						# 	if (n2,n) in instance.CO2DirectionalLinks:
						# 		half_shipping_time = ceil(value(instance.ShippingTimes[n,n2]) / 2)
						# 		if h < min(instance.FirstHoursOfPeakSeason): # We are in normal seasons
						# 			season_starting_hour = 1 + floor((h-1)/lengthRegSeason) * lengthRegSeason
						# 			starting_hour_receiving_ship = h - half_shipping_time
						# 			while starting_hour_receiving_ship < season_starting_hour:
						# 				starting_hour_receiving_ship += lengthRegSeason
						# 		else: # We are in peak seasons
						# 			season_starting_hour = 1 + floor((h - 1 - NoOfRegSeason * lengthRegSeason)/lengthPeakSeason) * lengthPeakSeason + NoOfRegSeason * lengthRegSeason
						# 			starting_hour_receiving_ship = h-half_shipping_time
						# 			while starting_hour_receiving_ship < season_starting_hour:
						# 				starting_hour_receiving_ship += lengthPeakSeason
						# 		ship_arrival += value(instance.CO2Shipped[n2,n,starting_hour_receiving_ship,i,w])
						# writer_string.extend([value(sum(instance.CO2Shipped[n,n2,h,i,w] for n2 in instance.OnshoreNodes if (n,n2) in instance.CO2DirectionalLinks))])
						# writer_string.extend([ship_arrival])
						# writer_string.extend([value(instance.CO2LiquidStorageCharge[n,h,i,w])])
						# writer_string.extend([value(instance.CO2LiquidStorageDischarge[n,h,i,w])])
						if n in instance.CO2SequestrationNodes:
							writer_string.extend([value(instance.CO2sequestered[n,h,i,w])])
						else:
							writer_string.extend([0])
						writer.writerow(writer_string)
		f.close()

		# print("{hour}:{minute}:{second}: Writing liquid CO2 ship investment results to results_liquid_CO2_ship_investments.csv...".format(
		# 		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_liquid_CO2_ship_investments.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ["From Node", "To Node", "Period", "Liquid CO2 ships bought [#]", "Added liquid CO2 ship capacity [tons]", "Total liquid CO2 ship capacity [tons]"]
		# writer.writerow(my_header)
		# for n in instance.OnshoreNodes:
		# 	for n2 in instance.OnshoreNodes:
		# 		if (n,n2) in instance.CO2DirectionalLinks:
		# 			for i in instance.Period:
		# 				write_str = [n,n2,inv_per[i-1],
		# 					   value(instance.CO2ShipBought[n,n2,i]),
		# 					   value(instance.CO2ShipBought[n,n2,i] * instance.CO2LiquidShipCapacity),
		# 					   value(instance.totalCO2ShipCapacity[n,n2,i])]
		# 				writer.writerow(write_str)
		# f.close()
		#
		# print("{hour}:{minute}:{second}: Writing liquid CO2 storage investment results to results_liquid_CO2_storage_investments.csv...".format(
		# 		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_liquid_CO2_storage_investments.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ["Node", "Period", "Storage built [tons]", "Total storage capacity [tons]"]
		# writer.writerow(my_header)
		# for n in instance.OnshoreNodes:
		# 	for i in instance.Period:
		# 		write_str = [n, inv_per[i-1],
		# 					 value(instance.CO2LiquidStorageBuilt[n,i]),
		# 					 value(instance.CO2LiquidTotalStorage[n,i])]
		# 		writer.writerow(write_str)
		# f.close()
		#
		# print("{hour}:{minute}:{second}: Writing CO2 liquefier investment results to results_liquid_CO2_liquefier_investments.csv...".format(
		# 		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_liquid_CO2_liquefier_investments.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ["Node", "Period", "Liquefier capacity built [tons/hr]", "Total liquefier capacity [tons/hr]"]
		# writer.writerow(my_header)
		# for n in instance.OnshoreNodes:
		# 	for i in instance.Period:
		# 		write_str = [n, inv_per[i-1],
		# 					 value(instance.CO2LiquefierCapacityBuilt[n,i]),
		# 					 value(instance.CO2LiquefierTotalCapacity[n,i])]
		# 		writer.writerow(write_str)
		# f.close()
		#
		# print("{hour}:{minute}:{second}: Writing liquid CO2 ship operational results to results_liquid_CO2_ship_operational.csv...".format(
		# 		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_liquid_CO2_ship_operational.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ["From Node", "To Node", "Period", "Hour", "Scenario", "CO2 sent via ships [tons]", "Arrival hour [hour]", "CO2 arrived via ships [tons]"]
		# writer.writerow(my_header)
		# for n in instance.OnshoreNodes:
		# 	for n2 in instance.OnshoreNodes:
		# 		if (n,n2) in instance.CO2DirectionalLinks:
		# 			half_shipping_time = ceil(value(instance.ShippingTimes[n,n2]) / 2)
		# 			for i in instance.Period:
		# 				for h in instance.Operationalhour:
		# 					if h < min(instance.FirstHoursOfPeakSeason): # We are in normal seasons
		# 						season_starting_hour = 1 + floor((h-1)/lengthRegSeason) * lengthRegSeason
		# 						starting_hour_receiving_ship = h - half_shipping_time
		# 						while starting_hour_receiving_ship < season_starting_hour:
		# 							starting_hour_receiving_ship += lengthRegSeason
		# 						for w in instance.Scenario:
		# 							write_str = [n,n2,inv_per[i-1],h,w,
		# 								   		value(instance.CO2Shipped[n,n2,h,i,w]),
		# 									 	(h - season_starting_hour + half_shipping_time) % lengthRegSeason + season_starting_hour,
		# 								   		value(instance.CO2Shipped[n2,n,starting_hour_receiving_ship,i,w])]
		# 							writer.writerow(write_str)
		# 					else:
		# 						season_starting_hour = 1 + floor((h - 1 - NoOfRegSeason * lengthRegSeason)/lengthPeakSeason) * lengthPeakSeason + NoOfRegSeason * lengthRegSeason
		# 						starting_hour_receiving_ship = h-half_shipping_time
		# 						while starting_hour_receiving_ship < season_starting_hour:
		# 							starting_hour_receiving_ship += lengthPeakSeason
		# 						for w in instance.Scenario:
		# 							write_str = [n,n2,inv_per[i-1],h,w,
		# 								   		value(instance.CO2Shipped[n,n2,h,i,w]),
		# 										(h - season_starting_hour + half_shipping_time) % lengthPeakSeason + season_starting_hour,
		# 								   		value(instance.CO2Shipped[n2,n,starting_hour_receiving_ship,i,w])]
		# 							writer.writerow(write_str)
		# f.close()
		#
		# print("{hour}:{minute}:{second}: Writing liquid CO2 storage operational results to results_liquid_CO2_storage_operational.csv...".format(
		# 		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_liquid_CO2_storage_operational.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ["Node", "Period", "Hour","Scenario", "Storage charge [tons]", "Storage discharge[tons]", "Storage operational level [tons]"]
		# writer.writerow(my_header)
		# for n in instance.OnshoreNodes:
		# 	for i in instance.Period:
		# 		for h in instance.Operationalhour:
		# 			for w in instance.Scenario:
		# 				my_str = [n, inv_per[i-1], h, w,
		# 						  value(instance.CO2LiquidStorageCharge[n,h,i,w]),
		# 						  value(instance.CO2LiquidStorageDischarge[n,h,i,w]),
		# 						  value(instance.CO2LiquidStorageOperational[n,h,i,w])]
		# 				writer.writerow(my_str)
		# f.close()
		#
		# print("{hour}:{minute}:{second}: Writing liquid CO2 liquefier & gasifier operational results to results_liquid_CO2_liquefier_regasifier_operational.csv...".format(
		# 		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_liquid_CO2_liquefier_regasifier_operational.csv', 'w', newline='')
		# writer = csv.writer(f)
		# my_header = ["Node", "Period", "Hour","Scenario", "CO2 liquefied [tons]", "CO2 regasified [tons]"]
		# writer.writerow(my_header)
		# for n in instance.OnshoreNodes:
		# 	for i in instance.Period:
		# 		for h in instance.Operationalhour:
		# 			for w in instance.Scenario:
		# 				my_str = [n, inv_per[i-1], h, w,
		# 						  value(instance.CO2Liquefied[n,h,i,w]),
		# 						  value(instance.CO2Regasified[n,h,i,w])]
		# 				writer.writerow(my_str)
		# f.close()

		# print("{hour}:{minute}:{second}: Writing hydrogen investments costs and NPV calculation to results_hydrogen_costs.csv...".format(
		# 	hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'results_hydrogen_costs.csv', 'w', newline='')
		# writer = csv.writer(f)
		# header = ['Period','Discounted electrolyzer cost [EUR]', 'Discounted Reformer cost [EUR]', 'Discounted pipeline cost [EUR]', 'Discounted storage cost [EUR]',
		# 		  'Total discounted cost [EUR]', 'Hydrogen sold [kg]','Non-discounted price for NPV = 0 [EUR]']
		# writer.writerow(header)
		# if h2storage is True:
		# 	for i in instance.Period:
		# 		electrolyzerCost = value(sum(instance.discount_multiplier[j] * instance.elyzerInvCost[j] * sum(instance.elyzerCapBuilt[n,j,j] for n in instance.HydrogenProdNode) for j in instance.Period if j<=i))
		# 		reformerCost = value(sum(instance.discount_multiplier[j] * sum(instance.ReformerPlantInvCost[p,j] * instance.ReformerCapBuilt[n,p,j] for n in instance.ReformerLocations for p in instance.ReformerPlants) for j in instance.Period if j<=i))
		# 		pipelineCost = value(sum(instance.discount_multiplier[j] * sum(instance.hydrogenPipelineInvCost[(n1,n2),j] * instance.hydrogenPipelineBuilt[(n1,n2),j] for (n1,n2) in instance.HydrogenBidirectionPipelines) for j in instance.Period if j<=i))
		# 		storageCost = value(sum(instance.discount_multiplier[j] * sum(instance.hydrogenStorageBuilt[n,j] for n in instance.HydrogenProdNode) for j in instance.Period if j<=i))
		# 		my_string = [inv_per[i-1],
		# 					 electrolyzerCost,
		# 					 reformerCost,
		# 					 pipelineCost,
		# 					 storageCost,
		# 					 electrolyzerCost + reformerCost + pipelineCost + storageCost]
		# 		if value(sum(instance.hydrogenDemand[n,i] for n in instance.HydrogenProdNode)) > 0:
		# 			my_string.extend([(electrolyzerCost + reformerCost + pipelineCost + storageCost)/value(sum(instance.discount_multiplier[i] * instance.hydrogenDemand[n,i] for n in instance.HydrogenProdNode))])
		# 		else:
		# 			my_string.extend([0])
		# 		writer.writerow(my_string)
		# else:
		# 	for i in instance.Period:
		# 		electrolyzerCost = value(sum(instance.discount_multiplier[j] * instance.elyzerInvCost[j] * sum(instance.elyzerCapBuilt[n,j,j] for n in instance.HydrogenProdNode) for j in instance.Period if j<=i))
		# 		reformerCost = value(sum(instance.discount_multiplier[j] * sum(instance.ReformerPlantInvCost[p,j] * instance.ReformerCapBuilt[n,p,j] for n in instance.ReformerLocations for p in instance.ReformerPlants) for j in instance.Period if j<=i))
		# 		pipelineCost = value(sum(instance.discount_multiplier[j] * sum(instance.hydrogenPipelineInvCost[(n1,n2),j] * instance.hydrogenPipelineBuilt[(n1,n2),j] for (n1,n2) in instance.HydrogenBidirectionPipelines) for j in instance.Period if j<=i))
		# 		storageCost = 0
		# 		my_string = [inv_per[i-1],
		# 					 electrolyzerCost,
		# 					 reformerCost,
		# 					 pipelineCost,
		# 					 storageCost,
		# 					 electrolyzerCost + reformerCost + pipelineCost + storageCost]
		# 		if value(sum(instance.hydrogenDemand[n,i] for n in instance.HydrogenProdNode)) > 0:
		# 			my_string.extend([(electrolyzerCost + reformerCost + pipelineCost + storageCost)/value(sum(instance.discount_multiplier[i] * instance.hydrogenDemand[n,i] for n in instance.HydrogenProdNode))])
		# 		else:
		# 			my_string.extend([0])
		# 		writer.writerow(my_string)
		# f.close()

		# This reporting was written to track the hydrogenForPower variable. To be deleted later.

		# print("{hour}:{minute}:{second}: Writing hydrogen use in as power fuel in hydrogen_for_power.csv...".format(
		# 	hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))
		# f = open(result_file_path + '/' + 'hydrogen_for_power.csv', 'w', newline='')
		# writer = csv.writer(f)
		# header = ['Node','Generator','Period','Hour','Scenario']
		# writer.writerow(header)
		# for (n,g) in instance.GeneratorsOfNode:
		# 	for i in instance.Period:
		# 		for (s,h) in instance.HoursOfSeason:
		# 			for w in instance.Scenario:
		# 				writer.writerow([n,g,inv_per[i-1],h,w,
		# 								 value(instance.hydrogenForPower[g,n,h,i,w])])

	endReporting = timeEnd = datetime.now()

	print("{hour}:{minute}:{second}: Writing time usage to time_usage.csv...".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))

	f = open(result_file_path + "/" + 'time_usage.csv', 'w', newline='')
	timeFrmt = "%H:%M:%S"
	dateFrmt = "%d.%m.%Y"
	timeDeltaFrmt = "{H}:{M}:{S}"
	writer = csv.writer(f)
	if (timeEnd - timeStart).days > 0:
		writer.writerow(["Process",
						 "Time started [HH:MM:SS]",
						 "Time ended [HH:MM:SS]",
						 "Time spent [HH:MM:SS]",
						 "Date started [DD.MM.YYYY]",
						 "Date finished [DD.MM.YYYY]"])
		writer.writerow(["Overall",
						 timeStart.strftime(timeFrmt),
						 timeEnd.strftime(timeFrmt),
						 strfdelta(timeEnd-timeStart,timeDeltaFrmt),
						 timeStart.strftime(dateFrmt),
						 timeEnd.strftime(dateFrmt)])
		writer.writerow(["Declaring and reading sets & parameters",
						 timeStart.strftime(timeFrmt),
						 stopReading.strftime(timeFrmt),
						 strfdelta(stopReading-timeStart,timeDeltaFrmt),
						 timeStart.strftime(dateFrmt),
						 stopReading.strftime(dateFrmt)])
		writer.writerow(["Declaring variables & constraints",
						 startConstraints.strftime(timeFrmt),
						 stopConstraints.strftime(timeFrmt),
						 strfdelta(stopConstraints-startConstraints,timeDeltaFrmt),
						 startConstraints.strftime(dateFrmt),
						 stopConstraints.strftime(dateFrmt)])
		writer.writerow(["Building model",
						 startBuild.strftime(timeFrmt),
						 endBuild.strftime(timeFrmt),
						 strfdelta(endBuild-startBuild,timeDeltaFrmt),
						 startBuild.strftime(dateFrmt),
						 endBuild.strftime(dateFrmt)])
		writer.writerow(["Optimizing model",
						 startOptimization.strftime(timeFrmt),
						 endOptimization.strftime(timeFrmt),
						 strfdelta(endOptimization-startOptimization,timeDeltaFrmt),
						 startOptimization.strftime(dateFrmt),
						 endOptimization.strftime(dateFrmt)])
		writer.writerow(["Reporting results",
						 StartReporting.strftime(timeFrmt),
						 endReporting.strftime(timeFrmt),
						 strfdelta(endReporting-StartReporting,timeDeltaFrmt),
						 StartReporting.strftime(dateFrmt),
						 endReporting.strftime(dateFrmt)])
	else:
		writer.writerow(["Process",
						 "Time started [HH:MM:SS]",
						 "Time ended [HH:MM:SS]",
						 "Time spent [HH:MM:SS]"])
		writer.writerow(["Overall",
						 timeStart.strftime(timeFrmt),
						 timeEnd.strftime(timeFrmt),
						 strfdelta(timeEnd - timeStart, timeDeltaFrmt)])
		writer.writerow(["Declaring and reading sets & parameters",
						 timeStart.strftime(timeFrmt),
						 stopReading.strftime(timeFrmt),
						 strfdelta(stopReading - timeStart, timeDeltaFrmt)])
		writer.writerow(["Declaring variables & constraints",
						 startConstraints.strftime(timeFrmt),
						 stopConstraints.strftime(timeFrmt),
						 strfdelta(stopConstraints - startConstraints, timeDeltaFrmt)])
		writer.writerow(["Building model",
						 startBuild.strftime(timeFrmt),
						 endBuild.strftime(timeFrmt),
						 strfdelta(endBuild - startBuild, timeDeltaFrmt)])
		writer.writerow(["Optimizing model",
						 startOptimization.strftime(timeFrmt),
						 endOptimization.strftime(timeFrmt),
						 strfdelta(endOptimization - startOptimization,
								   timeDeltaFrmt)])
		writer.writerow(["Reporting results",
						 StartReporting.strftime(timeFrmt),
						 endReporting.strftime(timeFrmt),
						 strfdelta(endReporting - StartReporting, timeDeltaFrmt)])
	f.close()



	print("{hour}:{minute}:{second} Finished writing results to files.".format(
		hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")))

	del results, instance, model

	# for i in instance.Period:
	#     custom_lines = [Line2D([0], [0], color='black', linewidth=0.5),
	#                     Line2D([0], [0], color='black', linewidth=2.5),
	#                     Line2D([0], [0], color='black', linewidth=5),
	#                     Line2D([0], [0], color='black', linewidth=7.5),
	#                     Line2D([0], [0], color='black', linewidth=10)]
	#     fig = plt.figure(figsize=(20,20))
	#     ax = plt.axes(projection=ccrs.Orthographic())
	#     ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':', alpha=.5)
	#     ax.add_feature(cartopy.feature.LAND.with_scale('50m'))
	#     ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
	#     for (n1,n2) in instance.BidirectionalArc:
	#         plt.plot([value(instance.Longitude[n1]),value(instance.Longitude[n2])], [value(instance.Latitude[n1]), value(instance.Latitude[n2])],
	#                   color = 'black',
	#                   linewidth = value(instance.transmissionInstalledCap[n1,n2,i]/2000),
	#                   marker = 'o',
	#                   transform=ccrs.Geodetic())
	#     ax.set_extent([-5, 11, 50, 65], crs=ccrs.PlateCarree())
	#     ax.legend(custom_lines, ['  1 GW', '  5 GW', '  10 GW','  15 GW','  20 GW'], borderpad=2, labelspacing=2, handlelength=7)
	#     plt.title('Transmission capacity North Sea ' + inv_per[int(i-1)])
	#     plt.savefig(result_file_path + "/" + 'TransmissionCap_NorthSea_' + inv_per[int(i-1)] + '.png')
	#
	#
	# for i in instance.Period:
	#     custom_lines = [Line2D([0], [0], color='black', linewidth=0.5),
	#                     Line2D([0], [0], color='black', linewidth=2.5),
	#                     Line2D([0], [0], color='black', linewidth=5),
	#                     Line2D([0], [0], color='black', linewidth=7.5),
	#                     Line2D([0], [0], color='black', linewidth=10)]
	#     fig = plt.figure(figsize=(20,20))
	#     ax = plt.axes(projection=ccrs.Orthographic())
	#     ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':', alpha=.5)
	#     ax.add_feature(cartopy.feature.LAND.with_scale('50m'))
	#     ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
	#     for (n1,n2) in instance.BidirectionalArc:
	#         plt.plot([value(instance.Longitude[n1]),value(instance.Longitude[n2])], [value(instance.Latitude[n1]), value(instance.Latitude[n2])],
	#                   color = 'black',
	#                   linewidth = value(instance.transmissionInstalledCap[n1,n2,i]/2000),
	#                   marker = 'o',
	#                   transform=ccrs.Geodetic())
	#     ax.set_extent([-9, 25, 36, 72], crs=ccrs.PlateCarree())
	#     ax.legend(custom_lines, ['  1 GW', '  5 GW', '  10 GW','  15 GW','  20 GW'], borderpad=2, labelspacing=2, handlelength=7)
	#     plt.title('Transmission capacity Europe ' + inv_per[int(i-1)])
	#     plt.savefig(result_file_path + "/" + 'TransmissionCap_Europe_' + inv_per[int(i-1)] + '.png')