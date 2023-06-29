import pandas as pd
import os

def read_file(filepath, excel, sheet, columns, tab_file_path):
    input_sheet = pd.read_excel(filepath + "/" +excel, sheet, skiprows=2)
    if "_noCCS" in excel:
        excel = excel.replace('_noCCS','')
    if '_sequestration_all' in excel:
        excel = excel.replace('_sequestration_all','')

    data_table = input_sheet.iloc[:, columns]
    data_table.columns = pd.Series(data_table.columns).str.replace(' ', '_')
    data_nonempty = data_table.dropna()

    save_csv_frame = pd.DataFrame(data_nonempty)

    save_csv_frame.replace('\s', '', regex=True, inplace=True)

    if not os.path.exists(tab_file_path):
        os.makedirs(tab_file_path)
    #excel = excel.replace(".xlsx", "_")
    #excel = excel.replace("Excel/", "")
    save_csv_frame.to_csv(tab_file_path + "/" + excel.replace(".xlsx", '_') + sheet + '.tab', header=True, index=None, sep='\t', mode='w')
    #save_csv_frame.to_csv(excel.replace(".xlsx", '_') + sheet + '.tab', header=True, index=None, sep='\t', mode='w')

def read_sets(filepath, excel, sheet, tab_file_path):
    input_sheet = pd.read_excel(filepath + "/" + excel, sheet)
    if "_noCCS" in excel:
        excel = excel.replace('_noCCS','')
    if '_sequestration_all' in excel:
        excel = excel.replace('_sequestration_all','')

    for ind, column in enumerate(input_sheet.columns):
        data_table = input_sheet.iloc[:, ind]
        data_nonempty = data_table.dropna()
        data_nonempty.replace(" ", "")
        save_csv_frame = pd.DataFrame(data_nonempty)
        save_csv_frame.replace('\s', '', regex=True, inplace=True)
        if not os.path.exists(tab_file_path):
            os.makedirs(tab_file_path)
        #excel = excel.replace(".xlsx", "_")
        #excel = excel.replace("Excel/", "")
        if 'Unnamed' in column:
            print('\n\n\nWARNING: Unnamed column found in sheet ' + sheet + '\n\n\n')
        save_csv_frame.to_csv(tab_file_path + "/" + excel.replace(".xlsx", '_') + column + '.tab', header=True, index=None, sep='\t', mode='w')
        #save_csv_frame.to_csv(excel.replace(".xlsx", '_') + column + '.tab', header=True, index=None, sep='\t', mode='w')

def generate_tab_files(filepath, tab_file_path, scenariogeneration=True, hydrogen=False, case = ''):
    # Function description: read column value from excel sheet and save as .tab file "sheet.tab"
    # Input: excel name, sheet name, the number of columns to be read
    # Output:  .tab file
    
    print("Generating .tab-files...")

    # Reading Excel workbooks using our function read_file

    if not os.path.exists(tab_file_path):
        os.makedirs(tab_file_path)

    read_sets(filepath, 'Sets.xlsx', 'Nodes',tab_file_path = tab_file_path)
    #read_sets(filepath, 'Sets.xlsx', 'Times', tab_file_path = tab_file_path)
    read_sets(filepath, 'Sets.xlsx', 'LineType', tab_file_path = tab_file_path)
    read_sets(filepath, 'Sets.xlsx', 'PipelineType', tab_file_path = tab_file_path)
    read_sets(filepath, 'Sets.xlsx', 'Technology', tab_file_path = tab_file_path)
    read_sets(filepath, 'Sets.xlsx', 'Storage', tab_file_path = tab_file_path)
    read_sets(filepath, 'Sets.xlsx', 'Generators', tab_file_path = tab_file_path)
    read_file(filepath, 'Sets.xlsx', 'StorageOfNodes', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Sets.xlsx', 'GeneratorsOfNode', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Sets.xlsx', 'GeneratorsOfTechnology', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Sets.xlsx', 'DirectionalLines', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Sets.xlsx', 'PipelineTypeOfLines', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Sets.xlsx', 'LineTypeOfDirectionalLines', [0, 1, 2], tab_file_path = tab_file_path)
    read_sets(filepath, 'Sets.xlsx', 'HydrogenGenerators', tab_file_path = tab_file_path)


    # Reading GeneratorPeriod
    if case == 'noCCS':
        generator_file_name = 'Generator_noCCS.xlsx'
    else:
        generator_file_name = 'Generator.xlsx'
    read_file(filepath, generator_file_name, 'FixedOMCosts', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'CapitalCosts', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'VariableOMCosts', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'FuelCosts', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'CCSCostTSVariable', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'Efficiency', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'RefInitialCap', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'ScaleFactorInitialCap', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'InitialCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'MaxBuiltCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'MaxInstalledCapacity', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'RampRate', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'GeneratorTypeAvailability', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'CO2Content', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'CO2Captured', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, generator_file_name, 'Lifetime', [0, 1], tab_file_path = tab_file_path)

    #Reading InterConnector
    read_file(filepath, 'Transmission.xlsx', 'lineEfficiency', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'MaxInstallCapacityRaw', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'MaxBuiltCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'Length', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'TypeConverterFixedCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'TypeConverterVariableCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'TypeCableFixedCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'TypeCableVariableCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'TypeFixedOMCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'InitialCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Transmission.xlsx', 'Lifetime', [0, 1, 2], tab_file_path = tab_file_path)

    #Reading Node
    read_file(filepath, 'Node.xlsx', 'ElectricAnnualDemand', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Node.xlsx', 'NodeLostLoadCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Node.xlsx', 'HydroGenMaxAnnualProduction', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Node.xlsx', 'Latitude', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Node.xlsx', 'Longitude', [0, 1], tab_file_path = tab_file_path)

    #Reading Season
    read_file(filepath, 'General.xlsx', 'seasonScale', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'General.xlsx', 'CO2Cap', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'General.xlsx', 'CO2Price', [0, 1], tab_file_path = tab_file_path)
    
    #Reading Storage
    read_file(filepath, 'Storage.xlsx', 'StorageBleedEfficiency', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'StorageChargeEff', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'StorageDischargeEff', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'StoragePowToEnergy', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'StorageInitialEnergyLevel', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'InitialPowerCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'PowerCapitalCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'PowerFixedOMCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'PowerMaxBuiltCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'EnergyCapitalCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'EnergyFixedOMCost', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'EnergyInitialCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'EnergyMaxBuiltCapacity', [0, 1, 2, 3], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'EnergyMaxInstalledCapacity', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'PowerMaxInstalledCapacity', [0, 1, 2], tab_file_path = tab_file_path)
    read_file(filepath, 'Storage.xlsx', 'Lifetime', [0, 1], tab_file_path = tab_file_path)

    #Reading CO2 files
    if case == 'sequestration_all':
        co2_file_name = 'CO2_sequestration_all.xlsx'
    else:
        co2_file_name = 'CO2.xlsx'

    read_sets(filepath, co2_file_name, 'CO2SequestrationNodes', tab_file_path = tab_file_path)
    read_file(filepath, co2_file_name, 'StorageSiteCapitalCost', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, co2_file_name, 'StorageSiteFixedOMCost', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, co2_file_name, 'PipelineCapitalCost', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, co2_file_name, 'PipelineFixedOM', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, co2_file_name, 'PipelineCapacity', [0, 1], tab_file_path = tab_file_path)
    read_file(filepath, co2_file_name, 'PipelineElectricityUsage', [0, 1], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquefierCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquefierFixedOMCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquefierElectricityUse', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidStorageCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidStorageFixedOM', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidShipCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidShipFixedOM', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidShipVariableCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidShipCapacity', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidShipLifetime', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'CO2.xlsx', 'LiquidShipLoadingDischargeTime', [0], tab_file_path = tab_file_path)

    #Reading ammonia files
    # read_file(filepath, 'Ammonia.xlsx', 'ShipCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'ShipFixedOMCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'ShipCapacity', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'ProducerCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'ProducerFixedOM', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'ProducerElectricityUse', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'CrackerCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'CrackerFixedOM', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'CrackerElectricityUse', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'StorageCapitalCost', [0], tab_file_path = tab_file_path)
    # read_file(filepath, 'Ammonia.xlsx', 'StorageFixedOM', [0], tab_file_path = tab_file_path)

    if hydrogen is True:
        read_sets(filepath, 'Hydrogen.xlsx', 'ProductionNodes', tab_file_path = tab_file_path)
        # read_file(filepath, 'Hydrogen.xlsx', 'Links', [0,1], tab_file_path = tab_file_path) # Depcreated; The links are now instead defined by the transmission links, but only between the production nodes
        read_sets(filepath, 'Hydrogen.xlsx', 'ReformerLocations', tab_file_path = tab_file_path)
        read_sets(filepath, 'Hydrogen.xlsx', 'ReformerPlants', tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerCapitalCost', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerFixedOMCost', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerVariableOMCost', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerEfficiency', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerElectricityUse', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerLifetime', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerEmissionFactor', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ReformerCO2CaptureFactor', [0,1,2], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ElectrolyzerPlantCapitalCost', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ElectrolyzerFixedOMCost', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ElectrolyzerStackCapitalCost', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ElectrolyzerLifetime', [0], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'ElectrolyzerMWhPerTon', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'PipelineCapitalCost', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'PipelineOMCostPerKM', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'PipelineCapacity', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'PipelineCompressorPowerUsage', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'StorageCapitalCost', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'StorageFixedOMCost', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'StorageMaxCapacity', [0,1], tab_file_path = tab_file_path)
        read_file(filepath, 'Hydrogen.xlsx', 'Demand', [0,1,2], tab_file_path = tab_file_path)
        # read_file(filepath, 'Hydrogen.xlsx', 'Distances', [0,1,2], tab_file_path = tab_file_path) # Depecreated; Distances are now copied from the transmission distances