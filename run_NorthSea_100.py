from reader import generate_tab_files
from Empire import run_empire
from scenario_random import generate_random_scenario
from datetime import datetime
import time
import gc
import numpy as np

########
##USER##
########

USE_TEMP_DIR = True #True/False
temp_dir = '/mnt/beegfs/users/gorand/TempDir'
version = 'north_sea'
NoOfPeriods = 6
NoOfNormalScenarios = 3
NoOfHydrogenScenarios = 3
NoOfScenarios = NoOfNormalScenarios * NoOfHydrogenScenarios
# NoOfScenarios = 4
NoOfRegSeason = 4
lengthRegSeason = 7*24
regular_seasons = ["winter", "spring", "summer", "fall"]
NoOfPeakSeason = 2
lengthPeakSeason = 24
discountrate = 0.05
WACC = 0.05
LeapYearsInvestment = 5
solver = "Gurobi" #"Gurobi" #"CPLEX" #"Xpress"
scenariogeneration = False #True #False
EMISSION_CAP = False #False
WRITE_LP = False #True
PICKLE_INSTANCE = False #True 
hydrogen = True
h2storage = True
TIME_LIMIT = 0

std_dev_percentages = [100]
hydrogen_demand_percentages = np.linspace(50,100,4, endpoint=True)
# hydrogen_demand_percentages = hydrogen_demand_percentages[::-1]
case=''

#######
##RUN##
#######
for std_dev in std_dev_percentages:
    for h2_demand_perc in hydrogen_demand_percentages:
        if NoOfHydrogenScenarios > 1:
                name = f'{version}_h2scen{NoOfHydrogenScenarios}_stddev{std_dev}_h2demandperc{h2_demand_perc:.2f}'
                # name = f'{version}_withoutCCS_h2scen{NoOfHydrogenScenarios}_stddev{std_dev}_h2demandperc{h2_demand_perc:.2f}'
                # name = f'{version}_sequestrationAll_h2scen{NoOfHydrogenScenarios}_stddev{std_dev}_h2demandperc{h2_demand_perc:.2f}'
        else:
            name = f'{version}_Deterministic_h2demandperc{h2_demand_perc:.2f}'
            # name = f'{version}_withoutCCS_Deterministic_h2demandperc{h2_demand_perc:.2f}'
            # name = f'{version}_sequestrationAll_Deterministic_h2demandperc{h2_demand_perc:.2f}'

        workbook_path = 'Data handler/' + version
        tab_file_path = 'Data handler/' + version + '/Tab_Files_' + name
        scenario_data_path = 'Data handler/' + version + '/ScenarioData'
        result_file_path = 'Results/' + version + '/' + name
        FirstHoursOfRegSeason = [lengthRegSeason*i + 1 for i in range(NoOfRegSeason)]
        FirstHoursOfPeakSeason = [lengthRegSeason*NoOfRegSeason + lengthPeakSeason*i + 1 for i in range(NoOfPeakSeason)]
        Period = [i + 1 for i in range(NoOfPeriods)]
        Scenario = ["scenario"+str(i + 1) for i in range(NoOfScenarios)]
        peak_seasons = ['peak'+str(i + 1) for i in range(NoOfPeakSeason)]
        Season = regular_seasons + peak_seasons
        Operationalhour = [i + 1 for i in range(FirstHoursOfPeakSeason[-1] + lengthPeakSeason - 1)]
        HoursOfRegSeason = [(s,h) for s in regular_seasons for h in Operationalhour \
                         if h in list(range(regular_seasons.index(s)*lengthRegSeason+1,
                                       regular_seasons.index(s)*lengthRegSeason+lengthRegSeason+1))]
        HoursOfPeakSeason = [(s,h) for s in peak_seasons for h in Operationalhour \
                             if h in list(range(lengthRegSeason*len(regular_seasons)+ \
                                                peak_seasons.index(s)*lengthPeakSeason+1,
                                                lengthRegSeason*len(regular_seasons)+ \
                                                    peak_seasons.index(s)*lengthPeakSeason+ \
                                                        lengthPeakSeason+1))]
        HoursOfSeason = HoursOfRegSeason + HoursOfPeakSeason
        dict_countries = {"BE": "Belgium", "DE": "Germany", "DK": "Denmark",
                          "GB": "GreatBrit.","NL": "Netherlands", "NO": "Norway",
                          "DB": "DoggerBank", "SEE": "SouthEastEngland", "BS": "Borssele",
                          "HK": "HollandseeKust", "HB": "HelgoländerBucht", "NS": "Nordsøen",
                          "UN": "UtsiraNord", "SN1": "SørligeNordsjøI", "SN2": "SørligeNordsjøII"}
        # offshoreNodesList = ["Energyhub Great Britain", "Energyhub Norway", "Energyhub EU"]
        windfarmNodes = ["Dogger Bank","South East England","Borssele","Hollandsee Kust","Helgoländer Bucht","Nordsøen","Utsira Nord","Sørlige Nordsjø I","Sørlige Nordsjø II"]

        print(f'{datetime.now().strftime("%A")}, {datetime.now().strftime("%d")}. {datetime.now().strftime("%B")}, {datetime.now().strftime("%Y")}')

        print('++++++++')
        print('+EMPIRE+')
        print('++++++++')
        print('Solver: ' + solver)
        print('Scenario Generation: ' + str(scenariogeneration))
        print('++++++++')
        print('ID: ' + name)
        print('++++++++')
        print('Hydrogen: ' + str(hydrogen))
        print('++++++++')


        if scenariogeneration:
            tick = time.time()
            generate_random_scenario(filepath = scenario_data_path,
                                     tab_file_path = tab_file_path,
                                     scenarios = NoOfScenarios,
                                     seasons = regular_seasons,
                                     Periods = NoOfPeriods,
                                     regularSeasonHours = lengthRegSeason,
                                     peakSeasonHours = lengthPeakSeason,
                                     dict_countries = dict_countries)
            tock = time.time()
            print("{hour}:{minute}:{second}: Scenario generation took [sec]:".format(
            hour=datetime.now().strftime("%H"), minute=datetime.now().strftime("%M"), second=datetime.now().strftime("%S")) + str(tock - tick))

        generate_tab_files(filepath = workbook_path, tab_file_path = tab_file_path,
                           scenariogeneration = scenariogeneration, hydrogen = hydrogen, case=case)

        run_empire(name = name,
                   tab_file_path = tab_file_path,
                   result_file_path = result_file_path,
                   scenariogeneration = scenariogeneration,
                   scenario_data_path = scenario_data_path,
                   solver = solver,
                   temp_dir = temp_dir,
                   FirstHoursOfRegSeason = FirstHoursOfRegSeason,
                   FirstHoursOfPeakSeason = FirstHoursOfPeakSeason,
                   lengthRegSeason = lengthRegSeason,
                   lengthPeakSeason = lengthPeakSeason,
                   Period = Period,
                   Operationalhour = Operationalhour,
                   Scenario = Scenario,
                   Season = Season,
                   HoursOfSeason = HoursOfSeason,
                   NoOfNormalScenarios = NoOfNormalScenarios,
                   NoOfHydrogenScenarios = NoOfHydrogenScenarios,
                   discountrate = discountrate,
                   WACC = WACC,
                   LeapYearsInvestment = LeapYearsInvestment,
                   WRITE_LP = WRITE_LP,
                   PICKLE_INSTANCE = PICKLE_INSTANCE,
                   EMISSION_CAP = EMISSION_CAP,
                   USE_TEMP_DIR = USE_TEMP_DIR,
                   NoOfRegSeason = NoOfRegSeason,
                   NoOfPeakSeason = NoOfPeakSeason,
                   verboseResultWriting = False,
                   hydrogen = hydrogen,
                   TIME_LIMIT = TIME_LIMIT,
                   h2storage = h2storage,
                   windfarmNodes = windfarmNodes,
                   hydrogen_demand_percentage = h2_demand_perc/100,
                   std_dev_percentage = std_dev/100)
    gc.collect()