from __future__ import print_function
import copy
import os
import re
import tempfile
import json

from COMPS.Data import Simulation, SimulationFile, QueryCriteria, Configuration, Experiment, Suite

from simtools.SetupParser import SetupParser

campaign_file_regex = re.compile("^campaign_.+\.json$")

# EDIT HERE #######################################################################


SetupParser.default_block = 'HPC'

# EDIT as needed
# The suite of experiments/simulations to copy, modify, and run
# original_suite_id = 'c79aa54e-2578-ea11-a2c5-c4346bcb1550'
original_suite_id = '6e6ff06d-fed6-e811-a2bd-c4346bcb1555'
EXPERIMENT_IDS_TO_RERUN = ['6f6ff06d-fed6-e811-a2bd-c4346bcb1555']  # baseline scenario, only


def modify_config_json(config_filename):
    with open(config_filename, 'r') as fp:
        config = json.load(fp)

    # EDIT
    # Change this to be the config update you actually want, if any
    # config['parameters']["Base_Population_Scale_Factor"] = 0.001  # ck4, testing only

    config['parameters']["Report_HIV_ByAgeAndGender_Add_Relationships"] = 1
    config['parameters']["Report_Transmission"] = 1
    config['parameters']["Report_HIV_ART"] = 1
    config['parameters']["Report_HIV_Mortality"] = 1
    config['parameters']["Report_Relationship_Start"] = 1

    with open(config_filename, 'w') as fp:
        json.dump(config, fp)
    return None


def modify_campaign_json(campaign_filename):
    with open(campaign_filename, 'r') as fp:
        campaign = json.load(fp)

    # EDIT
    # Change this to be the config update you actually want, if any

    # Reference Tracker ART Coverage

    ## FEMALE
    # Ages 45-100
    campaign['Events'][-4]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-4]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(1.0)

    # Ages 35-45
    campaign['Events'][-5]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-5]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.97)

    # Ages 25-35
    campaign['Events'][-6]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-6]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.894)

    # Ages 15-25
    campaign['Events'][-7]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-7]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.794)

    ## MALE
    # Ages 45-100
    campaign['Events'][-8]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-8]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.974)

    # Ages 35-45
    campaign['Events'][-9]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-9]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.902)

    # Ages 25-35
    campaign['Events'][-10]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-10]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.65)

    # Ages 15-25
    campaign['Events'][-11]['Event_Coordinator_Config']['Time_Value_Map']['Times'].append(2021)
    campaign['Events'][-11]['Event_Coordinator_Config']['Time_Value_Map']['Values'].append(0.832)

    with open(campaign_filename, 'w') as fp:
        json.dump(campaign, fp)
    return None


# END EDIT #######################################################################


def copy_simulation(simulation, to_experiment):
    simulation.refresh(query_criteria=QueryCriteria().select_children(['files', 'hpc_jobs', 'tags']))

    new_simulation = Simulation(simulation.name, description=simulation.description)
    new_simulation.experiment_id = to_experiment.id
    tags = copy.copy(simulation.tags)
    tags["CopiedFromSimulation"] = simulation.id
    tags["parameterization_id"] = tags["TPI"]
    tags["Run_Number"] = tags.get("Run_Number", 1)  # dummy run number if original sim does not have it as a tag
    new_simulation.set_tags(tags)

    job = simulation.hpc_jobs[-1]

    # override any fields here as necessary...
    if job and job.configuration:
        new_simulation.configuration = Configuration(
            environment_name=job.configuration.environment_name,
            simulation_input_args="--config config.json --input-path ./Assets",  # job.configuration.simulation_input_args,
            working_directory_root=job.configuration.working_directory_root,
            executable_path=job.configuration.executable_path,
            maximum_number_of_retries=SetupParser.get(parameter='num_retries'),
            priority=SetupParser.get(parameter='priority'),
            min_cores=job.configuration.min_cores,
            max_cores=job.configuration.max_cores,
            exclusive=job.configuration.exclusive,
            node_group_name=SetupParser.get(parameter='node_group'),
            asset_collection_id=job.configuration.asset_collection_id)

    with tempfile.TemporaryDirectory() as dir:
        files_to_add_last = {}
        for f in simulation.files:
            if f.file_name == 'config.json':
                # modify the config file from original simulation as desired/needed
                dest_file = os.path.join(dir, 'config.json')
                with open(dest_file, 'wb') as fp:
                    fp.write(f.retrieve())
                modify_config_json(config_filename=dest_file)
                filename = dest_file
                sf = SimulationFile(file_name=filename, file_type=f.file_type, description=f.description)
                files_to_add_last[filename] = sf
            elif campaign_file_regex.match(f.file_name) is not None:
                # modify the campaign file from original simulation as desired/needed
                dest_file = os.path.join(dir, f.file_name)
                with open(dest_file, 'wb') as fp:
                    fp.write(f.retrieve())
                modify_campaign_json(campaign_filename=dest_file)
                filename = dest_file
                sf = SimulationFile(file_name=filename, file_type=f.file_type, description=f.description)
                files_to_add_last[filename] = sf
            else:
                # use these files as-is from original simulation
                filename = f.file_name
                checksum = f.md5_checksum
                sf = SimulationFile(file_name=filename, file_type=f.file_type, description=f.description,
                                    md5_checksum=checksum)
                new_simulation.add_file(sf)

        new_simulation.save(return_missing_files=False)
        if len(files_to_add_last) > 0:
            for file_path, sf in files_to_add_last.items():
                new_simulation.add_file(sf, file_path=file_path)
            new_simulation.save(return_missing_files=False)

    print('new sim = ' + str(new_simulation.id))

    return new_simulation


def copy_experiment(experiment, to_suite):
    new_experiment = Experiment(name=experiment.name, suite_id=to_suite.id)
    new_experiment.set_tags({"CopiedFromExperiment": experiment.id})
    new_experiment.save()
    return new_experiment


if __name__ == "__main__":
    from simtools.Utilities.COMPSUtilities import exps_for_suite_id, sims_from_experiment_id
    from simtools.Utilities.Experiments import retrieve_experiment
    from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
    from simtools.ExperimentManager.BaseExperimentManager import BaseExperimentManager

    SetupParser.init()

    # load the original stuff
    original_suite = Suite.get(id=original_suite_id)
    original_experiments = exps_for_suite_id(original_suite_id)

    # start making the new stuff
    exp_manager = ExperimentManagerFactory.init()
    suite_id = exp_manager.create_suite(suite_name=original_suite.name)
    new_suite = Suite.get(id=suite_id)

    new_experiments = []
    for original_experiment in original_experiments:
        print(f"considering exp id: {original_experiment.id} ...")
        if str(original_experiment.id) in EXPERIMENT_IDS_TO_RERUN:
            new_experiment = copy_experiment(experiment=original_experiment, to_suite=new_suite)
            new_experiments.append(new_experiment)

            # simulation level items to set
            original_simulations = sims_from_experiment_id(exp_id=original_experiment.id)
            n_simulations = len(original_simulations)
            for sim_count, original_simulation in enumerate(original_simulations):
                new_simulation = copy_simulation(simulation=original_simulation, to_experiment=new_experiment)
                print(f"{sim_count+1}/{n_simulations} simulations created...")
