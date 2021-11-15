from __future__ import print_function
import sys
import copy
import tempfile
import os
import json
import hashlib

from COMPS import Client
from COMPS.Data import Simulation, SimulationFile, QueryCriteria, Configuration, Experiment, Suite

from simtools.SetupParser import SetupParser


#### EDIT HERE #######################################################################

SetupParser.default_block = 'HPC'

# The suite of experiments/simulations to copy, modify, and run
original_suite_id = 'c79aa54e-2578-ea11-a2c5-c4346bcb1550'


def modify_config_json(filename):
    import json
    with open(filename, 'r') as fp:
        config = json.load(fp)

    # EDIT
    # Change this to be the congig update you actually want
    config['parameters']["Report_HIV_ByAgeAndGender_Add_Relationships"] = 1
    config['parameters']["Report_Transmission"] = 1
    config['parameters']["Report_HIV_ART"] = 1
    config['parameters']["Report_HIV_Mortality"] = 1
    config['parameters']["Report_Relationship_Start"] = 1


    with open(filename, 'w') as fp:
        json.dump(config, fp)
    return None

#### END EDIT #######################################################################


def copy_simulation(simulation, to_experiment):
    simulation.refresh(query_criteria=QueryCriteria().select_children(['files', 'hpc_jobs', 'tags']))

    new_simulation = Simulation(simulation.name, description=simulation.description)
    new_simulation.experiment_id = to_experiment.id

    tags = copy.copy(simulation.tags)
    tags["CopiedFromSimulation"] = simulation.id
    new_simulation.set_tags(tags)

    job = simulation.hpc_jobs[-1]

    # override any fields here as necessary...
    if job and job.configuration:
        new_simulation.configuration = Configuration(
            environment_name=job.configuration.environment_name,
            simulation_input_args=job.configuration.simulation_input_args,
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
                dest_file = os.path.join(dir, 'config.json')
                with open(dest_file, 'wb') as fp:
                    fp.write(f.retrieve())
                modify_config_json(dest_file)
                # with open(dest_file, 'rb') as fp:
                #     data = fp.read()
                filename = dest_file
                # checksum = hashlib.md5(data).hexdigest()
                sf = SimulationFile(file_name=filename, file_type=f.file_type, description=f.description)
                files_to_add_last[filename] = sf
            else:
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
        new_experiment = copy_experiment(experiment=original_experiment, to_suite=new_suite)
        new_experiments.append(new_experiment)

        # simulation level items to set
        original_simulations = sims_from_experiment_id(exp_id=original_experiment.id)
        for original_simulation in original_simulations:
            new_simulation = copy_simulation(simulation=original_simulation, to_experiment=new_experiment)
