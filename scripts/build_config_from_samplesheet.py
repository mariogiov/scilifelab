#!/bin/env python
"""
Build a bcbio-nextgen run config using info from the SampleSheet.csv and
project info from the StatusDB.
"""
from __future__ import print_function

import argparse
import collections
import datetime
import functools
import os
import re
import shutil
import sys

from bcbio.solexa import samplesheet
from bcbio.workflow import template
from scilifelab.db import statusDB_utils as sdbu

## PLAN whether a sample is DNA or RNA controls which alignment Pipeline is executed
##      (i.e. RNASeq w/Bowtie2 or Standard DNA w/bwa)
## TODO this could be refined further to ChIPSeq, RADSeq, miRNA, etc.
SAMPLE_TYPE = { "genomic dna"       :   "dna",
				"amplified gdna"    :   "dna",
				"other type of dna" :   "dna",
				"total rna"         :   "rna",
				"mrna"              :   "rna",
				"small rna"         :   "rna",}
                ## problem values
				#"Finished Library",
                # TODO does this mean DNA?
				#"Amplicon with adapters",
				#"Amplicon without adapters",

## PLAN Anything in this list will hit a mapping pipeline; anything else will not.
SUPPORTED_GENOMES = set(["hg19",
                         "mm9",
                         "rn4",
                         "saccer2",
                         "dm3",
                         "tair9",
                         "xentro2",
                         "ws210",
                         "canfam3",])
                         ## problem values
                         #"multiple",
                         #"other",

### TEMPORARY ###
## TODO paths for testing, remove later as they will be loaded from a config file
qc_pipeline_template        = "/pica/h1/mario/bcbio-nextgen/config/templates/gatk-variant.yaml"
bwa_alignment_template      = "/pica/h1/mario/bcbio-nextgen/config/templates/gatk-variant.yaml"
rnaseq_alignment_template   = "/pica/h1/mario/bcbio-nextgen/config/templates/illumina-rnaseq.yaml"
### TEMPORARY ###

## TODO these template variables will have as values paths that are defined in the system config
PIPELINE_TEMPLATES  = { "dna"           : (qc_pipeline_template, bwa_alignment_template),
                        "rna"           : (qc_pipeline_template, rnaseq_alignment_template),
                        "unsupported"   : (qc_pipeline_template),}


def build_config_file_for_project(input_files, status_db_config=None, output_dir=None, upload_dir=None):
    ## TODO check the file naming convention, not sure on this one
    """
    Builds a bcbio-nextgen workflow template for a sample (or set of R1/R2 samples)
    using project information from the StatusDB.
    Project ID is pulled from the standard filename format, which is:
       <lane_num>_<date>_<fcid>_<project>_<sample_num>_<read>.fastq[.gz]
    The function needs information about how to access the StatusDB from status_db_config
    Writes the resulting configuration file to
        <sample_directory>/<project_name>/config/project_config.yaml
    """
    if not input_files:
        raise SyntaxError("Input files must be specified.")
    couch = sdbu.load_couch_server(status_db_config)
    if couch:
        proj_db = couch['projects']
    else:
        raise RuntimeError("Couldn't connect to StatusDB or "\
                           "config file lacked authentication information.")
    file_pairs = find_fastq_read_pairs(input_files)
    for sample_basename, sample_files in file_pairs.items():
        project_id = get_project_id_from_sample_id(sample_basename)
        if not output_dir:
            output_dir = os.path.join(os.path.dirname(os.path.abspath(sample_files[0])), "project_{}".format(sample_basename))
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        if not upload_dir:
            upload_dir = os.path.join(output_dir, "final")
        project_data_json = get_project_data_for_id(project_id, proj_db)
        ## TODO This has yet to be implemented in StatusDB so always is None,
        ##      which means the defaults in the templates will stay
        adapter_seqs      = project_data_json.get("adapter_sequences")
        reference_genome  = project_data_json.get("reference_genome")
        if reference_genome.lower() not in SUPPORTED_GENOMES:
            reference_genome = None
        sample_type       = SAMPLE_TYPE.get(project_data_json.get("details", None).get("sample_type", None).lower(), "unsupported")
        ## TODO need some way to merge the two template files
        ##      or have template.py modify both pipelines when it adds attributes.
        ##      Optionally we could just have multiple config files -- each pipeline runs separately;
        ##      the issue is then just giving each one a distinct filename.
        pipeline_list     = PIPELINE_TEMPLATES.get(sample_type, (qc_pipeline_template))
        for template_name in pipeline_list:
            namespace_args = argparse.Namespace(template=template_name,
                                                input_files=sample_files,
                                                out_dir=output_dir,
                                                upload_dir=upload_dir)
            if adapter_seqs or reference_genome:
                # Create the csv file that will be passed to the template creator
                project_csv = create_project_csv_from_dbinfo(sample_basename, output_dir,
                                                             adapter_seqs, reference_genome)
                namespace_args.__dict__["metadata"] = project_csv
            ## TODO We need to know the paths to the config file(s)
            ##      so we can kick off processing via the webserver.
            ##      Probably best to get this as a return value from template()
            template.setup(namespace_args)


def find_fastq_read_pairs(file_list):
    """
    Given a list of file names, finds read pairs (based on _R1_/_R2_ file naming)
    and returns a dict of {base_name: [ file_read_one, file_read_two ]}
    E.g.
        1_131129_BH7VPTADXX_P602_101_1.fastq.gz
        1_131129_BH7VPTADXX_P602_101_2.fastq.gz
    becomes
        { "1_131129_BH7VPTADXX_P602_101":
        [ "1_131129_BH7VPTADXX_P602_101_1.fastq.gz",
          "1_131129_BH7VPTADXX_P602_101_2.fastq.gz"]}
    """
    # Remove duplicates
    file_set = set(file_list)
    # Split on the read number
    split_pattern = re.compile(r'_\d\.fastq')
    matches_dict = collections.defaultdict(list)
    for file_name in file_list:
        file_basename = os.path.basename(file_name)
        try:
            base = split_pattern.split(file_basename)[0]
            matches_dict[base].append(file_name)
        except IndexError:
            print("Warning: file doesn't match expected file format, "
                  "cannot be paired: \"{}\"".format(file_name), file=sys.stderr)
            matches_dict[file_basename].append(file_name)
    return dict(matches_dict)


def create_project_csv_from_dbinfo(file_basename, output_dir, adapter_seqs, reference_genome):
    """
    Creates a file to be used by the template generator to modify the existing
    template. Output looks like:
        samplename,description,adapters,trim_reads,genome_build
        7_100326_FC6107FAAXX,7_100326_FC6107FAAXX,"[truseq,polya]",read_through,hg19
    Returns the path to the file that is created.
    """
    project_csv_file = os.path.join(output_dir, "{}.csv".format(file_basename))
    if os.path.exists(project_csv_file):
        shutil.move(project_csv_file,
                    project_csv_file + ".bak%s" % datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))
    with open(project_csv_file, 'w') as f:
        header_line = ["samplename","description"]
        sample_line = [file_basename,file_basename]
        if adapter_seqs:
            ## TODO Figure out how to deal with adapter sequences -- both fwd and revcom?
            ##      What about mixing default values (e.g. truseq, polya) and custom seqs?
            ##          There is some issue with this in the bcbio code I think
            ##      Maybe it's best just to use custom sequences in case they change upstream
            adapter_seqs = "\"[{}]\"".format(",".join(adapter_seqs))
            header_line.extend(["adapters","trim_reads"])
            sample_line.extend([adapter_seqs,"read_through"])
        if reference_genome:
            header_line.append("genome_build",)
            sample_line.append(reference_genome)
        f.write(",".join(header_line)+"\n")
        f.write(",".join(sample_line)+"\n")
    return project_csv_file


class memoized(object):
    """
    Decorator, caches results of function calls.
    """
    def __init__(self, func):
        self.func   = func
        self.cached = {}
    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            return self.func(*args)
        if args in self.cached:
            return self.cached[args]
        else:
            return_val = self.func(*args)
            self.cached[args] = return_val
            return return_val
    def __repr__(self):
        return self.func.__doc__
    # This ensures that attribute access (e.g. obj.attr)
    # goes through the __call__ function I defined above
    # functools is awesome
    # descriptors are the raddest
    # boy i love python
    def __get__(self, obj, objtype):
        return functools.partial(self.__call__, obj)


@memoized
def get_project_id_from_sample_id(sample_basename):
    """
    Project is pulled from the standard filename format, which is:
       <lane_num>_<date>_<fcid>_<project>_<sample_num>_<read>.fastq[.gz]
    returns the project portion or None if there is no match
    (which shouldn't generally happen).
    """
    try:
        project_id = re.match(r'\d_\d{6}_\w{10}_(P\d{3})_.*', sample_basename).groups()[0]
        return project_id
    except (IndexError, AttributeError):
        raise ValueError("Error: filename didn't match conventions, "
                         "couldn't find project id for sample "
                         "\"{}\"".format(sample_basename))


@memoized
def get_project_data_for_id(project_id, proj_db):
    """
    Pulls all the data about a project from the StatusDB
    given the project's id (e.g. "P602") and a couchdb view object.
    Returns a JSON string of the data.
    """
    db_view = proj_db.view('project/project_id')
    try:
        return proj_db.get([proj.id for proj in db_view if proj.key == project_id][0])
    except IndexError:
        # TODO this will be logged and should be caught on the calling side
        raise ValueError("Warning: project ID '{}' not found in Status DB".format(project_id))


if __name__ == "__main__":
    # This entire set of arguments is a stopgap
    parser = argparse.ArgumentParser("Get project information for a sample given "\
                                     "the SampleSheet.csv and a config file "\
                                     "with StatusDB auth info.")
    # TODO this will be a constant I suppose for all runs
    parser.add_argument("-c", "--dbconfig", dest="status_db_config", required=True,
            help="The config file (post_process.yaml) containing StatusDB info.")
    # TODO replace this with information based on where the samples are already (same dir)
    parser.add_argument("-o", "--output-dir",
            help="The directory to write the project config/work files to.")
    # TODO replace this with information based on where the samples are already
    parser.add_argument("-u", "--upload-dir",
            help="The directory to upload processed data to")
    parser.add_argument("input_files", nargs="+", metavar="input.fastq",
            help="The paths to the input files.")

    kwargs = vars(parser.parse_args())
    build_config_file_for_project(**kwargs)
