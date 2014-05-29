#!/usr/bin/env python

"""A module for building up the project objects that build up the project database on 
statusdb with lims as the main source of information.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""
import load_status_from_google_docs  ### Temorary solution untill 20158 implemented in LIMS!!!
import codecs
from scilifelab.google import _to_unicode, _from_unicode
from pprint import pprint
from genologics.lims import *
from helpers import *
from lims_utils import *
from scilifelab.db.statusDB_utils import *
import os
import couchdb
import bcbio.pipeline.config_utils as cl
import time
from datetime import date


config_file = os.path.join(os.environ['HOME'], 'opt/config/post_process.yaml')
db_conf = cl.load_config(config_file)['statusdb']
url = db_conf['username']+':'+db_conf['password']+'@'+db_conf['url']+':'+str(db_conf['port'])
samp_db = couchdb.Server("http://" + url)['samples']

class ProjectDB():
    """Instances of this class holds a dictionary formatted for building up the project database on statusdb. 
    Source of information come from different lims artifacts and processes. A detailed documentation of the 
    source of all values is found in: 
    https://docs.google.com/a/scilifelab.se/document/d/1OHRsSI9btaBU4Hb1TiqJ5wwdRqUQ4BAyjJR-Nn5qGHg/edit#"""
    def __init__(self, lims_instance, project_id):
        self.lims = lims_instance 
        self.lims_project = Project(self.lims,id = project_id)
        self.preps = ProcessInfo(self.lims , self.lims.get_processes(projectname = self.lims_project.name, type = AGRLIBVAL.values()))
        runs = self.lims.get_processes(projectname = self.lims_project.name, type = SEQUENCING.values())
        self.runs = ProcessInfo(self.lims, runs)
        project_summary = self.lims.get_processes(projectname = self.lims_project.name, type = SUMMARY.values())
        self.project = {'source' : 'lims',
            'samples':{},
            'open_date' : self.lims_project.open_date,
            'close_date' : self.lims_project.close_date,
            'entity_type' : 'project_summary',
            'contact' : self.lims_project.researcher.email,
            'project_name' : self.lims_project.name,
            'project_id' : self.lims_project.id}
        self.project = get_udfs('details', self.project, self.lims_project.udf.items(), PROJ_UDF_EXCEPTIONS)
        if dict(self.lims_project.researcher.lab.udf.items()).has_key('Affiliation'):
            self.project['affiliation'] = dict(self.lims_project.researcher.lab.udf.items())['Affiliation']
        if len(project_summary) == 1:
            self.project = get_udfs('project_summary', self.project, project_summary[0].udf.items())
        elif len(project_summary) > 1:
            print 'Warning. project summary process run more than once'

        #################Temporary solution untill 20158 implemented in lims >>>>>>>>>>>>>>>>>>>>>>>
        ## can be reooved when all project opened before 2014-01-27 have been closed for more than 60 days
        ## then we also need to block old projects so that they are not overwriten in case of someone manualy 
        ## updating it with the -p flagg
        opened = self.lims_project.open_date
        if opened:
            googledocs_status = {}
            #if comp_dates(opened, '2014-01-27'):
            try:
                googledocs_status = load_status_from_google_docs.get(self.lims_project.name)
            except:
                print 'issues finding status from 20158'
                pass
        seq_finished = None
        if self.lims_project.close_date and len(runs) > 0:
            d = '2000-10-10'
            for run in runs:
                try:
                    new_date = dict(run.udf.items())['Finish Date'].isoformat()
                    if comp_dates(d,new_date):
                        d = new_date
                    seq_finished = d
                except:
                    pass
        self.project['sequencing_finished'] = seq_finished
        #Temporary solution untill 20158 implemented in lims <<<<<<<<<<<<<<<<<<<<<<<


        samples = self.lims.get_samples(projectlimsid = self.lims_project.id)
        self.project['no_of_samples'] = len(samples)
        if len(samples) > 0:
            self.project['first_initial_qc'] = '3000-10-10'
            for samp in samples:
                sampDB = SampleDB(self.lims,
                                samp.id,
                                self.project['project_name'],
                                self.project['application'],
                                self.preps.info,
                                self.runs.info,
                                googledocs_status) #googledocs_status Temporary solution untill 20158 implemented in lims!!
                self.project['samples'][sampDB.name] = sampDB.obj
##### initial qc fixa
                try:
                    initial_qc_start_date = self.project['samples'][sampDB.name]['initial_qc']['start_date']
                    if comp_dates(initial_qc_start_date,self.project['first_initial_qc']):
                        self.project['first_initial_qc'] = initial_qc_start_date
                except:
                    pass
        self.project = delete_Nones(self.project)


class ProcessInfo():
    """This class takes a list of process type names. Eg 'Aggregate QC (Library Validation) 4.0'
    and forms  a dict with info about all processes of the type specified in runs which the 
    project has gon through.

    info = {24-8460:{'finish_date':'2013-04-20', 
              'start_date',
              'run_id':'24-8460',
              'samples':{'P424_111':{in_art_id1 : [in_art1, out_art1],
                         in_art_id2: [in_art2, out_art2]},
                     'P424_115': ...},
                       ...},
        '24-8480':...}"""

    def __init__(self, lims_instance, runs):
        self.lims = lims_instance
        self.info = self.get_run_info(runs)

    def get_run_info(self, runs):
        run_info = {}
        for run in runs:
            run_info[run.id] = {'type' : run.type.name ,'start_date': run.date_run,'samples' : {}}
            run_udfs = dict(run.udf.items())
            try:
                run_info[run.id]['run_id'] = run_udfs["Run ID"]
            except:
                pass
            try:
                run_info[run.id]['finish_date'] = run_udfs['Finish Date'].isoformat()
            except:
                run_info[run.id]['finish_date'] = None
                pass
            in_arts=[]
            for IOM in run.input_output_maps:
                in_art_id = IOM[0]['limsid']
                in_art = Artifact(self.lims, id= in_art_id)
                out_art_id = IOM[1]['limsid']
                out_art = Artifact(self.lims, id= out_art_id)
                samples = in_art.samples
                if in_art_id not in in_arts:
                    in_arts.append(in_art_id)
                    for samp in samples:
                        if not samp.name in run_info[run.id]['samples'].keys():
                            run_info[run.id]['samples'][samp.name] = {}
                        run_info[run.id]['samples'][samp.name][in_art_id] = [in_art, out_art]
        return run_info



class SampleDB():
    """
    Instances of this class holds a dictionary formatted for building up the samples in the project 
    database on status db. Source of information come from different lims artifacts and processes. 
    A detailed documentation of the source of all values is found in
    https://docs.google.com/a/scilifelab.se/document/d/1OHRsSI9btaBU4Hb1TiqJ5wwdRqUQ4BAyjJR-Nn5qGHg/edit#"""
    def __init__(self,lims_instance , sample_id, project_name, application = None, prep_info = [], run_info = [], googledocs_status = {}): # googledocs_status temporary solution untill 20158 implemented in lims!!
        self.lims = lims_instance
        self.lims_sample = Sample(self.lims, id = sample_id)
        self.name = self.lims_sample.name
        self.application = application
        self.outin, self.inout = make_sample_artifact_maps(self.name)
        self.obj = {'scilife_name' : self.name}
        self.obj = get_udfs('details', self.obj, self.lims_sample.udf.items(), SAMP_UDF_EXCEPTIONS)
        preps = self.get_initQC_preps_and_libval(prep_info)
        self.obj['first_prep_start_date'] = self.get_firts_day(self.name, PREPSTART.values() + PREPREPSTART.values())
        if self.application == 'Finished library':
            self.obj['first_initial_qc_start_date'] = self.get_firts_day(self.name, INITALQCFINISHEDLIB.values())
        else:
            self.obj['first_initial_qc_start_date'] = self.get_firts_day(self.name, INITALQC.values())
        if preps:
            runs = self.get_sample_run_metrics(run_info, preps)
            if preps.has_key('library_prep'):
                for prep in runs.keys():
                    if preps['library_prep'].has_key(prep):
                        preps['library_prep'][prep]['sample_run_metrics'] = runs[prep]
                self.obj['library_prep'] = self.get_prep_leter(preps['library_prep'])
            if preps.has_key('initial_qc'):
                self.obj['initial_qc'] = preps['initial_qc']
        try:
            if (googledocs_status is None)|(googledocs_status == {}):
                self.obj['status'] = 'doc_not_found'
                self.obj['m_reads_sequenced'] = 'doc_not_found'        
            else:
                # Temporary solution untill 20158 implemented in lims!!
                self.obj['status'] = googledocs_status[self.name][0]
                self.obj['m_reads_sequenced'] = googledocs_status[self.name][1]
        except:
            pass
        self.obj = delete_Nones(self.obj)

    def get_firts_day(self, sample_name ,process_list):
        """process_list is a list of process type names, 
        sample_name is a sample name :)"""
        arts = self.lims.get_artifacts(sample_name = sample_name, process_type = process_list)
        day = date.today().isoformat()
        for a in arts:
            new_day = a.parent_process.date_run
            if comp_dates(new_day, day):
                day = new_day
        if day==date.today().isoformat():
            day = None
        return day
 
        sample_runs['library_prep'] = delete_Nones(library_prep)
        return delete_Nones(sample_runs)

    def get_initQC_preps_and_libval(self, AgrLibQCs):
        """Input: AgrLibQCs - instance of the ProcessInfo class with AGRLIBVAL processes as argument.
        For each AGRLIBVAL process run on the sample, this function steps bacward in the artifact history of the 
        output artifact of the AGRLIBVAL process to find the folowing information:

        prep_status                     The qc_flag of the input artifact of process type AGRLIBVAL
        prep_start_date                 The date-run of the PREPSTART step 
        prep_finished_date              The date-run of a PREPEND step.
        pre_prep_start_date             The date-run of process 'Shear DNA (SS XT) 4.0'. Only for 
                                        'Exome capture' projects
                          
        Preps are  defined by the date of any PREPSTART step"""
        initial_qc = {}
        library_prep ={}
        top_level_agrlibval_steps = self.get_top_level_agrlibval_steps(AgrLibQCs)
        for AgrLibQC_id in top_level_agrlibval_steps.keys():
            AgrLibQC_info = AgrLibQCs[AgrLibQC_id]
            if AgrLibQC_info['samples'].has_key(self.name):
                inart, outart = AgrLibQC_info['samples'][self.name].items()[0][1]
                history = get_analyte_hist(outart.id, self.outin, self.inout)
                initial_qc = self.get_initial_qc_dates(history)
                prep_id, prep_info = self.get_lib_prep(history, inart)
                if not library_prep.has_key(prep_id):
                    library_prep[prep_id] = prep_info
                library_prep[prep_id]['library_validation'][AgrLibQC_id] = self.get_libval(history, AgrLibQC_info, inart)
        return {'library_prep':library_prep,'initial_qc':initial_qc}

    def get_lib_prep(self, history, inart):
        if self.application == 'Finished library':
            return 'Finished', {'prep_status':inart.qc_flag,'reagent_labels':self.lims_sample.artifact.reagent_labels,  'library_validation':{}}
        prep_info = {'prep_status' : inart.qc_flag, 'reagent_labels' : inart.reagent_labels, 'library_validation':{}}
        libPrep_id = None
        for step, info in history.items():
            if info['type'] in PREPREPSTART.keys():
                libPrep_id = info['id']
                prep_info['pre_prep_start_date'] = info['date']
            elif info['type'] in PREPSTART.keys():
                if self.application !='Exome capture':
                    libPrep_id = info['id']
                prep_info['prep_start_date'] = info['date']
            elif info['type'] in PREPEND.keys():
                prep_info['prep_finished_date'] = info['date']
                prep_info['prep_id'] = info['id']
            elif info['type'] in WORKSET.keys():
                prep_info['workset_setup'] = info['id']
        return libPrep_id, prep_info 

    def get_libval(self, history, AgrLibQC_info, inart):
        """library_validation/start_date   First of all LIBVAL steps found for in the artifact history 
                                        of the output artifact of one of the AGRLIBVAL step 
        library_validation/finish_date  date-run of AGRLIBVAL step 
        average_size_bp                 udf ('Size (bp)') of the input artifact to the process AGRLIBVAL"""
        library_validation = {'start_date' : self.get_lib_val_start_dates(history),
                         'finish_date' : AgrLibQC_info['start_date']}
        for key, val in inart.udf.items():
            key = key.replace(' ', '_').lower().replace('.','')
            if key=="size_(bp)":        ##remove when lims is updated with new key name....
                key="average_size_bp"
            library_validation[key] = val   
        return delete_Nones(library_validation)

    def get_top_level_agrlibval_steps(self, AgrLibQCs):
        topLevel_AgrLibQC={}
        for AgrLibQC_id, AgrLibQC_info in AgrLibQCs.items():
            if AgrLibQC_info['samples'].has_key(self.name):
                topLevel_AgrLibQC[AgrLibQC_id]=[]
                inart, outart = AgrLibQC_info['samples'][self.name].items()[0][1]
                history = get_analyte_hist(outart.id, self.outin, self.inout)
                for step, info in history.items():
                    if info['type'] in LIBVAL.keys():
                        topLevel_AgrLibQC[AgrLibQC_id].append(step)
        for AgrLibQC, LibQC in topLevel_AgrLibQC.items():
            LibQC=set(LibQC)
            for AgrLibQC_comp, LibQC_comp in topLevel_AgrLibQC.items():
                LibQC_comp=set(LibQC_comp)
                if AgrLibQC_comp != AgrLibQC:
                    if LibQC.issubset(LibQC_comp) and topLevel_AgrLibQC.has_key(AgrLibQC):
                        topLevel_AgrLibQC.pop(AgrLibQC) 
        return topLevel_AgrLibQC

    def get_prep_leter(self, prep_info):
        """Get preps and prep names; A,B,C... based on prep dates for sample_name. 
        Output: A dict where keys are prep_art_id and values are prep names."""
        dates = {}
        prep_info_new = {}
        preps_keys = map(chr, range(65, 65+len(prep_info)))
        if len(prep_info) == 1:
            prep_info_new['A'] = prep_info.values()[0]
        else:
            for key, val in prep_info.items():
                dates[key] = val['prep_start_date']
            for i, key in enumerate(sorted(dates,key= lambda x : dates[x])):
                prep_info_new[preps_keys[i]] = prep_info[key]
        return prep_info_new

    def get_sample_run_metrics(self, SeqRun_info, preps):
        """Input: SeqRun_info - instance of the ProcessInfo class with SEQUENCING processes as argument
        For each SEQUENCING process run on the sample, this function steps bacward in the artifact history of the 
        input artifact of the SEQUENCING process to find the folowing information:

        dillution_and_pooling_start_date    date-run of DILSTART step
        sequencing_start_date               date-run of SEQSTART step
        sequencing_run_QC_finished          date-run of SEQUENCING step
        sequencing_finish_date              udf ('Finish Date') of SEQUENCING step
        sample_run_metrics_id               The sample database (statusdb) _id for the sample_run_metrics 
                                            corresponding to the run, sample, lane in question.
        samp_run_met_id = lane_date_fcid_barcode            
            date and fcid: from udf ('Run ID') of the SEQUENCING step. 
            barcode: The reagent-lables of the input artifact of process type AGRLIBVAL
            lane: from the location of the input artifact to the SEQUENCING step    
        preps are defined as the id of the PREPSTART step in the artifact history. If appllication== Finished library, 
        prep is defined as "Finnished". These keys are used to connect the seqeuncing steps to the correct preps."""
        sample_runs = {}
        for id, run in SeqRun_info.items():
            if run['samples'].has_key(self.name) and run.has_key('run_id'):
                date = run['run_id'].split('_')[0]
                fcid = run['run_id'].split('_')[3]
                run_type = run['type']
                for id , arts in run['samples'][self.name].items():
                    lane_art = arts[0]
                    outart = arts[1]
                    if run_type == "MiSeq Run (MiSeq) 4.0":
                        lane = lane_art.location[1].split(':')[1]
                    else:
                        lane = lane_art.location[1].split(':')[0]
                    history = get_analyte_hist(lane_art.id, self.outin, self.inout)
                    key = None
                    dillution_and_pooling_start_date = None
                    for step , info in history.items():
                        if info['type'] in PREPSTART.keys():
                            key = info['id']
                        elif info['type'] in DILSTART.keys():
                            dillution_and_pooling_start_date = info['date']
                            type = info['type']
                        elif info['type'] in SEQSTART.keys():
                            sequencing_start_date = info['date']
                        if self.application == 'Finished library' :
                            key = 'Finished'
                    if key:
                        try:
                            barcode = self.get_barcode(preps['library_prep'][key]['reagent_labels'])
                            samp_run_met_id = '_'.join([lane, date, fcid, barcode])
                        except:
                            samp_run_met_id = None
                        dict = {'dillution_and_pooling_start_date': dillution_and_pooling_start_date,
                                'sequencing_start_date':sequencing_start_date,
                                'sequencing_run_QC_finished': run['start_date'],
                                'sequencing_finish_date': run['finish_date'],
                                'sample_run_metrics_id': find_sample_run_id_from_view(samp_db, samp_run_met_id) }
                        dict = delete_Nones(dict)
                        if not sample_runs.has_key(key): 
                            sample_runs[key] = {}
                        sample_runs[key][samp_run_met_id] = dict
        return sample_runs

    def get_sample_status():
        """ongoing,passed,aborted"""
        ##    Not yet implemented

    def get_barcode(self, reagent_lables):
        """Extracts barcode from list of artifact.reagent_labels"""
        if len(reagent_lables)>1:
            return None
        else: 
            try:
                index =reagent_lables[0].split('(')[1].strip(')')
            except:
                index = reagent_lables[0]
        return index
        

    def get_initial_qc_dates(self, history):
        """Extracts run dates for processes of type AGRINITQC 
        from a history dict."""
        initial_qc_finish_date = None
        for step , info in history.items():
            if (info['type'] in AGRINITQC.keys()) and info['date']:
                if initial_qc_finish_date is None:
                    initial_qc_finish_date = info['date']
                elif comp_dates(initial_qc_finish_date, info['date']):
                    initial_qc_finish_date = info['date']
        if initial_qc_finish_date:
            initial_qc_start_date = initial_qc_finish_date
            for step , info in history.items():
                if (info['type'] in INITALQC) and info['date']:
                    if comp_dates(info['date'], initial_qc_start_date):
                        initial_qc_start_date = info['date']
            return {'start_date' : initial_qc_start_date, 'finish_date' : initial_qc_finish_date}
        else: 
            return

    def get_lib_val_start_dates(self, history):
        """Extracts run dates for processes of type LIBVAL 
        from a history dict."""
        lib_val_start_date = None
        for step , info in history.items():
            if info['type'] in LIBVAL.keys():
                if lib_val_start_date == None:
                    lib_val_start_date = info['date']
                elif comp_dates(info['date'], lib_val_start_date):
                    lib_val_start_date = info['date']
        return lib_val_start_date
