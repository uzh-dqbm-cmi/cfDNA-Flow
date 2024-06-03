#! /usr/bin/env python3

import argparse
import filecmp
import glob
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import cfDNA
import configure
import scripts.bedpe_size_selection as szsel

# test constants
# should be the same as in then used config.yaml (test_cfDNA_pipeline*.yaml)
WPS_WINDOW = 120
WPS_MIN = 150
WPS_MAX = 180
WPS_BIGWIG_NORMALIZED = False
WPS_OUT_SUFFIX = f"_w{WPS_WINDOW}_frg{WPS_MIN}-{WPS_MAX}" if WPS_WINDOW != "120" or WPS_MIN != "120" or WPS_MAX != "180" else ""
FFT_BINSIZE = 1000000

# import shlex
# import sys

# for singularity image tests
# INSTALLATION_FOLDER = '/'
# Default: for local tests, pytest
# INSTALLATION_FOLDER = './'


class MainArgs:
    analysis: str
    singularity_images: str
    local: bool
    config: str
    j: int = 1
    np: bool = True
    unlock: bool = False


def assure_folder(folder_name):
    if not os.path.exists(folder_name):  # or not os.path.isdir(folder_name):
        try:
            os.makedirs(folder_name)
        except:
            print(f'Folder {folder_name} was created by other thread.')


def msg_diff_content(result, expected, log=None):
    diff = ""
    if log:
        with open(log, "r") as log_f:
            diff += "\n########### Log output:\n"
            diff += log_f.read()
            diff += "\n########### End of Log output:\n"
    diff += f"Unexpected content of file {result} \n########### Resulted:\n"
    with open(result, "r") as left, open(expected, "r") as right:
        diff += left.read()
        diff += "\n########### Expected:\n"
        diff += right.read()
    return diff


# @pytest.mark.skip(reason="deprecated command: use configure.py")
def test_main(inst_folder):
    cleanup_files(['results/test_cfDNA_pipeline_test001.yaml'])
    args = MainArgs()
    args.local = False
    args.analysis = 'do_preprocess'
    args.singularity_images = f'{inst_folder}test/singularity_images'
    args.config = f'{inst_folder}test/test_config.yaml'
    # args.j = 1
    # args.np = True
    # args.unlock = False
    cfDNA.main(args, test=True)
    result = 'results/test_cfDNA_pipeline_test001.yaml'
    expected = f'{inst_folder}test/expected/test_cfDNA_pipeline_test001.yaml'
    assert filecmp.cmp(result, expected, shallow=False)


def test_configure(inst_folder):
    config = f'{inst_folder}test/test_config.yaml'
    cleanup_files(['results/test_cfDNA_pipeline_test001.yaml'])
    result = configure.merge_config_and_settings(config)
    expected = f'{inst_folder}test/expected/test_cfDNA_pipeline_test001.yaml'
    assert filecmp.cmp(result, expected, shallow=False)


def test_configure2(inst_folder):
    config = f'{inst_folder}test/test_config.yaml'
    out_config = 'results/out_test_config.yaml'
    cleanup_files([out_config])
    result = configure.merge_config_and_settings(config, out_config)
    expected = f'{inst_folder}test/expected/test_cfDNA_pipeline_test001.yaml'
    assert filecmp.cmp(result, expected, shallow=False)


def test_size_selection_sh(inst_folder):
    # input_bed = '/Users/todor/data/cfdna/AmsterdamUMC/results/BED/FalseD25630/LP0020_02.bed'
    # result = '/Users/todor/data/cfdna/AmsterdamUMC/results/BED/FalseD25630/szsel_90_150_sh/LP0020_02.bed'
    input_bed = f'{inst_folder}test/data/BED/test_sample2.bed'
    result = 'results/test_sample2_szsel_sh.bed'
    temp = 'results/temp'
    cleanup_files([result])
    assure_folder(os.path.dirname(result))
    sz_min = 60
    sz_max = 150
    chunk_size = 3
    cores = 2
    szsel_cmd = f'{inst_folder}scripts/bedpe_size_selection.sh -i {input_bed} --min {sz_min} --max {sz_max} -c {chunk_size} -p {cores} --temp {temp} > {result}'
    run_cmd(szsel_cmd)
    # compare with all in 1 process/pass
    result_p1 = f'{result}.p1.bed'
    szsel_cmd = f'{inst_folder}scripts/bedpe_size_selection.sh -i {input_bed} --min {sz_min} --max {sz_max} -o {result_p1} --temp {temp}'
    run_cmd(szsel_cmd)
    assert filecmp.cmp(result_p1, result, shallow=False)
    # return
    expected = f'{inst_folder}test/expected/test_sample2_szsel.bed'
    assert filecmp.cmp(result, expected, shallow=False)
    # once more direct in output file
    os.remove(result)
    szsel_cmd = f'{inst_folder}scripts/bedpe_size_selection.sh -i {input_bed} --min {sz_min} --max {sz_max} -c {chunk_size} -o {result} --temp {temp}'
    run_cmd(szsel_cmd)
    expected = f'{inst_folder}test/expected/test_sample2_szsel.bed'
    assert filecmp.cmp(result, expected, shallow=False)


# ##########    snakemake tests    ##########################


def get_smk_output_status(smk_output_file, status_file):
    with open(smk_output_file, 'r') as outfile:
        with open(status_file, 'w') as resultfile:
            lines = outfile.readlines()
            status = False
            for line in lines:
                status = status or ("Job stats:" in line)
                if status:
                    if not line.strip():
                        break
                    resultfile.write(line)


def cleanup_files(file_dir_wildcards):
    for file_dir_wildcard in file_dir_wildcards:
        file_dir_list = glob.glob(file_dir_wildcard)
        for file_dir in file_dir_list:
            if os.path.exists(file_dir):
                if os.path.isdir(file_dir):
                    shutil.rmtree(file_dir)
                else:
                    os.remove(file_dir)


def run_test_smk(inst_folder, smk_rule, results, expecteds=None, expected_sizes=None, check_status_log=True,
                 config_file="test/test_cfDNA_pipeline.yaml"):
    cmd_np_out = f'results/test_smk_{smk_rule}_np.log'
    cmd_out = f'results/test_smk_{smk_rule}.log'
    status_log = f'results/test_smk_{smk_rule}_status.log'
    cleanup_files([cmd_out, status_log])
    cleanup_files(results)
    # run_cmd(f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml --unlock')
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}{config_file} -j 2 {smk_rule} -np > {cmd_np_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_np_out, status_log)
    if check_status_log:
        expected = f'{inst_folder}test/expected/test_smk_{smk_rule}_status.log'
        assert filecmp.cmp(status_log, expected, shallow=False), msg_diff_content(status_log, expected, cmd_np_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}{config_file} -j 2 --verbose {smk_rule} > {cmd_out}'
    run_cmd(smk_cmd)
    for idx, result in enumerate(results):
        result_size = os.path.getsize(result)
        if expected_sizes:
            assert (result_size >= expected_sizes[idx]), f'Unexpected file size {result_size} of file {result}!\n'
        elif expecteds:
            assert filecmp.cmp(result, expecteds[idx], shallow=False), msg_diff_content(result, expecteds[idx], cmd_out)
        else:
            print(f'Nothing to check so all good the result file {result} exists!')


def unlock(inst_folder, config_file="test/test_cfDNA_pipeline.yaml"):
    run_cmd(f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml --unlock')


def run_cmd(cmd):
    print(f'run_cmd: {cmd}')
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print("stdout:", result.stdout)
    print("stderr:", result.stderr)
    # os.system(cmd)
    # sleep(1)  # give a chance of the process to finish; normally should not be needed
    # other ways
    # print(os.popen(cmd).read())
    # cmd_args = shlex.split(cmd)
    # p = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) # somehow not working
    # p.wait()
    # print(p.returncode)
    # status, jobout = subprocess.getstatusoutput(cmd)
    # print(f" call: {cmd}\n\t {jobout}")
    # if status != 0:
    #    print(f"Error running: {cmd} \n\t status: {status}")


def sub_in_file(file, pattern, repl, bak=None):
    if bak:
        shutil.copyfile(file, file+bak)
    with open(file, 'r+') as f:
        text = f.read()
        text = re.sub(pattern, repl, text)
        f.seek(0)
        f.write(text)
        f.truncate()


@pytest.mark.skip(reason="manual local test")
def test_smk_do_bam_GCbias(inst_folder):
    smk_rule = 'do_bai'
    results = [f'results/BAM/TrueD25630/TEST_sample1.sortByCoord.bam',
               f'results/BAM/TrueD25630/TEST_sample1.sortByCoord.bai',
               f'results/BAM/TrueD25630/TEST_sample3.sortByCoord.bam',
               f'results/BAM/TrueD25630/TEST_sample3.sortByCoord.bai'
               ]
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[1, 3, 2, 4], check_status_log=False,
                 config_file="test/test_cfDNA_pipeline_GCbias.yaml")

@pytest.mark.skip(reason="manual local test")
def test_smk_do_bam_splitGCbias(inst_folder):
    smk_rule = 'do_bam'
    results = [f'results/BAM/TrueD25630/TEST_sample1.sortByCoord.bam',
               f'results/BAM/TrueD25630/TEST_sample1.sortByCoord.bai',
               f'results/BAM/TrueD25630/TEST_sample3.sortByCoord.bam',
               f'results/BAM/TrueD25630/TEST_sample3.sortByCoord.bai'
               ]
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[1, 3, 2, 4], check_status_log=False,
                 config_file="test/test_cfDNA_pipeline_GCbias_split.yaml")

# @pytest.fixture(scope='session')
@pytest.mark.order(1)
def test_smk_do_preprocess(inst_folder):
    # by problems <ou might need to run this command from the command line to unlock snakemake:
    # snakemake -s ./Snakefile --configfile ./test/test_cfDNA_pipeline.yaml --unlock
    cmd_out = 'results/test_smk_do_preprocess.log'
    result = 'results/test_smk_do_preprocess_status.log'
    result2 = 'results/BED/FalseD25630/TEST_sample1.bed'
    result3 = 'results/BED/FalseD25630/TEST_sample3.bed'
    result4 = 'results/QC/FalseD25630/multiqc_data/multiqc_report.html'
    cleanup_files([cmd_out, result, result2, result3, 'results/BAM', 'results/BED', 'results/logs',
                   'results/trimmed_reads', 'results/BED/FalseD25630/TEST_sample1.bed',
                   'results/BED/FalseD25630/TEST_sample3.bed', 'results/QC/FalseD25630/multiqc_data/multiqc_report.html'])
    assure_folder(os.path.dirname(result))

    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_preprocess -np > {cmd_out}'
    #                 f' | grep -wv "Multiple includes"'   # if needed to exclude lines
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    expected = f'{inst_folder}test/expected/test_smk_do_preprocess_status.log'
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_preprocess --verbose > {cmd_out}'
    run_cmd(smk_cmd)
    assert os.path.exists(result2)
    assert os.path.exists(result3)
    assert os.path.exists(result4)


# @pytest.fixture(scope='session')
# @pytest.mark.usefixtures('test_smk_do_preprocess')
@pytest.mark.order(2)
def test_smk_do_bedprocess(inst_folder):
    cmd_out = 'results/test_smk_do_bedprocess.log'
    result = 'results/test_smk_do_bedprocess_status.log'
    expected = f'{inst_folder}test/expected/test_smk_do_bedprocess_status.log'
    result2 = 'results/BED/FalseD25630/mergeddf.csv'
    result_lens = 'results/BED/FalseD25630/*_len.csv'
    expected2 = f'{inst_folder}test/expected/mergeddf.csv'
    cleanup_files([cmd_out, result, result2, result_lens])
    assure_folder(os.path.dirname(result))

    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_bedprocess -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_bedprocess --verbose'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False)


@pytest.mark.order(3)
def test_smk_do_uniq_counts(inst_folder):
    smk_cmd = 'do_uniq_counts'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result_outs = 'results/BED/FalseD25630/*lenuniqcount.csv'
    result2 = 'results/BED/FalseD25630/uniqcounts.tsv'
    expected2 = f'{inst_folder}test/expected/uniqcounts.tsv'
    #result3 = 'results/BED/FalseD25630/TEST_sample3_lenuniqcount.csv'
    #expected3 = f'{inst_folder}test/expected/TEST_sample3_lenuniqcount.csv'

    cleanup_files([cmd_out, result, result_outs, result2])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)
    #assert filecmp.cmp(result3, expected3, shallow=False), msg_diff_content(result3, expected3, cmd_out)


@pytest.mark.order(4)
def test_smk_do_global_length(inst_folder):
    cmd_out = 'results/test_smk_do_global_length.log'
    result = 'results/test_smk_do_global_length_status.log'
    expected = f'{inst_folder}test/expected/test_smk_do_global_length_status.log'
    result2 = 'results/feature/FalseD25630/global_length.tsv'
    expected2 = f'{inst_folder}test/expected/global_length.tsv'
    cleanup_files([cmd_out, result, result2])
    assure_folder(os.path.dirname(result))

    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_global_length -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_global_length --verbose'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False)


@pytest.mark.order(5)
def test_smk_do_bedprocess_szsel(inst_folder):
    szsel_suff = 'szsel_90_150/'
    cmd_out = 'results/test_smk_do_bedprocess_szsel.log'
    result = 'results/test_smk_do_bedprocess_status_szsel.log'
    expected = f'{inst_folder}test/expected/{szsel_suff}test_smk_do_bedprocess_status_szsel.log'
    result2 = f'results/BED/FalseD25630/{szsel_suff}mergeddf.csv'
    expected2 = f'{inst_folder}test/expected/{szsel_suff}mergeddf.csv'
    cleanup_files([cmd_out, result, f'results/BED/FalseD25630/{szsel_suff}'])
    assure_folder(os.path.dirname(result))

    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline_szsel.yaml -j 2 do_bedprocess -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline_szsel.yaml -j 2 do_bedprocess --verbose > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)


@pytest.mark.order(6)
def test_smk_do_mixed_model(inst_folder):
    # for Mac you need to add Rscript path to your PATH like that:
    # export PATH="/Library/Frameworks/R.framework/Resources:$PATH"
    # \.[\d-e]*\b  or better \.\d*(e-\d*)?\b
    # sed -i.bak 's/\.[\d-e]*\b/\./g' mixmod_results.tsv.txt
    smk_cmd = 'do_mixed_model'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    result2 = 'results/feature/FalseD25630/mixed_model/TEST_sample1_mixmod.pdf'
    expected2_min_size = 11000
    result3 = 'results/feature/FalseD25630/mixed_model/TEST_sample3_mixmod.pdf'
    expected3_min_size = 11000
    result4 = 'results/feature/FalseD25630/mixed_model/mixmod_results.tsv'
    expected4 = f'{inst_folder}test/expected/mixmod_results.tsv'
    cleanup_files([cmd_out, 'results/feature/FalseD25630/mixed_model'])  # result, result2, result3, result4])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    expected = f'{inst_folder}test/expected/test_smk_do_mixed_model_status.log'
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    result3_size = os.path.getsize(result3)
    print(f'result2_size={result2_size}, result3_size={result3_size}')
    sub_in_file(result4, r'\.(\d{5})\d*(e-\d*)?\b', r'.\1\2', '.bak')  # cut precision
    assert filecmp.cmp(result4, expected4, shallow=False), msg_diff_content(result4, expected4, cmd_out)
    assert result2_size >= expected2_min_size
    assert result3_size >= expected3_min_size
    # assert result4_size >= expected4_min_size


@pytest.mark.order(17)
def test_smk_do_hgmm(inst_folder):
    smk_rule = 'do_hgmm'
    cmd_out = f'results/test_smk_{smk_rule}.log'
    result = f'results/test_smk_{smk_rule}_status.log'
    #result2 = 'results/feature/FalseD25630/mixed_model/TEST_sample1_mixmod.pdf'
    #expected2_min_size = 11000
    #result3 = 'results/feature/FalseD25630/mixed_model/TEST_sample3_mixmod.pdf'
    #expected3_min_size = 11000
    result4 = 'results/feature/FalseD25630/mixed_model/hgmm_results.tsv'
    expected4 = f'{inst_folder}test/expected/hgmm_results.tsv'
    cleanup_files([cmd_out, result, result4])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_rule} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    expected = f'{inst_folder}test/expected/test_smk_{smk_rule}_status.log'
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_rule} > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result4, expected4, shallow=False), msg_diff_content(result4, expected4, cmd_out)


# @pytest.mark.usefixtures('test_smk_do_bedprocess')
@pytest.mark.order(7)
def test_smk_do_wps(inst_folder):
    cmd_out = 'results/test_smk_do_wps.log'
    result = 'results/test_smk_do_wps_status.log'
    result2 = f'results/feature/FalseD25630/wps/genome_wide/TEST_sample1/genome_wide{WPS_OUT_SUFFIX}_chr20.tsv.gz'
    result3 = f'results/feature/FalseD25630/wps/genome_wide/TEST_sample3/genome_wide{WPS_OUT_SUFFIX}_chr20.tsv.gz'
    expected2_size = 10424600
    expected3_size = 10436000
    cleanup_files([cmd_out, result, result2, result3,
                   'results/BAM/FalseD25630/TEST_sample1.sortByCoord.bam',
                   'results/BAM/FalseD25630/TEST_sample1.sortByCoord.bai',
                   'results/BAM/FalseD25630/TEST_sample3.sortByCoord.bam',
                   'results/BAM/FalseD25630/TEST_sample3.sortByCoord.bai',
                   ])
    assure_folder(os.path.dirname(result))
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_wps -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    expected = f'{inst_folder}test/expected/test_smk_do_wps_status.log'
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_wps --verbose > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    result3_size = os.path.getsize(result3)
    print(f'result2_size={result2_size}, result3_size={result3_size}')
    assert result2_size >= expected2_size, f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size}'
    assert result3_size >= expected3_size, f'Unexpected file size {result3_size} of file {result3}! Expected >= {expected3_size}'


@pytest.mark.order(8)
def test_smk_tmad_do_cal_blacklist(inst_folder):
    cmd_out = 'results/test_smk_tMAD_do_cal_blacklist.log'
    result = 'results/test_smk_tMAD_do_cal_blacklist_status.log'
    expected = f'{inst_folder}test/expected/test_smk_tMAD_do_cal_blacklist_status.log'
    result2 = 'results/BED/FalseD25630/control_blacklist_1000Kbp_norm_TRUE.RDS'
    expected2_size = 120
    cleanup_files([cmd_out, result, result2])
    assure_folder(os.path.dirname(result))
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_cal_blacklist -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_cal_blacklist --verbose > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    print(f'result2_size={result2_size}')
    assert result2_size >= expected2_size, f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size}'


@pytest.mark.order(9)
def test_smk_tmad_do_cal_refsample(inst_folder):
    cmd_out = 'results/test_smk_tMAD_do_cal_RefSample.log'
    result = 'results/test_smk_tMAD_do_cal_RefSample_status.log'
    expected = f'{inst_folder}test/expected/test_smk_tMAD_do_cal_RefSample_status.log'
    result2 = 'results/BED/FalseD25630/Ref_TEST_sample1_normTRUE_1000kbp_hg38_count.csv'
    expected2_size = 71100
    cleanup_files([cmd_out, result, result2])
    assure_folder(os.path.dirname(result))
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_cal_RefSample -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_cal_RefSample --verbose > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    print(f'result2_size={result2_size}')
    assert result2_size >= expected2_size, f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size}'


@pytest.mark.order(10)
def test_smk_do_cal_t_MAD_forall(inst_folder):
    cmd_out = 'results/test_smk_do_cal_t_MAD_forall.log'
    result = 'results/test_smk_do_cal_t_MAD_forall_status.log'
    expected = f'{inst_folder}test/expected/test_smk_do_cal_t_MAD_forall_status.log'
    result2 = 'results/BED/FalseD25630/tMAD/tMAD_results.tsv'
    expected2 = f'{inst_folder}test/expected/tMAD_results.tsv'
    # result3 = 'results/BED/FalseD25630/tMAD/tMAD_bs1000_control_Ref_TEST_sample1.tsv'
    # expected3 = f'{inst_folder}test/expected/tMAD_bs1000_control_Ref_TEST_sample1.tsv'
    cleanup_files([cmd_out, result, result2, 'results/BED/FalseD25630/tMAD/'])
    assure_folder(os.path.dirname(result))
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_cal_t_MAD_forall -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)

    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 do_cal_t_MAD_forall --verbose > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)
    # assert filecmp.cmp(result3, expected3, shallow=False), msg_diff_content(result3, expected3, cmd_out)


@pytest.mark.order(11)
def test_smk_do_circle_finder(inst_folder):
    smk_cmd = 'do_circle_finder'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result2 = 'results/BAM/microDNA/TEST_sample1/TEST_sample1.concordant_freqGr3.txt'
    expected2 = f'{inst_folder}test/expected/microDNA/TEST_sample1.concordant_freqGr3.txt'
    result3 = 'results/BAM/microDNA/TEST_sample1/TEST_sample1.microDNA-JT.txt'
    result4 = 'results/BAM/microDNA/TEST_sample3/TEST_sample3.disc.txt'
    expected4 = f'{inst_folder}test/expected/microDNA/TEST_sample3.disc.txt'
    result5 = 'results/BAM/microDNA/TEST_sample3/TEST_sample3.microDNA-JT.txt'
    cleanup_files([cmd_out, 'results/BAM/microDNA'])  # result, result2, result3, result4])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)
    assert os.path.exists(result3)
    assert filecmp.cmp(result4, expected4, shallow=False), msg_diff_content(result4, expected4, cmd_out)
    assert os.path.exists(result5)


# @pytest.mark.skip(reason="test data has not enough reads to produce output results") # go as further as possible
@pytest.mark.order(12)
def test_smk_do_ichorCNA(inst_folder):
    smk_rule = 'do_ichorCNA'
    cmd_out = f'results/test_smk_{smk_rule}.log'
    cmd_out_np = f'results/test_smk_{smk_rule}_np.log'
    result = f'results/test_smk_{smk_rule}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_rule}_status.log'
    # result2 = 'results/feature/FalseD25630/ichorCNA/TEST_sample1/TEST_sample1.cna.seg'
    # expected2_size = 300
    result2 = 'results/feature/FalseD25630/ichorCNA/TEST_sample1.wig'
    expected2 = f'{inst_folder}test/expected/ichorCNA/TEST_sample1.wig'
    # result3 = 'results/feature/FalseD25630/ichorCNA/TEST_sample1/TEST_sample3.cna.seg'
    # expected3_size = 300
    result3 = 'results/feature/FalseD25630/ichorCNA/TEST_sample3.wig'
    expected3 = f'{inst_folder}test/expected/ichorCNA/TEST_sample3.wig'
    result4 = 'results/feature/FalseD25630/ichorCNA/pon_wigs.txt'
    expected4 = f'{inst_folder}test/expected/ichorCNA/pon_wigs.txt'
    cleanup_files([cmd_out, result, 'results/feature/FalseD25630/ichorCNA',
                   'results/BAM/FalseD25630/*.sortByCoord.bam.bai'])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_rule} -np > {cmd_out_np}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out_np, result)
    # assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_rule} > {cmd_out}'
    run_cmd(smk_cmd)
    # assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)
    # assert filecmp.cmp(result3, expected3, shallow=False), msg_diff_content(result3, expected3, cmd_out)
    # assert filecmp.cmp(result4, expected4, shallow=False), msg_diff_content(result4, expected4, cmd_out)
    # due to not enough reads in all chromosomes to produce output results we take an finished esample
    #shutil.copytree(f'{inst_folder}test/expected/ichorCNA/TEST_sample1', 'results/feature/FalseD25630/ichorCNA/TEST_sample1')
    #shutil.copytree(f'{inst_folder}test/expected/ichorCNA/TEST_sample3', 'results/feature/FalseD25630/ichorCNA/TEST_sample3')
    #subprocess.run('touch results/feature/FalseD25630/ichorCNA/TEST_sample1/*', shell=True, capture_output=True, text=True)
    #subprocess.run('touch results/feature/FalseD25630/ichorCNA/TEST_sample3/*', shell=True, capture_output=True, text=True)
    # Path('results/feature/FalseD25630/ichorCNA/TEST_sample1/*').touch()
    # Path('results/feature/FalseD25630/ichorCNA/TEST_sample3/*').touch()


@pytest.mark.order(13)
def test_smk_do_ndetective_logo(inst_folder):
    smk_cmd = 'do_ndetective_logos'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result2 = 'results/BED/FalseD25630/logos/logo_TEST_sample1_all.pdf'
    expected2_size = 27000
    result3 = 'results/BED/FalseD25630/logos/logo_TEST_sample3_all.pdf'
    expected3_size = 27000
    result4 = 'results/BED/FalseD25630/logos/logo_TEST_sample1_point.csv.gz'
    expected4_size = 12300
    cleanup_files([cmd_out, 'results/BED/FalseD25630/logos'])  # result, result2)
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    assert result2_size >= expected2_size, f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size}'
    result3_size = os.path.getsize(result3)
    assert result3_size >= expected3_size, f'Unexpected file size {result3_size} of file {result3}! Expected >= {expected3_size}'
    result4_size = os.path.getsize(result4)
    assert result4_size >= expected4_size, f'Unexpected file size {result4_size} of file {result3}! Expected >= {expected4_size}'


@pytest.mark.order(14)
def test_smk_do_fft_bins(inst_folder):
    smk_cmd = 'do_fft_bins'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result2 = f'results/feature/FalseD25630/fft/TEST_sample1/TEST_sample1{WPS_OUT_SUFFIX}_fft_bins.tsv.gz'
    result2s = 'results/feature/FalseD25630/fft/TEST_sample1/*_fft_bins.ts*'
    expected2_size = 1500
    result3 = f'results/feature/FalseD25630/fft/TEST_sample3/TEST_sample3{WPS_OUT_SUFFIX}_fft_bins.tsv.gz'
    result3s = 'results/feature/FalseD25630/fft/TEST_sample3/*_fft_bins.ts*'
    expected3_size = 1200
    cleanup_files([cmd_out, result, result2s, result3s])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    result3_size = os.path.getsize(result3)
    assert (result2_size >= expected2_size) and (result3_size >= expected3_size), f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size} or \n' \
                                                                                  f'Unexpected file size {result3_size} of file {result3}! Expected >= {expected3_size}'


@pytest.mark.order(15)
def test_smk_do_fft_anno(inst_folder):
    smk_cmd = 'do_fft_anno'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result2 = f'results/feature/FalseD25630/fft/TEST_sample1/TEST_sample1_genes_body{WPS_OUT_SUFFIX}_fft.tsv.gz'
    result2s = 'results/feature/FalseD25630/fft/TEST_sample1/*_fft.ts*'
    expected2_size = 11000
    result3 = f'results/feature/FalseD25630/fft/TEST_sample3/TEST_sample3_genes_body{WPS_OUT_SUFFIX}_fft.tsv.gz'
    result3s = 'results/feature/FalseD25630/fft/TEST_sample3/*_fft.ts*'
    expected3_size = 22400
    cleanup_files([cmd_out, result, result2s, result3s])  # , result4, result5])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    result3_size = os.path.getsize(result3)
    assert (result2_size >= expected2_size) and (result3_size >= expected3_size), \
        f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size} or \n' \
        f'Unexpected file size {result3_size} of file {result3}! Expected >= {expected3_size}'


@pytest.mark.order(16)
def test_smk_do_fft_anno_plots(inst_folder):
    smk_cmd = 'do_fft_anno_plots'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result2 = f'results/feature/FalseD25630/fft/TEST_sample1/TEST_sample1_genes_body{WPS_OUT_SUFFIX}_fft_plots.pdf'
    expected2_size = 3000
    result3 = f'results/feature/FalseD25630/fft/TEST_sample3/TEST_sample3_genes_body{WPS_OUT_SUFFIX}_fft_plots.pdf'
    expected3_size = 3000
    cleanup_files([cmd_out, result, result2, result3])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    result2_size = os.path.getsize(result2)
    result3_size = os.path.getsize(result3)
    assert (result2_size >= expected2_size) and (result3_size >= expected3_size), \
        f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size} or \n' \
        f'Unexpected file size {result3_size} of file {result3}! Expected >= {expected3_size}'


@pytest.mark.order(18)
def test_smk_do_merge_results(inst_folder):
    smk_cmd = 'do_merge_results'
    cmd_out = f'results/test_smk_{smk_cmd}.log'
    result = f'results/test_smk_{smk_cmd}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_cmd}_status.log'
    result2 = 'results/feature/FalseD25630/all_results.tsv'
    expected2 = f'{inst_folder}test/expected/all_results.tsv'
    # expected2_size = 2700
    cleanup_files([cmd_out, result, result2])  # result, result2, result3, result4])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_cmd} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_cmd} > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)
    # result2_size = os.path.getsize(result2)
    # assert result2_size >= expected2_size, f'Unexpected file size {result2_size} of file {result2}! Expected >= {expected2_size}'


@pytest.mark.order(19)
def test_smk_do_ichorCNA_results(inst_folder):
    smk_rule = 'do_ichorCNA_results'
    cmd_out = f'results/test_smk_{smk_rule}.log'
    result = f'results/test_smk_{smk_rule}_status.log'
    expected = f'{inst_folder}test/expected/test_smk_{smk_rule}_status.log'
    result2 = 'results/feature/FalseD25630/ichorCNA/ichorCNA_results.tsv'
    expected2 = f'{inst_folder}test/expected/ichorCNA/ichorCNA_results.tsv'
    cleanup_files([cmd_out, result, result2])
    smk_cmd_np = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 {smk_rule} -np > {cmd_out}'
    run_cmd(smk_cmd_np)
    get_smk_output_status(cmd_out, result)
    assert filecmp.cmp(result, expected, shallow=False), msg_diff_content(result, expected, cmd_out)
    smk_cmd = f'snakemake -s {inst_folder}Snakefile --configfile {inst_folder}test/test_cfDNA_pipeline.yaml -j 2 --verbose {smk_rule} > {cmd_out}'
    run_cmd(smk_cmd)
    assert filecmp.cmp(result2, expected2, shallow=False), msg_diff_content(result2, expected2, cmd_out)


@pytest.mark.order(20)
def test_smk_do_hgmm_bins(inst_folder):
    smk_rule = 'do_hgmm_bins'
    results = ['results/feature/FalseD25630/mixed_model/hgmm_bins_TEST_sample1.tsv',
              'results/feature/FalseD25630/mixed_model/hgmm_bins_TEST_sample3.tsv']
    expecteds = None
    expected_sizes = [500, 300]
    run_test_smk(inst_folder, smk_rule, results, expecteds, expected_sizes)


@pytest.mark.order(21)
def test_smk_do_wps_to_bigWig(inst_folder):
    smk_rule = 'do_wps_to_bigWig'
    results = [f'results/feature/FalseD25630/wps/genome_wide/TEST_sample1_wps{WPS_OUT_SUFFIX}{"_norm" if WPS_BIGWIG_NORMALIZED else ""}.bw',
               f'results/feature/FalseD25630/wps/genome_wide/TEST_sample3_wps{WPS_OUT_SUFFIX}{"_norm" if WPS_BIGWIG_NORMALIZED else ""}.bw']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[10436000, 10487000])


@pytest.mark.order(22)
def test_smk_do_wps_to_bigWig_bins(inst_folder):
    smk_rule = 'do_wps_to_bigWig_bins'
    results = [f'results/feature/FalseD25630/wps/genome_wide/TEST_sample1_wps{WPS_OUT_SUFFIX}_bins1000_nmean.bw',
               f'results/feature/FalseD25630/wps/genome_wide/TEST_sample3_wps{WPS_OUT_SUFFIX}_bins1000_nmean.bw']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[60850, 63000])


@pytest.mark.order(23)
def test_smk_do_fft_to_bigWig(inst_folder):
    smk_rule = 'do_fft_to_bigWig'
    results = [f'results/feature/FalseD25630/fft/TEST_sample1{WPS_OUT_SUFFIX}_fft.bw',
               f'results/feature/FalseD25630/fft/TEST_sample3{WPS_OUT_SUFFIX}_fft.bw']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[19000, 19000])


# @pytest.mark.skip(reason="deprecated command: use configure.py")
def test_smk_do_fft_bins_py(inst_folder):
    smk_rule = 'do_fft_bins_py'
    results = [f'results/feature/FalseD25630/fft/TEST_sample1/TEST_sample1{WPS_OUT_SUFFIX}_fft_bins{FFT_BINSIZE}.tsv.gz',
               f'results/feature/FalseD25630/fft/TEST_sample3/TEST_sample3{WPS_OUT_SUFFIX}_fft_bins{FFT_BINSIZE}.tsv.gz']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[15000, 12000])


@pytest.mark.skip(reason="deprecated command: with this test data the computed CG.freq.txt is empty")
def test_smk_do_bam_split_GCBias(inst_folder):
    smk_rule = 'do_bam'
    cleanup_files(['bam_split_GCBiasCorrect_merge.log', 'results/BAM/TrueD25630'])
    results = [f'results/BAM/TrueD25630/TEST_sample1.sortByCoord.bam',
               f'results/BAM/TrueD25630/TEST_sample3.sortByCoord.bam']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[220000, 220000],
                 config_file="test/test_cfDNA_pipeline_GCbias.yaml")

@pytest.mark.order(24)
def test_smk_do_contig_summary(inst_folder):
    # prerequisites: normally done before
    # run_test_smk(inst_folder, 'do_bam', [], check_status_log=False)
    # run_test_smk(inst_folder, 'do_bai', [], check_status_log=False)
    # do
    smk_rule = 'do_contig_summary'
    results = [f'results/feature/FalseD25630/contig_summary.tsv']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[100], check_status_log=False)

@pytest.mark.order(25)
def test_smk_do_length_binning(inst_folder):
    # prerequisites: normally done before
    # run_test_smk(inst_folder, 'do_bam', [], check_status_log=False)
    # run_test_smk(inst_folder, 'do_bai', [], check_status_log=False)
    # do
    # if needed:
    # unlock(inst_folder)
    smk_rule = 'do_length_binning'
    results = [f'results/feature/FalseD25630/length/TEST_sample1_binned.tsv',
               f'results/feature/FalseD25630/length/TEST_sample3_binned.tsv']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[100, 100], check_status_log=False)

@pytest.mark.order(26)
def test_smk_do_length_hist(inst_folder):
    # prerequisites: normally done before
    # run_test_smk(inst_folder, 'do_bam', [], check_status_log=False)
    # run_test_smk(inst_folder, 'do_bai', [], check_status_log=False)
    # do
    smk_rule = 'do_length_hist'
    results = [f'results/feature/FalseD25630/length/TEST_sample1_hist.csv',
               f'results/feature/FalseD25630/length/TEST_sample3_hist.csv']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[1000, 800], check_status_log=False)

@pytest.mark.order(28)
def test_smk_do_counts_sl_ratio(inst_folder):
    # prerequisites: normally done before
    # run_test_smk(inst_folder, 'do_length_binning', [], check_status_log=False)
    # do
    # if needed:
    # unlock(inst_folder)
    smk_rule = 'do_counts_sl_ratio'
    results = [f'results/feature/FalseD25630/length/TEST_sample1_counts.tsv',
               f'results/feature/FalseD25630/length/TEST_sample3_counts.tsv',
               f'results/feature/FalseD25630/length/TEST_sample1_SLratio.tsv',
               f'results/feature/FalseD25630/length/TEST_sample3_SLratio.tsv']
    run_test_smk(inst_folder, smk_rule, results, expected_sizes=[15800, 15800, 16500, 16500], check_status_log=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Automated tests")
    parser.add_argument("--inst_folder", default="./",
                        help="Installation folder. Default: './' , used for testing from sources, "
                             "'/' would be for testing from a singularity image")

    args = parser.parse_args()
    # for singularity image tests
    inst_folder = args.inst_folder if args.inst_folder.endswith("/") else args.inst_folder + "/"

    test_configure(inst_folder)
    test_configure2(inst_folder)
    test_main(inst_folder)
    test_wps(inst_folder)
    test_size_selection_py(inst_folder)
    test_size_selection_sh(inst_folder)
    test_smk_do_preprocess(inst_folder)
    test_smk_do_bedprocess(inst_folder)
    test_smk_do_bedprocess_szsel(inst_folder)
    test_smk_do_mixed_model(inst_folder)
    test_smk_do_wps(inst_folder)  # takes too much time: > 7 min
    test_smk_tmad_do_cal_blacklist(inst_folder)
    test_smk_tmad_do_cal_refsample(inst_folder)
    test_smk_do_cal_t_MAD_forall(inst_folder)
    test_smk_do_global_length(inst_folder)
