import os
import csv
import glob
import pandas as pd
from . import get_reference_locations

def get_dataset_directory(dataset):
    """
    for a given dataset, get the directory where the Kallisto data is stored
    """
    kallisto_directory = get_reference_locations.get_kallisto_directory()
    return os.path.join(kallisto_directory, '%s_95' % dataset)

def get_fastq_locs():
    """
    get list of fastq files and their locations
    """
    
    direc = get_reference_locations.get_fastq_directory()
    fastqs = []
    
    for root, dirs, files in os.walk(direc):
        fastqs.append([os.path.join(root, file) for file in files if file.endswith('.gz')])
    
    fastqs = [file for files in fastqs for file in files]
    return fastqs

def get_fastq_cell(fastq_loc):
    """
    for a fastq file, get its corresponding cell
    """
    
    cell = fastq_loc.split('/')[-1].split('.fastq')[0]
    if cell.endswith('_1') or cell.endswith('_2'):
        return cell[:-2]
    return cell

def get_fastq_pairs(fastqs):
    """
    pair up fastq files
    """
    
    paired = set()
    fastqs = set(fastqs)
    
    while len(fastqs) > 0:
        item = fastqs.pop()
        if item.endswith('_1.fastq.gz'):
            item2 = item.replace('_1.fastq', '_2.fastq')
            fastqs.remove(item2)
            paired.add((item, item2))
        elif item.endswith('_2.fastq.gz'):
            item2 = item.replace('_2.fastq', '_1.fastq')
            fastqs.remove(item2)
            paired.add((item2, item))
        else:
            paired.add((item,))
    
    return list(paired)

def get_cell_locs():
    """
    for each cell, get the location of its fastq files
    """
    
    fastq_locs = get_fastq_locs()
    fastq_pairs = get_fastq_pairs(fastq_locs)
    cell_locs = {get_fastq_cell(fastq_pair[0]):fastq_pair for fastq_pair in fastq_pairs}
    
    return cell_locs

def get_target_cells(dataset, dataset_directory, skip_types=[]):
    """
    get list of cells for which kallisto needs to be run
    """
    
    # Step 1: Get all cells in a dataset
    fname = 'References/Cell_Labels_%s.tsv' % dataset
    df_cell = pd.read_csv(fname, sep='\t', header=0, index_col=0)

    # Step 1.5 If we want to skip some (len(skip_types) > 0), remove those
    if len(skip_types) > 0:
        df_cell = df_cell.loc[~df_cell.CellType.isin(skip_types),:]
    
    # Step 2: Remove cells that already have a Kallisto file
    already_generated = glob.glob(os.path.join(dataset_directory, '*'))
    already_generated = [cell.split('/')[-1] for cell in already_generated]
    df_cell = df_cell.loc[~df_cell.index.isin(already_generated)]
    
    # convert to a dictionary of targets
    target_cells = {cell:df_cell.loc[[cell], 'SRR'].tolist() for cell in set(df_cell.index)}
    
    return target_cells

def generate_kallisto_commands(dataset_directory, target_cells, umi=False):
    """
    generate commands to run kallisto for a list of cells
    """
    
    # initialize variables
    cell_locs = get_cell_locs()
    transcript = get_reference_locations.get_kallisto_reference_index()
    base_command_2 = 'kallisto quant -i %s -o %s -b 100 -t 8 %s'
    base_command_1 = 'kallisto quant -i %s -o %s --single -l 200 -s 20 -b 100 -t 8 %s'
    commands = []
    
    # generate command by cell
    for cell, SRRs in target_cells.items():
        # get target folder and fastq files
        cell_dir = '%s/%s' % (dataset_directory, cell)
        fastq_pairs = [cell_locs[SRR] for SRR in SRRs if SRR in cell_locs]
        if len(fastq_pairs) == 0:
            continue
        
        # asses if data is paired or not
        if len(fastq_pairs[0]) == 1:
            base_command = base_command_1
        else:
            base_command = base_command_2
            
        # split up fastq files into a connected string
        fastqs = [fastq for fastq_pair in fastq_pairs for fastq in fastq_pair]
        fastqs = ' '.join(fastqs)
        
        # generate command
        command = base_command % (transcript, cell_dir, fastqs)
        commands.append(command)
        
    return commands

def get_kallisto_commands(dataset, skip_types=[], umi=False):
    """
    generate commands to run kallisto for a dataset
    """
    
    # get list of target cells
    dataset_directory = get_dataset_directory(dataset)
    target_cells = get_target_cells(dataset, dataset_directory, skip_types=skip_types)
    
    # generate list of commands
    commands = generate_kallisto_commands(dataset_directory, target_cells, umi=umi)
    
    # save commands
    with open('Commands/kallisto_%s.sh' % dataset, 'w') as w:
        w.write('#!/bin/bash\n\n')
        w.write('\n'.join(commands))
    return

def get_kallisto_data(dataset, skip_types=[], umi=False):
    """
    For a given dataset, generate the kallisto data
    """
    
    # Step 1: Make sure the target directory exists
    dataset_directory = get_dataset_directory(dataset)
    if not os.path.isdir(dataset_directory):
        os.system('mkdir %s' % dataset_directory)
    
    # Step 2: Generate commands to file
    get_kallisto_commands(dataset, skip_types=skip_types, umi=umi)
    
    # Step 3: Execute bash file
    os.system('bash Commands/kallisto_%s.sh' % dataset)
    
    return
