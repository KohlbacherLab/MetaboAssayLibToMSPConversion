"""
author: Oliver Alka
date: 05.07.2019

Convert the metabolite assay library from OpenMS::AssayGeneratorMetabo [tsv, pqp, traML] to an msp file [msp].

Description:
1) rename headers
2) remove decoys (if available)
3) drop unused columns
4) set product charge to 1
5) use same adducts for product ions as for the precursor
6) recalculate rt in minutes
8) export
"""

# packages
import os
import tempfile
import click
import pandas as pd
from pyopenms import *

class entry:
    def __init__(self,
                 name = "",
                 precursor_mz = 0.0,
                 precursor_type = "", # Adduct
                 instrumenttype = "",
                 instrument = "",
                 authors = "",
                 license = "",
                 smiles = "",
                 inchi = "",
                 inchikey = "",
                 collisionenergy = "",
                 formula = "",
                 retentiontime = 0.0,
                 ionmode = "",
                 massbankaccession = "",
                 links = "",
                 comment = "",
                 numberoffragments = 0,
                 fragment_mz = [],
                 fragment_int = []):
        self.name = name
        self.precursor_mz = precursor_mz
        self.precursor_type = precursor_type 
        self.instrumenttype = instrumenttype
        self.instrument = instrument
        self.authors = authors
        self.license = license
        self.smiles = smiles
        self.inchi = inchi
        self.inchikey = inchikey
        self.collisionenergy = collisionenergy
        self.formula = formula
        self.retentiontime = retentiontime
        self.ionmode = ionmode
        self.massbankaccession = massbankaccession
        self.links = links
        self.comment = comment
        self.numberoffragments = numberoffragments
        self.fragment_mz = fragment_mz
        self.fragment_int = fragment_int

# function to reformat adduct string M+H+ (OpenMS) -> [M+H] (Skyline)
def reformatAdduct(adduct):
    charge = adduct[-1]
    adduct = adduct[:-1]
    adduct = '[' + adduct + ']' + charge
    return adduct

def fillTmpTSVWithValidTargetedExp(openmslib, tmpfile):
    # check file extension, validate, convert to tsv (if necessary)
    filename, extension = os.path.splitext(openmslib)
    targeted_exp = TargetedExperiment()
    
    if extension == '.pqp':
        TransitionPQPFile().convertPQPToTargetedExperiment(openmslib.encode(), targeted_exp, False)
    elif extension == ".traML" or extension == ".TraML" or extension == ".traml":
        TraMLFile().load(openmslib.encode(), targeted_exp)
    else:
        filetype = FileTypes().nameToType('TSV')
        TransitionTSVFile().convertTSVToTargetedExperiment(openmslib.encode(), filetype, targeted_exp)

    # check validity of OpenMS::TargetedExperiment
    TransitionTSVFile().validateTargetedExperiment(targeted_exp)
    TransitionTSVFile().convertTargetedExperimentToTSV(tmpfile.name.encode(), targeted_exp)
    
    return tmpfile

@click.command()
@click.option('--openmslib', '-in', multiple = False, type = click.Path(), help = 'Input assay library from OpenMS::AssayGeneratorMetbo (.tsv,.traML,.pqp)')
@click.option('--msp', '-out', multiple = False, type = click.Path(), help = 'Output spectral libray (.msp)')
@click.option('--removedecoys/--no-removedecoys', default=True, help = 'Removes the decoys from assay library (default: True)')

def main(openmslib, msp, removedecoys):

    tmpfile = tempfile.NamedTemporaryFile(suffix='.tsv')
    fillTmpTSVWithValidTargetedExp(openmslib, tmpfile)

    # read input 
    library = pd.read_csv(tmpfile.name, sep='\t')

    # remove decoys, since there are problems with the import into skyline! 
    if removedecoys:
    	indexDecoys = library[ library['Decoy'] == 1 ].index
    	# delete these row indexes from dataFrame
    	library.drop(indexDecoys, inplace=True)

    # drop unused columns (from AssayGeneratorMetabo output))
    library = library.drop([ 
        'PeptideSequence',
        'ModifiedPeptideSequence',
        'PeptideGroupLabel',
        'ProteinId',
        'UniprotId',
        'GeneName',
        'FragmentType',
        'FragmentSeriesNumber',
        'PrecursorIonMobility',
        'TransitionId',
        'DetectingTransition',
        'IdentifyingTransition',
        'QuantifyingTransition',
        'Decoy',
        'Peptidoforms'
        ], axis=1)

    # since currently only charge one features are used in OpenMS
    # the column 'ProductCharge" is set to 1 
    library['ProductCharge'] = 1

    # the adducts are represented differently OpenMS: M+H+ , Skyline [M+H]
    # reformat the adduct 
    library['Adducts'] = library['Adducts'].apply(reformatAdduct)

    # copy the precursor adducts into ProductAdduct column for calculation in Skyline
    library['Adducts'] = library['Adducts']

    # recalulate the rt in minutes
    library['NormalizedRetentionTime'] = library['NormalizedRetentionTime']/60

    # Go over the file and write msp
    # Construct MSP entries
    msp_entries = []
    msp_entry = entry()
    tgid = None
    for index, row in library.iterrows():
        current_tgid = row['TransitionGroupId']
        last_row = len(library.index) - 1
        if tgid == None or tgid == current_tgid and index != last_row :
            # add relevant entries
            msp_entry.name = row["CompoundName"]
            msp_entry.precursor_mz = row['PrecursorMz']
            msp_entry.precursor_type = row['Adducts']
            msp_entry.formula = row['SumFormula']
            msp_entry.retentiontime = row['NormalizedRetentionTime']
            # infer ion mode based on adduct charge
            adduct = row['Adducts']
            if adduct[-1] == "+":
                msp_entry.ionmode = "Positive"
            else:
                msp_entry.ionmode = "Negative"
            # add fragment mz / int
            msp_entry.fragment_mz.append(row['ProductMz'])
            msp_entry.fragment_int.append(row['LibraryIntensity'])
            # add fragment peak annotation as comment
            msp_entry.comment += str(row['Annotation'] + " ")
        elif index == last_row: #last row 
            msp_entry.fragment_mz.append(row['ProductMz'])
            msp_entry.fragment_int.append(row['LibraryIntensity'])
            msp_entry.comment += str(row['Annotation'] + " ")
            msp_entry.numberoffragments = len(msp_entry.fragment_mz)
            msp_entries.append(msp_entry)
            msp_entry = entry()
        else:
            current_tgid = row['TransitionGroupId']
            # fill number-of-fragments
            msp_entry.numberoffragments = len(msp_entry.fragment_mz)
            # save entry 
            msp_entries.append(msp_entry)
            # begin new fr fragment_int 
            msp_entry = entry()
            msp_entry.fragment_mz = [row['ProductMz']]
            msp_entry.fragment_int = [row['LibraryIntensity']]
            msp_entry.comment += str(row['Annotation'] + " ")
        tgid = current_tgid

    # Write msp file based on datastructure / msp_entry! 
    if os.path.exists(msp):
        print("File already exists, new spectral library entries will be added to the end./n")
    else:
        print("New file will be created: " + msp)

    with open(msp, 'a') as file:
        for element in msp_entries:
            file.write("NAME: " + element.name + "\n")
            file.write("PRECURSORMZ: " + str(element.precursor_mz) + "\n")
            file.write("PRECURSORTYPE: " + element.precursor_type + "\n")
            file.write("INSTRUMENTTYPE: " + element.instrumenttype + "\n")
            file.write("INSTRUMENT: " + element.instrument + "\n")
            file.write("Authors: " + element.authors + "\n")
            file.write("License: " + element.license + "\n")
            file.write("SMILES: " + element.smiles + "\n")
            file.write("INCHI: " + element.inchi + "\n")
            file.write("INCHIKEY: " + element.inchikey + "\n")
            file.write("COLLISIONENERGY: " + element.collisionenergy + "\n")
            file.write("FORMULA: " + element.formula + "\n")
            file.write("RETENTIONTIME: " + str(element.retentiontime) + "\n")
            file.write("IONMODE: " + str(element.ionmode) + "\n")
            file.write("MASSBANKACCESSION: " + element.massbankaccession + "\n")
            file.write("Links: " + element.links + "\n")
            file.write("Comment: " + element.comment + "\n")
            file.write("Num Peaks: " + str(element.numberoffragments) + "\n")
            for i in range(0,len(element.fragment_mz)):
                file.write(str(element.fragment_mz[i]) + " " + str(element.fragment_int[i]) + "\n")
            file.write("\n")

    # remove temporary file 
    tmpfile.close()

    print("Export successful")


if __name__ == "__main__":
    main()
