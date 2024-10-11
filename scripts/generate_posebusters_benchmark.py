import os
import gzip
import shutil
import string
from pathlib import Path

from Bio.PDB import PDBList, MMCIFParser, PDBIO, Select
from Bio.PDB import PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from rdkit import Chem

# Special thanks to ChatGPT for organizing the code!

# ======================= Configuration =======================

# Input file containing PDB codes and ligand codes
posebuster_file = Path('./posebuster_benchmark.txt')

# Directories
download_dir = Path('./structures/raw_from_pdb')           # Where PDB files are downloaded
cleaned_structures_dir = Path('./structures/all')          # Where cleaned PDB files are saved
monomer_structures_dir = Path('./structures/monomer')      # Where monomer PDB files are copied

protein_fasta_dir = Path('./sequences/all/no_ligand')      # Where protein-only FASTAs are saved
modified_fasta_dir = Path('./sequences/all/ligand')        # Where modified FASTAs with ligand SMILES are saved

protein_monomer_fasta_dir = Path('./sequences/monomer/no_ligand')  # Where monomer protein-only FASTAs are copied
modified_monomer_fasta_dir = Path('./sequences/monomer/ligand')    # Where monomer modified FASTAs are copied

master_fasta_path = Path('./sequences/monomer_combined/posebusters_monomers_combined.fasta')  # Master FASTA file

# Create directories if they don't exist
download_dir.mkdir(parents=True, exist_ok=True)
cleaned_structures_dir.mkdir(parents=True, exist_ok=True)
monomer_structures_dir.mkdir(parents=True, exist_ok=True)

protein_fasta_dir.mkdir(parents=True, exist_ok=True)
modified_fasta_dir.mkdir(parents=True, exist_ok=True)

protein_monomer_fasta_dir.mkdir(parents=True, exist_ok=True)
modified_monomer_fasta_dir.mkdir(parents=True, exist_ok=True)

master_fasta_path.parent.mkdir(parents=True, exist_ok=True)

# ======================= Functions =======================

class ProteinLigandSelect(Select):
    def __init__(self, ligand_code):
        self.ligand_code = ligand_code.upper()

    def accept_residue(self, residue):
        hetfield, resseq, icode = residue.id
        if hetfield == 'W':
            # Exclude water molecules
            return False
        elif hetfield.strip() == '' and residue.resname == 'HOH':
            # Exclude water molecules that might not be labeled with 'W'
            return False
        elif hetfield.strip() != '' and residue.resname != self.ligand_code:
            # Exclude heteroatoms that are not the specified ligand
            return False
        else:
            # Include protein residues and the specified ligand
            return True

def clean_structure(structure):
    for model in structure:
        for chain in model:
            for residue in chain:
                # Clean residue id
                hetfield, resseq, icode = residue.id
                if icode is None:
                    icode = ' '
                residue.id = (hetfield, resseq, icode)
                for atom in residue:
                    # Clean altloc
                    if not atom.altloc:
                        atom.altloc = ' '
    return structure

def shorten_chain_ids(structure):
    # Create a list of single-character chain IDs (uppercase, lowercase letters, digits)
    chain_id_chars = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    # Create a mapping from original chain IDs to single-character IDs
    chain_id_mapping = {}
    i = 0
    for model in structure:
        for chain in model:
            original_chain_id = chain.id
            if original_chain_id not in chain_id_mapping:
                if i < len(chain_id_chars):
                    chain_id_mapping[original_chain_id] = chain_id_chars[i]
                    i += 1
                else:
                    # If we run out of single-character IDs, raise an error
                    raise ValueError("Too many chains to assign unique single-character IDs.")
            # Update the chain ID
            chain.id = chain_id_mapping[original_chain_id]
    return chain_id_mapping

def read_pdb_ligand_list(posebuster_file):
    """
    Reads the PDB codes and ligand codes from the input file.
    """
    # Read input file
    with open(posebuster_file, 'r') as f:
        content = f.read()

    # Replace newlines with spaces and split on commas
    content = content.replace('\n', ' ')
    entries = content.split(',')

    pdb_ligand_list = []
    for entry in entries:
        entry = entry.strip()
        if entry:
            parts = entry.split()
            if len(parts) == 2:
                pdb_code, ligand_code = parts
                pdb_code = pdb_code.strip().lower()
                ligand_code = ligand_code.strip().upper()
                pdb_ligand_list.append((pdb_code, ligand_code))
            else:
                print(f"Invalid entry: {entry}")

    return pdb_ligand_list

def download_and_process_structures(pdb_ligand_list, download_dir, cleaned_structures_dir):
    """
    Downloads PDB structures and processes them to include only the protein and specified ligand.
    """
    # Create PDBList instance
    pdbl = PDBList()

    for pdb_code, ligand_code in pdb_ligand_list:
        print(f"Processing {pdb_code.upper()} with ligand {ligand_code}")
        # Download CIF file
        cif_file = pdbl.retrieve_pdb_file(
            pdb_code, file_format='mmCif', pdir=str(download_dir), overwrite=True
        )
        cif_file = Path(cif_file)
        # Decompress if needed
        if cif_file.suffix == '.gz':
            with gzip.open(cif_file, 'rb') as f_in:
                with open(cif_file.with_suffix(''), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            cif_file.unlink()
            cif_file = cif_file.with_suffix('')

        # Parse structure
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_code, str(cif_file))
        # Clean structure
        clean_structure(structure)
        # Shorten chain IDs
        chain_id_mapping = shorten_chain_ids(structure)
        print("Chain ID mapping:")
        for original_id, new_id in chain_id_mapping.items():
            print(f"  Original: {original_id} -> New: {new_id}")
        # Prepare to save
        io = PDBIO()
        io.set_structure(structure)
        # Create the select object
        select = ProteinLigandSelect(ligand_code)
        # Output file name
        output_filename = cleaned_structures_dir / f"{pdb_code}_{ligand_code}.pdb"
        # Save the file
        io.save(str(output_filename), select=select)
        print(f"Saved to {output_filename}")

def is_monomer(pdb_file_path):
    """
    Determines if the PDB file contains only one chain.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('structure', str(pdb_file_path))
        chains = set()
        for model in structure:
            for chain in model:
                chains.add(chain.id)
        if len(chains) == 1:
            return True
        else:
            return False
    except Exception as e:
        print(f"Error parsing {pdb_file_path}: {e}")
        return False

def copy_monomer_pdbs(input_dir, output_dir):
    """
    Copies PDB files with only one chain from input_dir to output_dir.
    """
    pdb_files = [f for f in input_dir.iterdir() if f.is_file() and f.suffix == '.pdb']
    print(f"Found {len(pdb_files)} PDB files in {input_dir}")

    monomer_count = 0
    for pdb_file in pdb_files:
        pdb_file_path = pdb_file
        if is_monomer(pdb_file_path):
            shutil.copy(pdb_file_path, output_dir / pdb_file.name)
            monomer_count += 1
            print(f"Copied monomer PDB: {pdb_file.name}")

    print(f"Total monomer PDB files copied: {monomer_count}")

def generate_fastas(input_dir, output_dir1, output_dir2):
    """
    Processes PDB files to generate protein FASTAs and modified FASTAs with ligand SMILES.
    """
    # Get list of PDB files
    pdb_files = [f for f in input_dir.iterdir() if f.is_file() and f.suffix == '.pdb']

    for pdb_file in pdb_files:
        pdb_path = pdb_file
        base_name = pdb_file.stem  # Remove '.pdb' extension
        if '_' not in base_name:
            print(f"Unexpected file name format: {pdb_file.name}")
            continue
        pdb_code, lig_code = base_name.rsplit('_', 1)

        # Extract protein sequences for each chain
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_code, str(pdb_path))
        ppb = PPBuilder()

        chain_records = []
        for model in structure:
            for chain in model:
                chain_id = chain.id
                # Build peptides for this chain
                peptides = ppb.build_peptides(chain)
                if peptides:
                    # Concatenate sequences from all peptides in this chain
                    chain_sequence = ''.join([str(pp.get_sequence()) for pp in peptides])
                    chain_record = SeqRecord(Seq(chain_sequence),
                                             id=f"{pdb_code}_{chain_id}",
                                             description=f"Chain {chain_id}")
                    chain_records.append(chain_record)
                    # For the modified FASTA
                    chain_record_mod = SeqRecord(Seq(chain_sequence),
                                                 id=f"protein|{pdb_code}_{chain_id}",
                                                 description=f"Chain {chain_id}")
                else:
                    print(f"No peptides found for chain {chain_id} in {pdb_code}")
                    continue

        if not chain_records:
            print(f"No protein chains found in {pdb_file.name}")
            continue

        # Save individual protein FASTA (all chains)
        output_path1 = output_dir1 / f"{pdb_code}.fasta"
        SeqIO.write(chain_records, str(output_path1), 'fasta')

        # Extract ligand atoms
        ligand_atoms = []
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    resname = line[17:20].strip()
                    if resname == lig_code:
                        ligand_atoms.append(line)
        if ligand_atoms:
            temp_ligand_pdb = Path(f"temp_{pdb_code}_{lig_code}.pdb")
            with open(temp_ligand_pdb, 'w') as f:
                f.write(''.join(ligand_atoms))

            # Generate SMILES using RDKit
            mol = Chem.MolFromPDBFile(str(temp_ligand_pdb), removeHs=False)
            if mol is None:
                print(f"Could not parse ligand from {pdb_file.name}")
                ligand_smiles = ''
            else:
                ligand_smiles = Chem.MolToSmiles(mol)

            # Remove temporary file
            temp_ligand_pdb.unlink()
        else:
            print(f"No ligand {lig_code} found in {pdb_file.name}")
            ligand_smiles = ''

        # Prepare modified records for the second output
        modified_records = []
        for chain_record in chain_records:
            chain_id = chain_record.id.split('_')[1]
            mod_record = SeqRecord(chain_record.seq,
                                   id=f"protein|{pdb_code}_{chain_id}",
                                   description=chain_record.description)
            modified_records.append(mod_record)

        # Create ligand record
        ligand_record = SeqRecord(Seq(ligand_smiles),
                                  id=f"ligand|{lig_code}",
                                  description='')

        # Write to the second output directory
        output_path2 = output_dir2 / f"{pdb_code}.fasta"
        with open(output_path2, 'w') as f:
            SeqIO.write(modified_records, f, 'fasta')
            SeqIO.write(ligand_record, f, 'fasta')

    print("Processing complete!")

def copy_monomer_fastas(protein_fasta_dir, protein_monomer_fasta_dir, modified_fasta_dir, modified_monomer_fasta_dir):
    """
    Copies monomer FASTA files to monomer directories.
    """
    # Process protein-only FASTAs
    for fasta_file in protein_fasta_dir.iterdir():
        if fasta_file.suffix in ['.fasta', '.fa']:
            fasta_path = fasta_file
            sequences = list(SeqIO.parse(str(fasta_path), 'fasta'))
            if len(sequences) == 1:
                # Copy to monomer directory
                shutil.copy(fasta_path, protein_monomer_fasta_dir / fasta_file.name)
                print(f"Copied monomer protein FASTA: {fasta_file.name}")

    # Process modified ligand FASTAs
    for fasta_file in modified_fasta_dir.iterdir():
        if fasta_file.suffix in ['.fasta', '.fa']:
            fasta_path = fasta_file
            sequences = list(SeqIO.parse(str(fasta_path), 'fasta'))
            # Count protein sequences (exclude ligand entries)
            protein_sequences = [seq for seq in sequences if seq.id.startswith('protein|')]
            if len(protein_sequences) == 1:
                # Copy to monomer directory
                shutil.copy(fasta_path, modified_monomer_fasta_dir / fasta_file.name)
                print(f"Copied monomer modified ligand FASTA: {fasta_file.name}")

    print("Monomer FASTA files have been copied to the respective directories.")

def combine_monomer_fastas(protein_monomer_fasta_dir, master_fasta_path):
    """
    Combines monomer protein FASTA files into a master FASTA file.
    """
    # List to hold all sequences
    all_sequences = []

    # Iterate over each FASTA file in the monomer directory
    for fasta_file in protein_monomer_fasta_dir.iterdir():
        if fasta_file.suffix in ['.fasta', '.fa']:
            fasta_path = fasta_file
            # Parse the FASTA file and append sequences to the list
            sequences = list(SeqIO.parse(str(fasta_path), 'fasta'))
            all_sequences.extend(sequences)

    # Write all sequences to the master FASTA file
    with open(master_fasta_path, 'w') as output_handle:
        SeqIO.write(all_sequences, output_handle, 'fasta')

    print(f"Master FASTA file created at: {master_fasta_path}")

# ======================= Main Execution =======================

def main():
    # Step 1: Read PDB codes and ligand codes
    pdb_ligand_list = read_pdb_ligand_list(posebuster_file)

    # Step 2: Download and process structures
    download_and_process_structures(pdb_ligand_list, download_dir, cleaned_structures_dir)

    # Step 3: Copy monomer PDB files
    copy_monomer_pdbs(cleaned_structures_dir, monomer_structures_dir)

    # Step 4: Generate FASTAs from cleaned structures
    generate_fastas(cleaned_structures_dir, protein_fasta_dir, modified_fasta_dir)

    # Step 5: Copy monomer FASTAs to monomer directories
    copy_monomer_fastas(protein_fasta_dir, protein_monomer_fasta_dir, modified_fasta_dir, modified_monomer_fasta_dir)

    # Step 6: Combine monomer protein FASTAs into a master FASTA file
    combine_monomer_fastas(protein_monomer_fasta_dir, master_fasta_path)

if __name__ == "__main__":
    main()
