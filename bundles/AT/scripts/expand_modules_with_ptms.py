#!/usr/bin/env python3
"""
expand_modules_with_ptms.py
===========================
Expand EpiProfile-PLANTS H2A/H2B/H1 MATLAB modules to include multiple
peptidoforms (PTM variants), following the same pattern as H3 modules.

Chemistry rules:
- N-terminal (position 0): always 'pr' (propionylation)
- Internal K (not C-terminal): 'acD3' if unmodified (deuterated propionylation),
  or the specific PTM (ac, me1, me2, me3) if modified
- C-terminal K or R: NO modification (trypsin cleavage site)
- N-terminal K (position 0): gets BOTH 'pr' AND 'acD3'/'ac' etc.
- mod_type uses 0-indexed PEPTIDE positions
- rt_ref = 0 for all entries

Author: Pelayo Gonzalez de Lena (@pelamovic)
Date: 2026-03-26
"""

import csv
import os
import re
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ── Paths ──────────────────────────────────────────────────────────────────
MANIFEST_PATH = Path(r"D:/Antigravity/repos/epiprofile_yuan_ecosistem/audit/plant_module_manifest.tsv")
GRAUBOVE_PATH = Path(r"D:/Antigravity/TESIS_PHD/docs_fuente/hPTMs_GrauBove_AT_PP_CR.tsv")
TIER3_DIR = Path(r"D:/Antigravity/revision_tesis/epiprofile_plants_expanded/bundles/AT/src/TIER3")


# ── Known PTM sites from Grau-Bove + literature ───────────────────────────
# Format: {(histone_class, protein_position): [list of PTM types]}
# PTM types: 'ac' (acetylation), 'me1', 'me2', 'me3' (methylation)
# We focus on K residues only (propionylation chemistry)

def parse_graubove_ptms(path: Path) -> Dict[Tuple[str, int], set]:
    """Parse Grau-Bove table to extract known PTMs per histone/position."""
    ptm_map = {}  # (histone_type, position) -> set of PTMs

    # Map from Grau-Bove modification names to EpiProfile mod codes
    mod_map = {
        'Ace': 'ac',
        'Me1': 'me1',
        'Me2': 'me2',
        'Me3': 'me3',
        # Cro (crotonylation) and Pho (phosphorylation) are NOT handled by
        # propionylation chemistry - skip them
    }

    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            species = row.get('Species', '').strip()
            if species != 'AT':
                continue

            histone_type = row.get('Histone_Type', '').strip()
            modification = row.get('Modification', '').strip()
            residue = row.get('Residue', '').strip()
            pos_ref = row.get('Pos_Reference', '').strip()

            # Only K modifications relevant for propionylation chemistry
            if residue != 'K':
                continue

            # Map modification
            if modification not in mod_map:
                continue
            mod_code = mod_map[modification]

            # Parse position - handle complex cases like '3-27-I', '10-11-12'
            # For ambiguous positions, we extract each individual position
            positions = set()
            for part in pos_ref.replace('-I', '').split('-'):
                part = part.strip()
                if part.isdigit():
                    positions.add(int(part))

            # Classify histone
            ht = histone_type.split()[0]  # 'H2A legacy' -> 'H2A'

            for pos in positions:
                key = (ht, pos)
                if key not in ptm_map:
                    ptm_map[key] = set()
                ptm_map[key].add(mod_code)

    return ptm_map


def parse_manifest(path: Path) -> List[dict]:
    """Parse the plant module manifest TSV."""
    modules = []
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            species = row.get('species', '').strip()
            if species != 'AT':
                continue

            histone_type = row.get('histone_type', '').strip()
            # Only process H2A, H2B, H1 (not H3, H4)
            if histone_type not in ('H2A', 'H2B', 'H1'):
                continue

            modules.append({
                'file_name': row['file_name'].strip(),
                'species': species,
                'histone_name': row.get('histone_name', '').strip(),
                'histone_type': histone_type,
                'pep_seq': row.get('pep_seq', '').strip(),
                'pep_len': int(row.get('pep_len', 0)),
                'pos_start': int(row.get('pos_start', 0)),
                'pos_end': int(row.get('pos_end', 0)),
                'n_mods': int(row.get('n_mods', 1)),
                'mod_shorts': row.get('mod_shorts', 'unmod').strip(),
                'charges': row.get('charges', '1 2').strip(),
                'variant_code': row.get('variant_code', '').strip(),
            })

    return modules


def find_modifiable_k_positions(pep_seq: str) -> List[Tuple[int, int]]:
    """
    Find K residues that can be modified (not C-terminal).

    Returns list of (peptide_0idx_position, is_n_terminal) tuples.
    C-terminal K and R are excluded.
    """
    modifiable = []
    last_idx = len(pep_seq) - 1

    for i, aa in enumerate(pep_seq):
        if aa == 'K' and i < last_idx:  # Not C-terminal
            modifiable.append(i)

    return modifiable


def map_peptide_to_protein_pos(pep_pos: int, pos_start: int) -> int:
    """Convert 0-indexed peptide position to protein position."""
    return pos_start + pep_pos


def determine_ptms_for_module(module: dict, graubove_ptms: dict) -> Dict[int, List[str]]:
    """
    Determine which PTMs each modifiable K can have.

    Returns: {peptide_0idx_pos: [list of PTM codes]}
    """
    pep_seq = module['pep_seq']
    pos_start = module['pos_start']
    histone_type = module['histone_type']

    # Map histone variants to Grau-Bove categories
    # H2A variants: canonical -> H2A, H2AX -> H2A, H2AZ -> H2AZ, H2AW -> H2A
    histone_name = module['histone_name']

    graubove_types = set()
    if 'H2AZ' in histone_name or 'H2A.Z' in histone_name:
        graubove_types = {'H2AZ', 'H2A'}
    elif 'H2AX' in histone_name or 'H2A.X' in histone_name:
        graubove_types = {'H2A'}
    elif 'H2AW' in histone_name or 'H2A.W' in histone_name:
        graubove_types = {'H2A'}
    elif 'H2A' in histone_name:
        graubove_types = {'H2A'}
    elif 'H2B' in histone_name:
        graubove_types = {'H2B'}
    elif 'H1' in histone_name:
        graubove_types = {'H1'}

    modifiable_ks = find_modifiable_k_positions(pep_seq)
    ptm_map = {}

    for pep_pos in modifiable_ks:
        protein_pos = map_peptide_to_protein_pos(pep_pos, pos_start)

        ptms = set()
        for ht in graubove_types:
            key = (ht, protein_pos)
            if key in graubove_ptms:
                ptms.update(graubove_ptms[key])

        # If no PTMs found in Grau-Bove, still offer acetylation as default
        # for K residues (very common in histones, conservative approach)
        # BUT only if the manifest already hints at it
        if not ptms:
            # Check manifest mod_shorts for this position
            manifest_mods = module['mod_shorts'].split('|')
            for mm in manifest_mods:
                if mm == 'unmod':
                    continue
                # Parse: K13ac, K5me2, etc.
                m = re.match(r'K(\d+)(ac|me1|me2|me3)', mm)
                if m and int(m.group(1)) == protein_pos:
                    ptms.add(m.group(2))

        if ptms:
            ptm_map[pep_pos] = sorted(ptms)

    return ptm_map


def generate_peptidoforms(pep_seq: str, ptm_map: Dict[int, List[str]],
                          pos_start: int) -> List[Tuple[str, str]]:
    """
    Generate all combinatorial peptidoforms.

    Returns list of (mod_short_name, mod_type_string) tuples.
    First entry is always 'unmod'.
    """
    modifiable_ks = sorted(ptm_map.keys())

    if not modifiable_ks:
        # No modifiable K with known PTMs - just return unmod
        mod_type = build_mod_type_string(pep_seq, {})
        return [('unmod', mod_type)]

    # Build all combinations: each K site can be unmodified or have any of its PTMs
    site_options = []
    for pep_pos in modifiable_ks:
        options = [None]  # None = unmodified (acD3)
        options.extend(ptm_map[pep_pos])
        site_options.append(options)

    peptidoforms = []

    for combo in product(*site_options):
        # combo is a tuple of (None or 'ac'/'me1'/'me2'/'me3') for each K site
        mod_dict = {}
        for pep_pos, ptm in zip(modifiable_ks, combo):
            if ptm is not None:
                mod_dict[pep_pos] = ptm

        mod_type = build_mod_type_string(pep_seq, mod_dict)
        mod_short = build_mod_short_name(modifiable_ks, combo, pos_start)

        peptidoforms.append((mod_short, mod_type))

    return peptidoforms


def build_mod_type_string(pep_seq: str, mod_dict: Dict[int, str]) -> str:
    """
    Build the mod_type string for EpiProfile.

    Rules:
    - Position 0 (N-terminal): always '0,pr;'
    - K at position 0: add '0,acD3;' (unmod) or '0,ac;' (acetylated) etc.
    - Internal K (not C-terminal): 'pos,acD3;' (unmod) or 'pos,ac;' etc.
    - C-terminal residue: nothing
    """
    parts = ['0,pr;']  # N-terminal propionylation always first

    last_idx = len(pep_seq) - 1

    for i, aa in enumerate(pep_seq):
        if aa == 'K' and i < last_idx:  # Modifiable K (not C-terminal)
            if i in mod_dict:
                ptm = mod_dict[i]
                parts.append(f'{i},{ptm};')
            else:
                parts.append(f'{i},acD3;')

    return ''.join(parts)


def build_mod_short_name(modifiable_ks: List[int], combo: tuple,
                         pos_start: int) -> str:
    """Build the mod_short name like 'unmod', 'K13ac', 'K5me2K13ac'."""
    active_mods = []
    for pep_pos, ptm in zip(modifiable_ks, combo):
        if ptm is not None:
            protein_pos = pos_start + pep_pos
            active_mods.append(f'K{protein_pos}{ptm}')

    if not active_mods:
        return 'unmod'

    return ''.join(active_mods)


def determine_charges(pep_len: int, charges_str: str) -> str:
    """Return the charge state vector string from manifest."""
    # Use the charges from manifest as-is
    return charges_str


def generate_matlab_module(module: dict, peptidoforms: List[Tuple[str, str]]) -> str:
    """Generate the complete MATLAB module file content."""
    func_name = module['file_name'].replace('.m', '')
    pep_seq = module['pep_seq']
    charges = module['charges'].split()
    n_charges = len(charges)
    charge_vec = ' '.join(charges)
    n_peps = len(peptidoforms)

    # Build mod_short cell array
    mod_short_lines = []
    for i, (ms, mt) in enumerate(peptidoforms):
        prefix = "His.mod_short = {" if i == 0 else "    "
        suffix = "};" if i == n_peps - 1 else ""
        mod_short_lines.append(f"{prefix}'{ms}'{suffix}")

    # Build mod_type cell array
    mod_type_lines = []
    for i, (ms, mt) in enumerate(peptidoforms):
        prefix = "His.mod_type = {" if i == 0 else "    "
        suffix = "};" if i == n_peps - 1 else ""
        mod_type_lines.append(f"{prefix}'{mt}'{suffix}")

    # Build rt_ref
    rt_ref_lines = []
    for i in range(n_peps):
        prefix = "His.rt_ref = [" if i == 0 else "    "
        suffix = "];" if i == n_peps - 1 else ""
        rt_ref_lines.append(f"{prefix}0{suffix}")

    # Build the calculate_layout section for modified peptidoforms
    calc_mods_section = ""
    if n_peps > 1:
        # Group by mod name for comments
        mod_comments = []
        for i, (ms, mt) in enumerate(peptidoforms):
            if i == 0:
                continue  # skip unmod
            mod_comments.append(f"% {ms}")

        calc_mods_section = f"""
% {peptidoforms[1][0] if n_peps > 1 else ''}
for hno=2:{n_peps}
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;"""
    else:
        calc_mods_section = f"""
for hno=2:1
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;"""

    # Build relocate section
    relocate_body = ""
    relocate2_body = ""

    if n_peps > 1:
        reloc_parts = []
        reloc2_parts = []
        for i in range(1, n_peps):
            ms = peptidoforms[i][0]
            reloc_parts.append(f"""
% {ms}
hno = {i+1};
t1 = His.rt_ref(1)-delta;
t2 = His.rt_ref(1)+30;
[rts{i+1},top1_rt{i+1}] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts{i+1})
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt{i+1};
end;""")
            reloc2_parts.append(f"""
% {ms}
hno = {i+1};
t1 = His.rt_ref(1)-delta;
t2 = His.rt_ref(1)+30;
[rts{i+1},top1_rt{i+1}] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts{i+1})
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt{i+1};
end;""")

        relocate_body = ''.join(reloc_parts)
        relocate2_body = ''.join(reloc2_parts)

    # Assemble the complete file
    mod_short_str = '\n'.join(mod_short_lines)
    mod_type_str = '\n'.join(mod_type_lines)
    rt_ref_str = '\n'.join(rt_ref_lines)

    matlab = f"""function {func_name}(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = '{func_name}';
fprintf(1,'%s..',out_filename(2:end));
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone();

% calculate
unitdiff = 1.0032;
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% draw
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Get PSM
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone()
%%

His.pep_seq = '{pep_seq}';
{mod_short_str}
{mod_type_str}

His.pep_ch = repmat([{charge_vec}],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
{rt_ref_str}
His.display = ones(length(His.mod_type),1);

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% unmod
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        [His.rt_ref(1),special.ndebug] = check_ref(special.raw_path,[His.pep_seq,His.mod_type{{1}}],His.rt_ref(1),special.ndebug);
    else
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{{1}}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end;
        hno = 1;% unmod
        [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,1,t1,t2,nhmass);%#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end;
    end;
end;

hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% calibrate the rt_ref
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end;
end;
{calc_mods_section}

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;
nsplit = 1;
{relocate_body}

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;
{relocate2_body}
"""
    return matlab


def main():
    print("=" * 70)
    print("EpiProfile-PLANTS Module Expander: H2A/H2B/H1 PTM Peptidoforms")
    print("=" * 70)

    # Parse Grau-Bove PTM reference
    print("\n[1/4] Parsing Grau-Bove PTM reference table...")
    graubove_ptms = parse_graubove_ptms(GRAUBOVE_PATH)
    print(f"  Found {len(graubove_ptms)} unique (histone, position) PTM entries")

    # Show summary
    for ht in ['H2A', 'H2AZ', 'H2B', 'H1']:
        entries = {k: v for k, v in graubove_ptms.items() if k[0] == ht}
        if entries:
            positions = sorted([k[1] for k in entries])
            print(f"  {ht}: K positions {positions}")
            for k, v in sorted(entries.items()):
                print(f"    K{k[1]}: {sorted(v)}")

    # Parse manifest
    print("\n[2/4] Parsing plant module manifest...")
    modules = parse_manifest(MANIFEST_PATH)
    print(f"  Found {len(modules)} AT H2A/H2B/H1 modules")

    # Process each module
    print("\n[3/4] Generating expanded peptidoforms...")
    stats = {'expanded': 0, 'unchanged': 0, 'total_peptidoforms': 0}
    module_results = []

    for mod in modules:
        ptm_map = determine_ptms_for_module(mod, graubove_ptms)
        peptidoforms = generate_peptidoforms(
            mod['pep_seq'], ptm_map, mod['pos_start']
        )

        n_pf = len(peptidoforms)
        stats['total_peptidoforms'] += n_pf

        if n_pf > 1:
            stats['expanded'] += 1
            mod_names = [p[0] for p in peptidoforms]
            print(f"  {mod['file_name']:30s} {mod['pep_seq']:30s} -> {n_pf} peptidoforms: {mod_names}")
        else:
            stats['unchanged'] += 1

        module_results.append((mod, peptidoforms))

    print(f"\n  Summary: {stats['expanded']} expanded, {stats['unchanged']} unchanged, "
          f"{stats['total_peptidoforms']} total peptidoforms")

    # Write MATLAB modules
    print(f"\n[4/4] Writing MATLAB modules to {TIER3_DIR}...")
    for mod, peptidoforms in module_results:
        matlab_code = generate_matlab_module(mod, peptidoforms)

        output_path = TIER3_DIR / mod['file_name']
        with open(output_path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(matlab_code)

    print(f"  Written {len(module_results)} modules")

    # Verification summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    print(f"Total modules processed: {len(module_results)}")
    print(f"Modules with multiple peptidoforms: {stats['expanded']}")
    print(f"Modules with only unmod: {stats['unchanged']}")
    print(f"Total peptidoforms across all modules: {stats['total_peptidoforms']}")

    # Show detailed breakdown for expanded modules
    print("\nExpanded modules detail:")
    for mod, peptidoforms in module_results:
        if len(peptidoforms) > 1:
            print(f"\n  {mod['file_name']} ({mod['histone_name']}):")
            print(f"    Peptide: {mod['pep_seq']} (pos {mod['pos_start']}-{mod['pos_end']})")
            for ms, mt in peptidoforms:
                print(f"    {ms:30s} -> {mt}")


if __name__ == '__main__':
    main()
